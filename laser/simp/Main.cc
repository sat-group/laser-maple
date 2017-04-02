/*****************************************************************************************[Main.cc]
Copyright (c) 2003-2006, Niklas Een, Niklas Sorensson
Copyright (c) 2007,      Niklas Sorensson

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and
associated documentation files (the "Software"), to deal in the Software without restriction,
including without limitation the rights to use, copy, modify, merge, publish, distribute,
sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or
substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT
NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM,
DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT
OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
**************************************************************************************************/

#include <errno.h>

#include <random>  // EDXXX for rand
#include <vector>  // EDXXX
#include <algorithm> //EDXXX

#include <iostream>
#include <time.h> //EDXXX
#include <signal.h>
#include <zlib.h>
#include <sys/resource.h>

#include "utils/System.h"
#include "utils/ParseUtils.h"
#include "utils/Options.h"
#include "core/Dimacs.h"
#include "simp/SimpSolver.h"

using namespace Minisat;

//=================================================================================================


void printStats(Solver& solver)
{
    double cpu_time = cpuTime();
    double mem_used = memUsedPeak();
    printf("restarts              : %"PRIu64"\n", solver.starts);
    printf("conflicts             : %-12"PRIu64"   (%.0f /sec)\n", solver.conflicts   , solver.conflicts   /cpu_time);
    printf("decisions             : %-12"PRIu64"   (%4.2f %% random) (%.0f /sec)\n", solver.decisions, (float)solver.rnd_decisions*100 / (float)solver.decisions, solver.decisions   /cpu_time);
    printf("propagations          : %-12"PRIu64"   (%.0f /sec)\n", solver.propagations, solver.propagations/cpu_time);
    printf("conflict literals     : %-12"PRIu64"   (%4.2f %% deleted)\n", solver.tot_literals, (solver.max_literals - solver.tot_literals)*100 / (double)solver.max_literals);
    long double total_actual_rewards = 0;
    long double total_actual_count = 0;
    for (int i = 0; i < solver.nVars(); i++) {
        total_actual_rewards += solver.total_actual_rewards[i];
        total_actual_count += solver.total_actual_count[i];
    }
    printf("actual reward         : %Lf\n", total_actual_rewards / total_actual_count);
    if (mem_used != 0) printf("Memory used           : %.2f MB\n", mem_used);
    printf("CPU time              : %g s\n", cpu_time);
}


static Solver* solver;
// Terminate by notifying the solver and back out gracefully. This is mainly to have a test-case
// for this feature of the Solver as it may take longer than an immediate call to '_exit()'.
static void SIGINT_interrupt(int signum) { solver->interrupt(); }

// Note that '_exit()' rather than 'exit()' has to be used. The reason is that 'exit()' calls
// destructors and may cause deadlocks if a malloc/free function happens to be running (these
// functions are guarded by locks for multithreaded use).
static void SIGINT_exit(int signum) {
    printf("\n"); printf("*** INTERRUPTED ***\n");
    if (solver->verbosity > 0){
        printStats(*solver);
        printf("\n"); printf("*** INTERRUPTED ***\n"); }
    _exit(1); }



void weighted_choice(vec<int>& num_occs, vec<int>& chosen, int target_size, int initial_sum){
	std::random_device rd;
	std::mt19937 gen(rd());
	int curr_sum = initial_sum;
	vec<bool> picked;
	bool flag;
	picked.growTo(num_occs.size(), false);
	for(int iter = 0; iter < target_size; iter++){
		if(curr_sum == 0){
			chosen.clear();
			return;
			//printf("ERROR in weighted choice -- degenerate curr_sum\n");
			//exit(1);
		}
		std::uniform_int_distribution<unsigned long long> dis(1, curr_sum);
		int val = dis(gen);
		//printf("val: %d %d\n", val, initial_sum);
		int upto = 0;
		flag = false;
		for(int i = 0; i < num_occs.size(); i++){
			if(picked[i])
				continue;
			if(upto + num_occs[i] >= val){
				//hit, add to chosen, makes sure it can't be picked again
				picked[i] = true;
				curr_sum -= num_occs[i];
				chosen.push(i);
				flag = true;
				break;
			}
			upto += num_occs[i];
		}
		if(!flag){
			printf("ERROR in weighted choice\n");
			exit(1);
		}

	}
	if(target_size != chosen.size()){
		printf("ERROR, wrong chosen size!\n");
	}
	//printf("Chosen(%d, %d): ", chosen.size(), target_size);

}


double gini(vec<double>& vals){
	// compute gini coefficient of normalized picks
	for(int i = 0; i < vals.size() - 1; i++){
		for(int j = i + 1; j < vals.size(); j++){
			if(vals[i] > vals[j]){
				double temp = vals[i];
				vals[i] = vals[j];
				vals[j] = temp;
			}
		}
	}
	double height = 0;
	double area = 0;
	double fair_area = 0;
	for(int i = 0; i < vals.size(); i++){
		height += vals[i];
		area += height - (vals[i] / 2);
	}
	fair_area = height * (float(vals.size()) / 2);
	if(fair_area != 0)
		return (fair_area - area) / fair_area;
	else
		return -1;
}

//=================================================================================================
// Main:

int main(int argc, char** argv)
{
    try {
        setUsageHelp("USAGE: %s [options] <input-file> <result-output-file>\n\n  where input may be either in plain or gzipped DIMACS.\n");
        // printf("This is MiniSat 2.0 beta\n");
        
#if defined(__linux__)
        fpu_control_t oldcw, newcw;
        _FPU_GETCW(oldcw); newcw = (oldcw & ~_FPU_EXTENDED) | _FPU_DOUBLE; _FPU_SETCW(newcw);
        //printf("WARNING: for repeatability, setting FPU to use double precision\n");
#endif
        // Extra options:
        //
        IntOption    verb   ("MAIN", "verb",   "Verbosity level (0=silent, 1=some, 2=more).", 1, IntRange(-1, 2));
        BoolOption   pre    ("MAIN", "pre",    "Completely turn on/off any preprocessing.", true);
        StringOption dimacs ("MAIN", "dimacs", "If given, stop after preprocessing and write the result to this file.");
        StringOption assumptions ("MAIN", "assumptions", "If given, use the assumptions in the file.");
        IntOption    cpu_lim("MAIN", "cpu-lim","Limit on CPU time allowed in seconds.\n", INT32_MAX, IntRange(0, INT32_MAX));
        IntOption    mem_lim("MAIN", "mem-lim","Limit on memory usage in megabytes.\n", INT32_MAX, IntRange(0, INT32_MAX));

        StringOption decision_vars ("MAIN", "decision-vars", "Only branch on the listed vars (zero-based).");

        // LASER options:
        StringOption lsr_file("LASER","lsr-out","Write LSR backdoor to a file (zero-based).\n");
        StringOption all_decisions_file("LASER","all-dec-out","Write all unique decision vars to a file (same as LS paper) (zero-based).\n");

        BoolOption   lsr_num("LASER","lsr-num","Number of LSR backdoor variables.\n",false);

        // certificate generation
        StringOption lsr_file_in("LASER","lsr-in","Used to create a certificate for an lsr, generate with -lsr-out (zero-based).\n");
        StringOption lsr_certificate_clauses_file_out("LASER","lsr-cert-cls-out","File to output the clauses sequence witnessing a backdoor (one-based).\n");

        // certificate validation
        StringOption lsr_certificate_clauses_file_in("LASER","lsr-cert-cls-in","File containing the clause sequence witnessing a backdoor (one-based).\n");

        // structure logging
        StringOption cmty_file("LASER","cmty-file", "var+cmty pairs, where vars should be 0 based");
        StringOption backbone_file("LASER","backbone-file", "backbone literals (one-based)");


        parseOptions(argc, argv, true);
        
        SimpSolver  S;
        double      initial_time = cpuTime();

        if (!pre) S.eliminate(true);

        S.verbosity = verb;
        solver = &S;
        // Use signal handlers that forcibly quit until the solver will be able to respond to
        // interrupts:
        signal(SIGINT, SIGINT_exit);
        signal(SIGXCPU,SIGINT_exit);

        // Set limit on CPU-time:
        if (cpu_lim != INT32_MAX){
            rlimit rl;
            getrlimit(RLIMIT_CPU, &rl);
            if (rl.rlim_max == RLIM_INFINITY || (rlim_t)cpu_lim < rl.rlim_max){
                rl.rlim_cur = cpu_lim;
                if (setrlimit(RLIMIT_CPU, &rl) == -1)
                    printf("WARNING! Could not set resource limit: CPU-time.\n");
            } }

        // Set limit on virtual memory:
        if (mem_lim != INT32_MAX){
            rlim_t new_mem_lim = (rlim_t)mem_lim * 1024*1024;
            rlimit rl;
            getrlimit(RLIMIT_AS, &rl);
            if (rl.rlim_max == RLIM_INFINITY || new_mem_lim < rl.rlim_max){
                rl.rlim_cur = new_mem_lim;
                if (setrlimit(RLIMIT_AS, &rl) == -1)
                    printf("WARNING! Could not set resource limit: Virtual memory.\n");
            } }
        
        if (argc == 1)
            printf("Reading from standard input... Use '--help' for help.\n");

        gzFile in = (argc == 1) ? gzdopen(0, "rb") : gzopen(argv[1], "rb");
        if (in == NULL)
            printf("ERROR! Could not open file: %s\n", argc == 1 ? "<stdin>" : argv[1]), exit(1);
        


        if (S.verbosity > 0){
            printf("============================[ Problem Statistics ]=============================\n");
            printf("|                                                                             |\n"); }

        S.setDecisionVarsList(decision_vars);

        parse_DIMACS(in, S);


        // maps vars to cmtys
		if(cmty_file){
			S.structure_logging = true;
			S.cmty_logging = true;
			S.var_cmty.growTo(S.nVars(), -1);
			FILE* cfile = fopen(cmty_file, "r");
			if (cfile == NULL)
				fprintf(stderr, "could not open file %s\n", (const char*) cfile), exit(1);
			int v;
			int c;
			//Zero-based
			int largest_cmty_index = -1;
			while (fscanf(cfile, "%d %d\n", &v, &c) == 2){
				S.var_cmty[v] = c;
				if(c > largest_cmty_index){
					largest_cmty_index = c;
					S.cmty_picks.growTo(largest_cmty_index + 1, 0);
					S.cmty_size.growTo(largest_cmty_index + 1, 0);
					S.cmty_clauses.growTo(largest_cmty_index + 1, 0);
				}
				S.cmty_size[c] = S.cmty_size[c] + 1;
			}
			fclose(cfile);
		}

        // lsr verification files
        if(backbone_file){
        	S.structure_logging = true;
			S.backbone_logging = true;
        	S.backbone.growTo(S.nVars(), 0);
        	S.backbone_prev_polarity.growTo(S.nVars(), 2);
			const char* file_name = backbone_file;
			FILE* backbone_in = fopen (file_name, "r");
			if (backbone_in == NULL)
				printf("ERROR! Could not open file: %s\n", file_name), exit(1);
			int i = 0;
			while (fscanf(backbone_in, "%d", &i) == 1) {
				Var v = abs(i) - 1;
				if(i > 0)
					S.backbone[v] = 1;
				else
					S.backbone[v] = -1;
			}
			fclose(backbone_in);
        }

        // lsr verification files
        if(lsr_file_in){
        	S.lsr_in.growTo(S.nVars(), 0);
			const char* file_name = lsr_file_in;
			FILE* lsr_in = fopen (file_name, "r");
			if (lsr_in == NULL)
				printf("ERROR! Could not open file: %s\n", file_name), exit(1);
			int i = 0;
			while (fscanf(lsr_in, "%d", &i) == 1) {
				S.lsr_in[i] = true;
			}
			fclose(lsr_in);
        }

		if(lsr_file_in && lsr_certificate_clauses_file_out){
			S.generate_certificate = true;
			const char* file_name3 = lsr_certificate_clauses_file_out;
			S.certificate_clauses_out = fopen(file_name3, "wb");
		}
        
		if(lsr_certificate_clauses_file_in){
			int i;
			S.verification_mode = true;
			const char* file_name2 = lsr_certificate_clauses_file_in;
			FILE* cls_in = fopen (file_name2, "r");
			if (cls_in == NULL)
				printf("ERROR! Could not open file: %s\n", file_name2), exit(1);
			vec<Lit>* c = new vec<Lit>;
			while (fscanf(cls_in, "%d", &i) == 1) {
				if(i == 0){
					S.cert.push(c);
					c = new vec<Lit>;
				}
				else{
					Var v = abs(i) - 1;
					Lit l = i > 0 ? mkLit(v) : ~mkLit(v);
					c->push(l);
				}
			}
			S.curr_cert_index = 0;
			fclose(cls_in);

		}


        if(S.nVars() == 0 && S.nClauses() == 0){ // EDXXX
        	printf("TRIVIAL\n");
        	exit(0);
        }
        gzclose(in);

        FILE* res = (argc >= 3) ? fopen(argv[2], "wb") : NULL;

        if (S.verbosity > 0){
            printf("|  Number of variables:  %12d                                         |\n", S.nVars());
            printf("|  Number of clauses:    %12d                                         |\n", S.nClauses()); }
        
        double parsed_time = cpuTime();
        if (S.verbosity > 0)
            printf("|  Parse time:           %12.2f s                                       |\n", parsed_time - initial_time);

        // Change to signal-handlers that will only notify the solver and allow it to terminate
        // voluntarily:
        signal(SIGINT, SIGINT_interrupt);
        signal(SIGXCPU,SIGINT_interrupt);



        S.eliminate(true);


        double simplified_time = cpuTime();
        if (S.verbosity > 0){
            printf("|  Simplification time:  %12.2f s                                       |\n", simplified_time - parsed_time);
            printf("|                                                                             |\n"); }

        if (!S.okay()){
            if (res != NULL) fprintf(res, "UNSAT\n"), fclose(res);
            if (S.verbosity > 0){
                printf("===============================================================================\n");
                printf("Solved by simplification\n");
                printStats(S);
                printf("\n"); }
            printf("UNSATISFIABLE\n");
            exit(20);
        }

        if (dimacs){
            if (S.verbosity > 0)
                printf("==============================[ Writing DIMACS ]===============================\n");
            S.toDimacs((const char*)dimacs);
            if (S.verbosity > 0)
                printStats(S);
            exit(0);
        }

        vec<Lit> dummy;
        if (assumptions) {
            const char* file_name = assumptions;
            FILE* assertion_file = fopen (file_name, "r");
            if (assertion_file == NULL)
                printf("ERROR! Could not open file: %s\n", file_name), exit(1);
            int i = 0;
            while (fscanf(assertion_file, "%d", &i) == 1) {
                Var v = abs(i) - 1;
                Lit l = i > 0 ? mkLit(v) : ~mkLit(v);
                dummy.push(l);
            }
            fclose(assertion_file);
        }


        for( int i = 0; i < dummy.size(); i++) {
            printf("%s%d\n", sign(dummy[i]) ? "-" : "", var(dummy[i]));
        }

        S.setLSR(lsr_num);
        S.setFilename(lsr_file);
        S.setAllDecisionsFilename(all_decisions_file);

        lbool ret = S.solveLimited(dummy);
        
        if(S.structure_logging){
        	if(S.cmty_logging){
        		// compute normalized pick values for each cmty
        		vec<double> normalized_picks;
        		vec<double> normalized_cmty_clauses;
        		for(int i = 0; i < S.cmty_size.size(); i++){
        			if(S.cmty_size[i] != 0){
        				normalized_picks.push((float(S.cmty_picks[i])) / S.cmty_size[i]);
        			}
        		}
        		for(int i = 0; i < S.cmty_size.size(); i++){
					if(S.cmty_size[i] != 0){
						normalized_cmty_clauses.push((float(S.cmty_clauses[i])) / S.cmty_size[i]);
					}
        		}
        		double gini_normalized_picks = gini(normalized_picks);
        		if(gini_normalized_picks >= 0)
        			printf("GiniNormalizedPicks %f\n", gini_normalized_picks);
        		else
        			printf("Gini failed\n");

        		double gini_normalized_clauses = gini(normalized_cmty_clauses);
				if(normalized_cmty_clauses >= 0)
					printf("GiniNormalizedClauses %f\n", gini_normalized_clauses);
				else
					printf("Gini failed\n");

        	}
        	if(S.backbone_logging){
        		int backbone_size = 0;
        		for(int i = 0; i < S.backbone.size(); i++){
        			if(S.backbone[i] != 0)
        				backbone_size++;
        		}
        		printf("NormalizedBackboneFlips %f\n", S.num_backbone_flips / float(backbone_size));
        		printf("NormalizedBackboneSubsumedClauses %f\n", S.num_backbone_subsumed_clauses / float(backbone_size));

        	}
        }
        if(S.compute_avg_clause_lsr){
        	printf("AvgClauseLSR %f\n", S.total_clause_lsr_weight / S.all_learnts);
        }

        if (S.verbosity > 0){
            printStats(S);
            printf("\n"); }

        if(S.verbosity > -1)
        	printf(ret == l_True ? "SATISFIABLE\n" : ret == l_False ? "UNSATISFIABLE\n" : "INDETERMINATE\n");
        if (res != NULL){
            if (ret == l_True){
                fprintf(res, "SAT\n");
                for (int i = 0; i < S.nVars(); i++)
                    if (S.model[i] != l_Undef)
                        fprintf(res, "%s%s%d", (i==0)?"":" ", (S.model[i]==l_True)?"":"-", i+1);
                fprintf(res, " 0\n");
            }else if (ret == l_False) {
                fprintf(res, "UNSAT\n");
                for (int i = 0; i < S.conflict.size(); i++) {
                    // Reverse the signs to keep the same sign as the assertion file.
                    fprintf(res, "%s%d\n", sign(S.conflict[i]) ? "" : "-", var(S.conflict[i]) + 1);
                }
            } else
                fprintf(res, "INDET\n");
            fclose(res);
        }

#ifdef NDEBUG
        exit(ret == l_True ? 10 : ret == l_False ? 20 : 0);     // (faster than "return", which will invoke the destructor for 'Solver')
#else
        return (ret == l_True ? 10 : ret == l_False ? 20 : 0);
#endif
    } catch (OutOfMemoryException&){
        printf("===============================================================================\n");
        printf("INDETERMINATE\n");
        exit(0);
    }
}
