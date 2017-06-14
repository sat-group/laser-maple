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
#include <queue>

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


// returns true if we should re-enqueue the state (there are more branches to consider)
bool generateSolverFromState(Solver& s2, Solver::SolverState* state){
	// add learnts from the state
	if(s2.lsr_verb)
		printf("Generating\n");
	s2.printState(state);
	for(unsigned i = 0; i < state->mLearnts->size(); i++){
		vec<Lit>* clause = new vec<Lit>;
		std::vector<Lit> v = state->mLearnts->at(i);
		for (unsigned j = 0; j < v.size(); j++)
			clause->push(v[j]);
		s2.addClause(*clause);
		s2.prev_learnts.push_back(clause);
	}
	// add the found units (search will immediately propagate)
	for(unsigned i = 0; i < state->mFoundUnits->size(); i++){
		vec<Lit> clause;
		clause.push(state->mFoundUnits->at(i));
		if(s2.lsr_verb)
			printf("Adding unit: %s%d\n", sign(state->mFoundUnits->at(i))?"-":"", var(state->mFoundUnits->at(i)));
		s2.addClause(clause);
		s2.found_units.push_back(state->mFoundUnits->at(i));
	}

	// if lit_Undef is the next var to add, it represents restarts, also the last lit added to mNext (ret false)
	if(state->mNextIndex < state->mNext->size() && state->mNext->at(state->mNextIndex) == lit_Undef){
		return false;
	}


	s2.curr_replay_trail_index = 0; //s2.nAssigns();
	// reconstruct the trail -- the branching heuristic will simply follow this ordering
	for(unsigned i = 0; i < state->mTrail->size(); i++){
		s2.replay_trail.push_back(state->mTrail->at(i));
	}

	if(!state->mNext->empty()){
		// PRE -- mNextIndex < mNext.size()
		s2.replay_trail.push_back(state->mNext->at(state->mNextIndex));
		state->mNextIndex = state->mNextIndex + 1;
		if(state->mNextIndex < state->mNext->size())
			return true;
		else
			return false;
	}
	return false;
}

// copies the initial problem from S to s2, and sets the decision vars
void initializeSolver(SimpSolver& S, Solver& s2, vec<Var>& vars){
	s2.set_decision_vars = true;
	s2.lsr_verb = S.lsr_verb;
	s2.growDecision(S.nVars());
	for(int i = 0; i < S.nVars(); i++)
		s2.decision[i] = 0;
	for(int i =	 0; i < vars.size(); i++)
		s2.decision[vars[i]] = 1;

	while(s2.nVars() < S.nVars()){
		s2.newVar();
	}

	for(int i = 0; i < S.nClauses(); i++){
		Clause* cl = S.get_clause(i);
		vec<Lit> clause;
		for (int j = 0; j < cl->size(); j++){
			clause.push((*cl)[j]);
		}
		s2.addClause(clause);
	}
}

bool processPotentialBD(SimpSolver& S, vec<Var>& vars){
	if(S.lsr_verb){
		printf("Processing: ");
		for(int i = 0; i < vars.size(); i++)
			printf("%d ", vars[i] + 1);
		printf("\n");
	}
	printf("Processing: ");
	for(int i = 0; i < vars.size(); i++)
		printf("%d ", vars[i] + 1);
	printf("\n");

	std::queue<Solver::SolverState*> Q;
	Solver::SolverState* initial_state = new Solver::SolverState;
    std::vector<Solver::SolverState*> seen_states;

    for(int i = 0; i < S.trivial_units.size(); i++){
    	initial_state->mFoundUnits->push_back(S.trivial_units[i]);
    }

    int largest_trail = 0;

	Q.push(initial_state);
	while(!Q.empty()){
		Solver s2;
		initializeSolver(S, s2, vars);
		if(s2.nLearnts() > 0){
			printf("Should not hit\n");
			exit(1);
		}
		//s2.nuke();
		//printf("Nlearnts: %d \n", s2.nLearnts());
		if(s2.nLearnts() > 0){
			for(int i = 0; i < s2.nLearnts(); i++){
				printf("%d", s2.get_learnt(i));
			}
			exit(1);
		}
		Solver::SolverState* state = Q.front();

		int tsize = state->mLearnts->size();
		if(tsize > largest_trail)
			largest_trail = tsize;

		bool retainState = generateSolverFromState(s2, state);

		vec<Lit> dummy;
		/*
		printf("Solving with trail: (");
		for(unsigned i = 0; i < state->mFoundUnits->size(); i++)
			printf("%s%d ", sign(state->mFoundUnits->at(i)) ? "-" : "", var(state->mFoundUnits->at(i)));
		printf(") ");
		for(unsigned i = 0; i < s2.replay_trail.size(); i++)
			printf("%s%d ", sign(s2.replay_trail[i]) ? "-" : "", var(s2.replay_trail[i]));
		printf("\n");
		 */
		lbool res = s2.solveLimited(dummy);

		// copy mAllDecisions from previous state to new state, and add the most recent decision var
		for(unsigned i = 0; i < state->mAllDecisions->size(); i++){
			//printf("test %d\n", s2.final_state->mAllDecisions->size());
			s2.final_state->mAllDecisions->push_back(state->mAllDecisions->at(i));
		}
		if(!s2.skipped_last_lit && s2.replay_trail.size() > 0)
			s2.final_state->mAllDecisions->push_back(s2.replay_trail[s2.replay_trail.size() - 1]);
		//if(s2.skipped_last_lit)
		//	printf("CAUGHT %d\n", s2.replay_trail[s2.replay_trail.size() - 1]);
		if(res != l_Undef){
			printf("SUCCESS!\n");
			Solver::SolverState s = *(s2.final_state);
			printf("Decisions:\n");
			for(unsigned i = 0; i < s.mAllDecisions->size(); i++){
				printf("%s%d\n", sign(s.mAllDecisions->at(i)) ? "-" : "", var(s.mAllDecisions->at(i))+1);
			}
			printf("\n");
			if(S.replay_bd_out_file){
				FILE* f = fopen(S.replay_bd_out_file, "wr");
				for(unsigned i = 0; i < s.mAllDecisions->size(); i++){
					fprintf(f, "%s%d\n", sign(s.mAllDecisions->at(i)) ? "-" : "", var(s.mAllDecisions->at(i))+1);
				}
				fclose(f);
			}
			exit(1);
		}
		// printf("Res (0T/1F/2U): %d\n", res);

		//printf("Next Vars: ");
		//for(int i = 0; i < s2.next_vars.size(); i++){
		//	 printf("%d ", s2.next_vars[i]);
		//}
		//printf("\n");
		if(!Solver::existingState(s2.final_state, seen_states)){
			seen_states.push_back(s2.final_state);
			Q.push(s2.final_state);
		}
		else{
			if(s2.lsr_verb)
				printf("Not exploring:\n");
			s2.printState(s2.final_state);
			delete s2.final_state;
		}
		// todo figure out right way to clear unit variables (currently initting solver in loop).

		// we've exhausted all branches for this state
		if(!retainState)
			Q.pop();
	}
	if(S.lsr_verb){
		printf("Len Seen %d\n", seen_states.size());
		printf("%d\n", largest_trail);
	}
	for(int i = 0; i < seen_states.size(); i++){
		delete seen_states[i];
	}
	delete initial_state;
	return false;
}

bool expand(SimpSolver& S, vec<Var>& vars, int level)
{
  vec<Var> tmp;
  vars.copyTo(tmp);

  if (level == 1 && processPotentialBD(S, vars))
    return true;

  if (level <= 1)
    return false;

  for (int i = vars[vars.size() - 1] + 1; i < S.nVars(); i++) {
	  if(S.value(i) == l_Undef){
		tmp.push(i);
		expand(S, tmp, level - 1);
		tmp.pop();
	  }
  }
  return false;
}


bool findMinLS(SimpSolver& S, int minSize)
{
	vec<vec<Var> > expand_list;
	for (int i = 0; i < S.nVars(); i++) {
		if (S.value(i) == l_Undef) {
		  expand_list.push();
		  expand_list[expand_list.size() - 1].push(i);
		}
	}
	for (int i = 0; i < expand_list.size(); i++) {
		expand(S, expand_list[i], minSize);
	}
	return false;
}


bool findMinLS(SimpSolver& S, vec< vec<Var>* >& var_sets)
{
	for(int i = 0; i < var_sets.size(); i++){
		vec<Var>* vars = var_sets[i];
		processPotentialBD(S, *vars);
	}
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
        BoolOption   lsr_num("LASER","lsr-num","Number of LSR backdoor variables.\n",false);

        // minimum LSR options
        BoolOption   lsr_verb("LASER","lsr-verb", "Debugging for MinLS tool.\n",false);
        IntOption    min_lsr_size("LASER", "min-lsr-size", "Minimum size LSR set to try.\n", 1, IntRange(1, INT32_MAX));
        IntOption    max_lsr_size("LASER", "max-lsr-size", "Maximum size LSR set to try.\n", 5, IntRange(1, INT32_MAX));
        BoolOption   replay_bd("LASER","replay-bd", "Replay the decisions that witness the backdoor.\n",false);
        StringOption replay_bd_file ("MAIN", "replay-bd-file", "Infile of decisions that witness the backdoor (one-based).");
        StringOption var_sets_file ("MAIN", "var-sets-file", "List of potential backdoor sets to try (dimacs format without header).");
        StringOption replay_bd_out_file ("MAIN", "replay-bd-out", "Outfile of decisions that witness the backdoor (one-based).");



        parseOptions(argc, argv, true);
        
        SimpSolver  S;
        double      initial_time = cpuTime();

        if (!pre) S.eliminate(true);

        S.verbosity = verb;
        S.lsr_verb = lsr_verb;
        S.replay_bd_out_file = (const char*) replay_bd_out_file;
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
        //signal(SIGINT, SIGINT_interrupt);
        //signal(SIGXCPU,SIGINT_interrupt);



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

        vec< vec<Var>* > var_sets;
        if (var_sets_file) {
			const char* file_name = var_sets_file;
			FILE* var_sets_file = fopen (file_name, "r");
			if (var_sets_file == NULL)
				printf("ERROR! Could not open file: %s\n", file_name), exit(1);
			int i = 0;
			vec<Var>* vars = new vec<Var>;
			while (fscanf(var_sets_file, "%d", &i) == 1) {
				if(i == 0){
					var_sets.push(vars);
					vars = new vec<Var>;
				}
				else{
					Var v = abs(i) - 1;
					vars->push(v);
				}
			}
			//printf("VSIZE:%d\n", var_sets.size());
			fclose(var_sets_file);
		}



        for( int i = 0; i < dummy.size(); i++) {
            printf("%s%d\n", sign(dummy[i]) ? "-" : "", var(dummy[i]));
        }

        S.setLSR(lsr_num);
        S.setFilename(lsr_file);

        if(replay_bd){
        	std::vector<Lit> replay_lits;
        	const char* replay_bd_file_name = replay_bd_file;
			FILE* replay_file = fopen(replay_bd_file_name, "r");
			if (replay_file == NULL)
				printf("ERROR! Could not open replay bd file: %s\n", replay_bd_file_name), exit(1);
			int i = 0;
			int v;
			while (fscanf(replay_file, "%d", &i) == 1) {
				v = abs(i) - 1;
				Lit l = (i > 0) ? mkLit(v) : ~mkLit(v);
				replay_lits.push_back(l);
			}
			fclose(replay_file);
			// create replay_trail
			printf("Full trail: ");
			for(i = 0; i < replay_lits.size(); i++){
				S.replay_trail.push_back(replay_lits[i]);
				printf("%s%d ", sign(replay_lits[i])?"-":"", var(replay_lits[i]));
			}
			printf("\n");
	        lbool ret = S.solveLimited(dummy);
	        printf("Result: %d\n", ret);
        	exit(1);
        }

        if(var_sets.size() != 0)
        	findMinLS(S, var_sets);
        else
			for(int level = min_lsr_size; level <= max_lsr_size; level++)
				findMinLS(S, level);
        exit(1);


        lbool ret = S.solveLimited(dummy);
        
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
