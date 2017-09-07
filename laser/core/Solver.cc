/***************************************************************************************[Solver.cc]
Copyright (c) 2003-2006, Niklas Een, Niklas Sorensson
Copyright (c) 2007-2010, Niklas Sorensson

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

#include <math.h>

#include "mtl/Sort.h"
#include "core/Solver.h"


using namespace Minisat;

//=================================================================================================
// Options:


static const char* _cat = "CORE";

#if BRANCHING_HEURISTIC == CHB || BRANCHING_HEURISTIC == LRB || BRANCHING_HEURISTIC == SGDB
static DoubleOption  opt_step_size         (_cat, "step-size",   "Initial step size",                             0.40,     DoubleRange(0, false, 1, false));
static DoubleOption  opt_step_size_dec     (_cat, "step-size-dec","Step size decrement",                          0.000001, DoubleRange(0, false, 1, false));
static DoubleOption  opt_min_step_size     (_cat, "min-step-size","Minimal step size",                            0.06,     DoubleRange(0, false, 1, false));
#endif
#if BRANCHING_HEURISTIC == VSIDS
static DoubleOption  opt_var_decay         (_cat, "var-decay",   "The variable activity decay factor",            0.95,     DoubleRange(0, false, 1, false));
#endif
#if BRANCHING_HEURISTIC == SGDB
static DoubleOption  opt_regularization    (_cat, "regularization", "L2 Regularization", 0.95, DoubleRange(0, true, 1, true));
#endif
#if ! LBD_BASED_CLAUSE_DELETION
static DoubleOption  opt_clause_decay      (_cat, "cla-decay",   "The clause activity decay factor",              0.999,    DoubleRange(0, false, 1, false));
#endif
static DoubleOption  opt_random_var_freq   (_cat, "rnd-freq",    "The frequency with which the decision heuristic tries to choose a random variable", 0, DoubleRange(0, true, 1, true));
static DoubleOption  opt_random_seed       (_cat, "rnd-seed",    "Used by the random variable selection",         91648253, DoubleRange(0, false, HUGE_VAL, false));
static IntOption     opt_ccmin_mode        (_cat, "ccmin-mode",  "Controls conflict clause minimization (0=none, 1=basic, 2=deep)", 2, IntRange(0, 2));
static IntOption     opt_phase_saving      (_cat, "phase-saving", "Controls the level of phase saving (0=none, 1=limited, 2=full)", 2, IntRange(0, 2));
static BoolOption    opt_rnd_init_act      (_cat, "rnd-init",    "Randomize the initial activity", false);
static BoolOption    opt_rnd_pol     (_cat, "rnd-pol",    "Randomize the polarity selection", false);
static BoolOption    opt_luby_restart      (_cat, "luby",        "Use the Luby restart sequence", true);
static IntOption     opt_restart_first     (_cat, "rfirst",      "The base restart interval", 100, IntRange(1, INT32_MAX));
static DoubleOption  opt_restart_inc       (_cat, "rinc",        "Restart interval increase factor", 2, DoubleRange(1, false, HUGE_VAL, false));
static DoubleOption  opt_garbage_frac      (_cat, "gc-frac",     "The fraction of wasted memory allowed before a garbage collection is triggered",  0.20, DoubleRange(0, false, HUGE_VAL, false));
#if BRANCHING_HEURISTIC == CHB
static DoubleOption  opt_reward_multiplier (_cat, "reward-multiplier", "Reward multiplier", 0.9, DoubleRange(0, true, 1, true));
#endif
static BoolOption    opt_always_restart      (_cat, "always-restart",        "Restart after every conflict.", false);
static BoolOption    opt_never_restart      (_cat, "never-restart",        "Restart never.", false);
static IntOption     opt_uniform_restart     (_cat, "uniform-restart",      "Uniform k-restarts", -1, IntRange(-1, INT32_MAX));

static BoolOption    opt_never_gc      (_cat, "never-gc",        "Never remove clauses.", false);

// lsr computation types
static BoolOption    opt_clause_and_conflict_side_lsr      (_cat, "conf-side-lsr",        "Dependencies of a clause are the clause itself and the dependents on the conflict side.", true);
static IntOption     opt_lsr_reduce     (_cat, "lsr-red",      "Use an LSR-centric clause deletion policy. 0) off; 1) sort based on LSR; 2) sort based on LSR-spread", 0, IntRange(0, 2));
static IntOption     opt_initial_max_learnts     (_cat, "init-max-learnts", "The initial maximum learnt clause database size", 2000, IntRange(1, INT32_MAX));
static IntOption     opt_learnt_db_bump     (_cat, "db-bump", "How much to bump the clause database size at each reduction", 500, IntRange(1, INT32_MAX));

static BoolOption    opt_embed_lsr      ("EMBED_LSR", "embed-lsr", "Given a decision vars file, add enough clauses to make them an LSR-backdoor.", false);
static IntOption     opt_embed_lsr_target_clause_size ("EMBED_LSR", "embed-lsr-cls-size", "How large to try to make each learnt clause.", 3, IntRange(1, INT32_MAX));
static IntOption     opt_embed_lsr_clause_type ("EMBED_LSR", "embed-lsr-cls-type", "When branching fails, how should the clause be learnt using the trail? (0 - random, 1 - beginning, 2 -end)", 1, IntRange(0, 2));


// structure logging
static StringOption   opt_average_clause_lsr_out("LASER","avg-clause-lsr-out","For each learnt, record its lsr size, compute the average. Dump to given file.");
static StringOption   opt_lsr_frequency("LASER","lsr-frequency-out","Record how many times each variable is a dependent of a clause.");


//lsr maxsat for SAT final deps
static StringOption opt_lsr_final_deps("LASER", "lsr-final-deps", "Smaller LSR for SAT case. Don't look at deps of propagated vars, just add it.");



//=================================================================================================
// Constructor/Destructor:


Solver::Solver() :

    // Parameters (user settable):
    //
    verbosity        (0)
#if BRANCHING_HEURISTIC == CHB || BRANCHING_HEURISTIC == LRB || BRANCHING_HEURISTIC == SGDB
  , step_size        (opt_step_size)
  , step_size_dec    (opt_step_size_dec)
  , min_step_size    (opt_min_step_size)
#endif
#if BRANCHING_HEURISTIC == VSIDS
  , var_decay        (opt_var_decay)
#endif
#if BRANCHING_HEURISTIC == SGDB
  , regularization   (opt_regularization)
#endif
#if ! LBD_BASED_CLAUSE_DELETION
  , clause_decay     (opt_clause_decay)
#endif
  , random_var_freq  (opt_random_var_freq)
  , random_seed      (opt_random_seed)
  , luby_restart     (opt_luby_restart)
  , ccmin_mode       (opt_ccmin_mode)
  , phase_saving     (opt_phase_saving)
  , rnd_pol          (opt_rnd_pol)
  , rnd_init_act     (opt_rnd_init_act)
  , garbage_frac     (opt_garbage_frac)
  , restart_first    (opt_restart_first)
  , restart_inc      (opt_restart_inc)

    // Parameters (the rest):
    //
  , learntsize_factor((double)1/(double)3), learntsize_inc(1.1)

  , uniform_restarts (opt_uniform_restart)
  , always_restart (opt_always_restart)
  , never_restart (opt_never_restart)
  , clause_and_conflict_side_lsr (opt_clause_and_conflict_side_lsr)
    // Parameters (experimental):
    //
  , learntsize_adjust_start_confl (100)
  , learntsize_adjust_inc         (1.5)

    // Statistics: (formerly in 'SolverStats')
    //
  , solves(0), starts(0), decisions(0), rnd_decisions(0), propagations(0), conflicts(0)
#if BRANCHING_HEURISTIC == SGDB
  , in_conflicts(0), sampleds(0)
#endif
  , dec_vars(0), clauses_literals(0), learnts_literals(0), max_literals(0), tot_literals(0)

  , lbd_calls(0)
#if BRANCHING_HEURISTIC == CHB
  , action(0)
  , reward_multiplier(opt_reward_multiplier)
#endif

  , ok                 (true)
#if ! LBD_BASED_CLAUSE_DELETION
  , cla_inc            (1)
#endif
#if BRANCHING_HEURISTIC == VSIDS
  , var_inc            (1)
#endif
  , watches            (WatcherDeleted(ca))
  , qhead              (0)
  , simpDB_assigns     (-1)
  , simpDB_props       (0)
  , order_heap         (VarOrderLt(activity))
  , progress_estimate  (0)
  , remove_satisfied   (true)

    // Resource constraints:
    //
  , conflict_budget    (-1)
  , propagation_budget (-1)
  , asynch_interrupt   (false)
  , generate_certificate (false)
  , verification_mode (false)
  , final_sat_trail_mode (false)
  , final_sat_trail_index (0)
  , restart_now(false)
  , restart_immediately(false)
  , just_restarted(false)
  , never_gc(opt_never_gc)
  , lsr_deletion_policy(opt_lsr_reduce)
  , initial_max_learnts(opt_initial_max_learnts)
  , db_bump (opt_learnt_db_bump)
  , structure_logging(false)
  , cmty_logging(false)

  , popsim_branching(false)
  , focused_branching_failure_limit(0)
  , focused_branching_remaining_failures(0)

  , all_learnts(0)
  , total_clause_lsr_weight(0)
  , num_backbone_flips(0)
  , num_backbone_subsumed_clauses(0)
  , backbone_logging(false)
  , embed_lsr(opt_embed_lsr)
  , embed_lsr_target_clause_size(opt_embed_lsr_target_clause_size)
  , embed_lsr_clause_type(opt_embed_lsr_clause_type)
{
  //strcpy(lsr_filename,"");
  lsr_filename = NULL;
  all_decisions_filename = NULL;
  //, compute_avg_clause_lsr(opt_average_clause_lsr)
  if(opt_average_clause_lsr_out){
	  const char* fname = opt_average_clause_lsr_out;
	  avg_clause_lsr_out = fopen(fname, "wb");
  }
  else
	  avg_clause_lsr_out = NULL;

  if(opt_lsr_final_deps){
  	  const char* fname = opt_lsr_final_deps;
  	lsr_final_deps_file = fname;
    }
    else
    	lsr_final_deps_file = NULL;
  if(opt_lsr_frequency){
	  record_lsr_frequency = true;
	  lsr_frequency_file = opt_lsr_frequency;
  }
  else{
	  record_lsr_frequency = false;
	  lsr_frequency_file = NULL;
  }
  lsr_num = false;
} 


Solver::~Solver()
{
	if(generate_certificate){
		fclose(certificate_clauses_out);

	}
}

double sigmoid(double z) {
    return 1.0 / (1.0 + exp(-z));
}


//=================================================================================================
// Minor methods:


// Creates a new SAT variable in the solver. If 'decision' is cleared, variable will not be
// used as a decision variable (NOTE! This has effects on the meaning of a SATISFIABLE result).
//
Var Solver::newVar(bool sign, bool dvar)
{
    int v = nVars();
    watches  .init(mkLit(v, false));
    watches  .init(mkLit(v, true ));
    assigns  .push(l_Undef);
    vardata  .push(mkVarData(CRef_Undef, 0));
    //activity .push(0);
    activity .push(rnd_init_act ? drand(random_seed) * 0.00001 : 0);
    seen     .push(0);
    polarity .push(sign);
    if(!set_decision_vars)
    	decision .push();
    trail    .capacity(v+1);
    lbd_seen.push(0);
    picked.push(0);
    conflicted.push(0);
    lsr_seen.push(0);
    all_decisions.push(0);
    unit_lsr.push(CRef_Undef);
#if ALMOST_CONFLICT
    almost_conflicted.push(0);
#endif
#if ANTI_EXPLORATION || BRANCHING_HEURISTIC == SGDB
    canceled.push(0);
#endif
#if BRANCHING_HEURISTIC == CHB || BRANCHING_HEURISTIC == SGDB
    last_conflict.push(0);
#endif
    total_actual_rewards.push(0);
    total_actual_count.push(0);
    if(set_decision_vars) dvar = decision[v] ? true : false; // EDXXX could probably be done better...
    if(popsim_branching)
    	dvar = false;
    setDecisionVar(v, dvar);
    return v;
}


bool Solver::addClause_(vec<Lit>& ps)
{
    assert(decisionLevel() == 0);
    if (!ok) return false;

    // Check if clause is satisfied and remove false/duplicate literals:
    sort(ps);
    Lit p; int i, j;
    for (i = j = 0, p = lit_Undef; i < ps.size(); i++)
        if (value(ps[i]) == l_True || ps[i] == ~p)
            return true;
        else if (value(ps[i]) != l_False && ps[i] != p)
            ps[j++] = p = ps[i];
    ps.shrink(i - j);

    if (ps.size() == 0)
        return ok = false;
    else if (ps.size() == 1){
        uncheckedEnqueue(ps[0]);
        return ok = (propagate() == CRef_Undef);
    }else{
        CRef cr = ca.alloc(ps, false);
        clauses.push(cr);
        attachClause(cr);
    }

    return true;
}


void Solver::attachClause(CRef cr) {
    const Clause& c = ca[cr];
    assert(c.size() > 1);
    watches[~c[0]].push(Watcher(cr, c[1]));
    watches[~c[1]].push(Watcher(cr, c[0]));
    if (c.learnt()) learnts_literals += c.size();
    else            clauses_literals += c.size(); }


void Solver::detachClause(CRef cr, bool strict) {
    const Clause& c = ca[cr];
    assert(c.size() > 1);
    
    if (strict){
        remove(watches[~c[0]], Watcher(cr, c[1]));
        remove(watches[~c[1]], Watcher(cr, c[0]));
    }else{
        // Lazy detaching: (NOTE! Must clean all watcher lists before garbage collecting this clause)
        watches.smudge(~c[0]);
        watches.smudge(~c[1]);
    }

    if (c.learnt()) learnts_literals -= c.size();
    else            clauses_literals -= c.size(); }


void Solver::removeClause(CRef cr) {
    Clause& c = ca[cr];
    detachClause(cr);
    // Don't leave pointers to free'd memory!
    if (locked(c)) vardata[var(c[0])].reason = CRef_Undef;
    c.mark(1); 
    ca.free(cr);
}


bool Solver::satisfied(const Clause& c) const {
    for (int i = 0; i < c.size(); i++)
        if (value(c[i]) == l_True)
            return true;
    return false; }


// Revert to the state at given level (keeping all assignment at 'level' but not beyond).
//
void Solver::cancelUntil(int level) {
    if (decisionLevel() > level){
        for (int c = trail.size()-1; c >= trail_lim[level]; c--){
            Var      x  = var(trail[c]);
            uint64_t age = conflicts - picked[x];
            if (age > 0) {
                double reward = ((double) conflicted[x]) / ((double) age);
#if BRANCHING_HEURISTIC == LRB
#if ALMOST_CONFLICT
                double adjusted_reward = ((double) (conflicted[x] + almost_conflicted[x])) / ((double) age);
#else
                double adjusted_reward = reward;
#endif
                double old_activity = activity[x];
                activity[x] = step_size * adjusted_reward + ((1 - step_size) * old_activity);
                if (order_heap.inHeap(x)) {
                    if (activity[x] > old_activity)
                        order_heap.decrease(x);
                    else
                        order_heap.increase(x);
                }
#endif
                total_actual_rewards[x] += reward;
                total_actual_count[x] ++;
            }
#if ANTI_EXPLORATION
            canceled[x] = conflicts;
#endif
            assigns [x] = l_Undef;

            if (phase_saving > 1 || (phase_saving == 1) && c > trail_lim.last())
                polarity[x] = sign(trail[c]);
            insertVarOrder(x); }
        qhead = trail_lim[level];
        trail.shrink(trail.size() - trail_lim[level]);
        trail_lim.shrink(trail_lim.size() - level);
    } }


//=================================================================================================
// Major methods:


Lit Solver::pickBranchLit()
{
	if(verification_mode){
		// need to absorb the current clause at the given asserting literal
		Lit a = checkAbsorptionStatus();
		if(a == lit_Undef){
			assert(final_sat_trail_mode);
		}
		else{
			return a;
		}
	}

    Var next = var_Undef;

    // Random decision:
    if (drand(random_seed) < random_var_freq && !order_heap.empty()){
        next = order_heap[irand(random_seed,order_heap.size())];
        if (value(next) == l_Undef && decision[next])
            rnd_decisions++; }

    // Activity based decision:

    while (next == var_Undef || value(next) != l_Undef || !decision[next])
        if (order_heap.empty()){
            next = var_Undef;
            break;
        } else {
#if ANTI_EXPLORATION
            next = order_heap[0];
            uint64_t age = conflicts - canceled[next];
            while (age > 0) {
                double decay = pow(0.95, age);
                activity[next] *= decay;
                if (order_heap.inHeap(next)) {
                    order_heap.increase(next);
                }
                canceled[next] = conflicts;
                next = order_heap[0];
                age = conflicts - canceled[next];
            }
#endif
#if BRANCHING_HEURISTIC == SGDB
			next = order_heap[0];
			uint64_t age = conflicts - canceled[next];
			while (age > 0 && value(next) == l_Undef) {
				double decay = pow(regularization, age);
				activity[next] *= decay;
				if (order_heap.inHeap(next)) {
					order_heap.increase(next);
					order_heap.decrease(next);
				}
				canceled[next] = conflicts;
				next = order_heap[0];
				age = conflicts - canceled[next];
			}
#endif
            next = order_heap.removeMin();
        }
    if(next != var_Undef)
    	all_decisions[next] = 1;
    if(verification_mode){
    	if(next != var_Undef){
    		printf("Verification failed, needed to branch on a variable after exhausting final trail\n");
    		exit(1);
    	}
    	else
    		printf("Verification succeeded\n");
    }


    return next == var_Undef ? lit_Undef : mkLit(next, rnd_pol ? drand(random_seed) < 0.5 : polarity[next]);
}

/*_________________________________________________________________________________________________
|
|  analyze : (confl : Clause*) (out_learnt : vec<Lit>&) (out_btlevel : int&)  ->  [void]
|  
|  Description:
|    Analyze conflict and produce a reason clause.
|  
|    Pre-conditions:
|      * 'out_learnt' is assumed to be cleared.
|      * Current decision level must be greater than root level.
|  
|    Post-conditions:
|      * 'out_learnt[0]' is the asserting literal at level 'out_btlevel'.
|      * If out_learnt.size() > 1 then 'out_learnt[1]' has the greatest decision level of the 
|        rest of literals. There may be others from the same level though.
|  
|________________________________________________________________________________________________@*/
void Solver::analyze(CRef confl, vec<Lit>& out_learnt, vec<Lit>& lsr_conflict_side, int& out_btlevel)
{
    int pathC = 0;
    Lit p     = lit_Undef;

    for (int i = 0; i < lsr_seen.size(); i++) assert(lsr_seen[i] == 0);


    // Generate conflict clause:
    //
    out_learnt.push();      // (leave room for the asserting literal)
    int index   = trail.size() - 1;

    do{
        assert(confl != CRef_Undef); // (otherwise should be UIP)
        Clause& c = ca[confl];

        // lsr -- grab the vars that any learnt depends upon
        if(c.learnt()){
        	for (int i = 0; i < c.rsize(); i++){
			  if (!lsr_seen[var(c[i])]){
				lsr_seen[var(c[i])] = 1;
				lsr_conflict_side.push(c[i]);
			  }
			}
        }

#if LBD_BASED_CLAUSE_DELETION
        if (c.learnt() && c.activity() > 2)
            c.activity() = lbd(c);
#else
        if (c.learnt())
            claBumpActivity(c);
#endif

        for (int j = (p == lit_Undef) ? 0 : 1; j < c.size(); j++){
            Lit q = c[j];

            if (!seen[var(q)] && level(var(q)) > 0){
#if BRANCHING_HEURISTIC == CHB
                last_conflict[var(q)] = conflicts;
#elif BRANCHING_HEURISTIC == VSIDS
                varBumpActivity(var(q));
#elif BRANCHING_HEURISTIC == SGDB
                in_conflict.push(var(q));
                last_conflict[var(q)] = conflicts;
#endif
                conflicted[var(q)]++;
                seen[var(q)] = 1;
                if (level(var(q)) >= decisionLevel())
                    pathC++;
                else
                    out_learnt.push(q);
            }
        }
        
        // Select next clause to look at:
        while (!seen[var(trail[index--])]);
        p     = trail[index+1];
        confl = reason(var(p));
        seen[var(p)] = 0;
        pathC--;

    }while (pathC > 0);
    out_learnt[0] = ~p;

    // Simplify conflict clause:
    //
    int i, j;
    out_learnt.copyTo(analyze_toclear);
    if (ccmin_mode == 2){
    	//printf("Exiting -- ccmin still broken for lsr\n");
    	//exit(1);
        uint32_t abstract_level = 0;
        for (i = 1; i < out_learnt.size(); i++)
            abstract_level |= abstractLevel(var(out_learnt[i])); // (maintain an abstraction of levels involved in conflict)

        for (i = j = 1; i < out_learnt.size(); i++)
            if (reason(var(out_learnt[i])) == CRef_Undef || !litRedundant(out_learnt[i], abstract_level, lsr_conflict_side))
                out_learnt[j++] = out_learnt[i];
        
    }else if (ccmin_mode == 1){
    	//printf("Exiting -- did not handle ccmin mode for lsr\n");
    	//exit(1);
        for (i = j = 1; i < out_learnt.size(); i++){
            Var x = var(out_learnt[i]);

            if (reason(x) == CRef_Undef)
                out_learnt[j++] = out_learnt[i];
            else{
                Clause& c = ca[reason(var(out_learnt[i]))];
                // lsr -- grab the vars that any learnt depends upon
				if(c.learnt()){
					for (int i = c.size(); i < c.rsize(); i++){
					  if (!lsr_seen[var(c[i])]){
						lsr_seen[var(c[i])] = 1;
						lsr_conflict_side.push(c[i]);
					  }
					}
				}
                for (int k = 1; k < c.size(); k++)
                    if (!seen[var(c[k])] && level(var(c[k])) > 0){
                        out_learnt[j++] = out_learnt[i];
                        break; }
            }
        }
    }else
        i = j = out_learnt.size();

    max_literals += out_learnt.size();
    out_learnt.shrink(i - j);
    tot_literals += out_learnt.size();

    // Find correct backtrack level:
    //
    if (out_learnt.size() == 1)
        out_btlevel = 0;
    else{
        int max_i = 1;
        // Find the first literal assigned at the next-highest level:
        for (int i = 2; i < out_learnt.size(); i++)
            if (level(var(out_learnt[i])) > level(var(out_learnt[max_i])))
                max_i = i;
        // Swap-in this literal at index 1:
        Lit p             = out_learnt[max_i];
        out_learnt[max_i] = out_learnt[1];
        out_learnt[1]     = p;
        out_btlevel       = level(var(p));
    }

#if ALMOST_CONFLICT
    seen[var(p)] = true;
    for(int i = out_learnt.size() - 1; i >= 0; i--) {
        Var v = var(out_learnt[i]);
        CRef rea = reason(v);
        if (rea != CRef_Undef) {
            Clause& reaC = ca[rea];
            for (int i = 0; i < reaC.size(); i++) {
                Lit l = reaC[i];
                if (!seen[var(l)]) {
                    seen[var(l)] = true;
                    almost_conflicted[var(l)]++;
                    analyze_toclear.push(l);
                }
            }
        }
    }
#endif
    for (i = 0; i < lsr_conflict_side.size(); i++) lsr_seen[var(lsr_conflict_side[i])] = 0;
    for (int j = 0; j < analyze_toclear.size(); j++) seen[var(analyze_toclear[j])] = 0;    // ('seen[]' is now cleared)
}


void Solver::getDecisionsFinalUnsat(CRef cr, vec<Lit>& decisions)
{
	Clause& c = ca[cr];
	vec<Lit> clause;
	// add the literals in the final conflicting clause itself
	for (int i = 0; i < c.size(); i++){
		clause.push(c[i]);
		if(!seen[var(c[i])]){
			seen[var(c[i])] = 1;
			decisions.push(c[i]);
		}
	}

  for (int i = 0; i < lsr_seen.size(); i++) assert(lsr_seen[i] == 0);

   // don't re-add vars from the conflict side learnts to decisions
   // (but still explore them by not adding to lsr_seen)
   for (int i = 0; i < decisions.size(); i++){
	  seen[var(decisions[i])] = 1;
  }

   bool print_flag = false;

  vec<Lit> workpool; clause.copyTo(workpool);

  vec<Lit> cert_decisions; // only used for certificate gen

  while (workpool.size() > 0){
    Var x = var(workpool.last());
    Lit p = workpool.last();
    if (lsr_seen[x]){
      workpool.pop();
      continue;
    }

    if(print_flag)
		printf("workpool %s%d\n", sign(p)?"-":"", var(p) + 1);


    assert(lsr_seen[x] == 0);
    lsr_seen[x] = 1;
    lsr_toclear.push(x);
    workpool.pop();

    if (assumption(x)){
    	if(print_flag)
			printf("assumption: %d\n", x + 1);
    	// i think this should only happen on the conflicting unit literal
    	if(reason(x) != CRef_Undef && print_flag){
    		printf("Assumption with reason\n");
    	}

      // Unit clause
      assert (unit_lsr[x] != CRef_Undef);
      Clause &cu = ca_lsr[unit_lsr[x]];
      //printf("size: %d\n",cu.size());
      for (int i = 0; i < cu.size(); i++){
        if (!seen[var(cu[i])]){
          seen[var(cu[i])] = 1;
          decisions.push(cu[i]);
        }
        // add the unit itself
        if(!seen[x]){
        	seen[x] = 1;
        	decisions.push(p);
        }

      }
    }
    else if (reason(x) == CRef_Undef){
      // Decision variable
    	if(level(x) > 0){
			printf("Should never have a decision in final conflict (unless level 0 due to no-pre). Error\n");
			printf("%s%d ", sign(p)?"-":"", x);
			printf("Level %d\n", level(x));
			exit(1);
    	}
    	continue;
    }
    // can hit this case for an assumption now
    if(reason(x) != CRef_Undef){
      // Reason is a clause

    	Clause &c = ca[reason(x)];
		if (c.learnt()){
			if(print_flag){
				  printf("learnt: ");
				  for (int i = 0; i < c.size(); i++){
					  printf("%s%d ", sign(c[i])?"-":"", var(c[i])+1);
				  }
				  printf("\n");
				  printf("ldeps: ");
				  for (int i = c.size(); i < c.rsize(); i++){
					  printf("%s%d ", sign(c[i])?"-":"", var(c[i])+ 1);
				  }
				  printf("\n");
			  }

			//the whole clause and dependencies
			for (int i = 0; i < c.rsize(); i++){
				if (!seen[var(c[i])]){
					seen[var(c[i])] = 1;
					decisions.push(c[i]);
				}
			}
		}
		 if(print_flag && !c.learnt()){
			  printf("clause!\n");
				  for (int i = 0; i < c.size(); i++){
					  printf("%s%d ", sign(c[i])?"-":"", var(c[i]) + 1);
				  }
			  printf("\n");
		  }
		// for all types of clauses, add the clause itself to the workpool
		for (int i = 0; i < c.size(); i++){
			if (!lsr_seen[var(c[i])]){
				workpool.push(c[i]);
		}
		}
    }
  }

  for (int i = 0; i < lsr_toclear.size(); i++) lsr_seen[lsr_toclear[i]] = 0;
  for (int i = 0; i < decisions.size(); i++){
    seen[var(decisions[i])] = 0;
    //printf("dec %d\n",var(decisions[i])+1);
  }
}


void Solver::getDecisions(vec<Lit>& clause, vec<Lit>& decisions, bool print_flag, bool final_sat_call)
{
  for (int i = 0; i < lsr_seen.size(); i++) assert(lsr_seen[i] == 0);

  vec<Lit> final_decisions; // only used to generate_certificate for final sat trail

  // now that decisions may start with lits, don't re-add var already in decisions
  for (int i = 0; i < decisions.size(); i++){
	  seen[var(decisions[i])] = 1;
	  lsr_toclear.push(var(decisions[i]));
  }
  if(clause_and_conflict_side_lsr){
	  if(!final_sat_call){
		for(int i = 0; i < clause.size(); i++){
		  Var v = var(clause[i]);
		  if(!seen[v]){
			  seen[v] = 1;
			  lsr_toclear.push(v);
			  decisions.push(clause[i]);
		  }
		}
		for (int i = 0; i < lsr_toclear.size(); i++) lsr_seen[lsr_toclear[i]] = 0;
		for (int i = 0; i < decisions.size(); i++){
			seen[var(decisions[i])] = 0;
			//printf("dec %d\n",var(decisions[i])+1);
		}

		lsr_toclear.clear();
		if(record_lsr_frequency){
			for (int i = 0; i < decisions.size(); i++){
				lsr_frequency[var(decisions[i])] += 1;
			}
		}
	  }
	  else{ // for the final call to the sat case TODO: make separate function
		assert(decisions.size() == 0);
		for (int i = 0; i < lsr_toclear.size(); i++) lsr_seen[lsr_toclear[i]] = 0;
		for (int i = 0; i < trail.size(); i++){
			Lit p = trail[i];
			Var x = var(p);
			if (assumption(x)){
				if(print_flag){
					printf("assumption: %d\n", x);
				}
				// Unit clause -- include its dependencies
				assert (unit_lsr[x] != CRef_Undef);
				Clause &cu = ca_lsr[unit_lsr[x]];
				for (int i = 0; i < cu.size(); i++){
					if (!seen[var(cu[i])]){
						seen[var(cu[i])] = 1;
						decisions.push(cu[i]);
					}
				}
				// also include the unit literal itself
				if(!seen[x]){
					seen[x] = true;
					decisions.push(p);
				}
			}
			else if (reason(x) == CRef_Undef){
				if(print_flag){
					printf("decision: %s%d\n", sign(p)?"-":"", var(p));
				}
				if(generate_certificate && final_sat_call)
					final_decisions.push(p);
				// Decision variable on the final trail
				if (!seen[x]){
					seen[x] = 1;
					decisions.push(p);
				}

			}
			else {
				// Reason is a clause

				Clause &c = ca[reason(x)];
				if (c.learnt()){
					if(print_flag){
						printf("learnt: ");
						for (int i = 0; i < c.size(); i++){
							printf("%s%d ", sign(c[i])?"-":"", var(c[i]));
						}
						printf("\n");
						printf("ldeps: ");
						for (int i = c.size(); i < c.rsize(); i++){
							printf("%s%d ", sign(c[i])?"-":"", var(c[i]));
						}
						printf("\n");
					}
					// include the learnt clause itself, as well as its dependencies
					// note: for the non-absorption approach below, i starts at c.size()
					for (int i = 0; i < c.rsize(); i++){
						if (!seen[var(c[i])]){
							seen[var(c[i])] = 1;
							decisions.push(c[i]);
						}
					}
				}
				else { // with the absorption result, I don't need to consider original clauses
					continue;
				}
			}
		}
	  }
	  if(generate_certificate){
		  // separate last sat call with an extra 0
		  if(final_sat_call){
			  fprintf(certificate_clauses_out, " 0\n");
			  for(int i = 0; i < final_decisions.size(); i++)
				  fprintf(certificate_clauses_out, "%s%d ", sign(final_decisions[i])?"-":"", var(final_decisions[i]) + 1);
			  fprintf(certificate_clauses_out, "0\n");
		  }
		  else{
			  bool all_lsr_clause = true;
			  // check if all decisions are in lsr_in
			  for (int i = 0; i < decisions.size(); i++){
				  if(!lsr_in[var(decisions[i])]){
					  all_lsr_clause = false;
					  break;
				  }
			  }


			  // if so, add the current clause to the certificate, which will be used to generate the current clause
			  if(all_lsr_clause){
				  for (int i = 0; i < clause.size(); i++){
					  fprintf(certificate_clauses_out, "%s%d ", sign(clause[i])?"-":"", var(clause[i]) + 1);
					  //printf("%s%d ", sign(clause[i])?"-":"", var(clause[i]) + 1);
				  }
				  fprintf(certificate_clauses_out, "0\n");
			  }
		  }
	  }
	  for (int i = 0; i < lsr_toclear.size(); i++) lsr_seen[lsr_toclear[i]] = 0;
	  for (int i = 0; i < decisions.size(); i++){
			seen[var(decisions[i])] = 0;
			//printf("dec %d\n",var(decisions[i])+1);
	  }
	  lsr_toclear.clear();
	  return;
  }

  printf("Only conf-side-lsr for now\n");
  exit(1);
}

Lit Solver::checkAbsorptionStatus(){
	Lit a = lit_Undef;
	// need to be careful with cancelUntil due to units
	int base_level = assumptions.size() + unit_assumptions.size();
	while(!final_sat_trail_mode && curr_cert_index < cert.size()){
		vec<Lit>* c = cert[curr_cert_index];
		if(c->size() == 0){
			//printf("Final sat decisions!\n");
			cancelUntil(base_level);
			final_sat_trail_mode = true;
			final_sat_trail = cert[curr_cert_index+1];
			break;
		}
		/*printf("Checking:");
		for(int i = 0; i < c->size(); i++){
			printf(" %s%d", sign((*c)[i])?"-":"", var((*c)[i])+1);
		}
		printf("\n");
		*/

		cancelUntil(base_level);
		// need to take care of units through assumptions
		vec<Lit> cl;
		bool clause_absorbed = false;

		for(int i = 0; i < c->size(); i++){
			Lit l = (*c)[i];
			if(value(l) == l_Undef){
				cl.push(l);
			}
			else if(value(l) == l_False){
				//printf("hit\n");
			}
			else
				clause_absorbed = true;
		}

		if(clause_absorbed){
			// printf("Clause is already implied by a unit assumption\n");
			curr_cert_index++;
			continue;
		}

		clause_absorbed = true;
		// find an asserting literal
		for(int i = 0; i < cl.size(); i++){
			cancelUntil(base_level);
			a = cl[i];
			for(int j = 0; j < cl.size(); j++){
				if(j == i)
					continue;
				Lit l = cl[j];

				if(value(l) == l_Undef){
					newDecisionLevel();
					uncheckedEnqueue(~l);
				}
			}

			CRef cr = propagate();
			if(cr == CRef_Undef && value(a) == l_Undef){
				newDecisionLevel();
				uncheckedEnqueue(~a);
				CRef cr = propagate();
				// should not happen
				if(cr == CRef_Undef){
					printf("Clause at index %d is not 1-provable? Missing clause bug?\n", curr_cert_index);
					printf("Units(%d): ", unit_assumptions.size());
					for(int i = 0; i < unit_assumptions.size(); i++){
						printf(" %s%d", sign(unit_assumptions[i])?"-":"", var(unit_assumptions[i])+1);
					}
					printf("\n");
					printf("Curr trail(%d): ", trail.size());
					for(int i = 0; i < trail.size(); i++){
						printf(" %s%d@%d", sign(trail[i])?"-":"", var(trail[i])+1, level(var(trail[i])));
					}
					printf("\n");
					exit(1);
				}
				else{
					clause_absorbed = false;
					// printf("%s%d is asserting\n", sign(a)?"-":"", var(a)+1);
					cancelUntil(base_level);

					for(int j = 0; j < cl.size(); j++){
						if(j == i)
							continue;
						Lit l = cl[j];
						if(value(l) == l_Undef){
							newDecisionLevel();
							uncheckedEnqueue(~l);
						}
					}
					propagate();

					assert(value(a) == l_Undef);
					break;
				}
			}
		}
		if(clause_absorbed){
			//printf("Clause already absorbed, trying next.\n");
			curr_cert_index++;
		}
		else
			break;
	}
	if(final_sat_trail_mode){
		/*printf("Curr trail(%d): ", trail.size());
		for(int i = 0; i < trail.size(); i++){
			printf(" %s%d", sign(trail[i])?"-":"", var(trail[i])+1);
		}
		printf("\n");
		*/
		while(true){
			if(final_sat_trail_index == final_sat_trail->size()){
				return lit_Undef;
			}
			Lit l = (*final_sat_trail)[final_sat_trail_index++];
			if(value(l) == l_Undef)
				return l;
		}
	}
	// final unsat case
	else if(curr_cert_index >= cert.size()){
		printf("this case should never be hit, as the analyzer will reach unsat before here\n");
		exit(1);
		/*
		printf("Attempting to derive final conflict\n");
		cancelUntil(0);
		for(int i = 0; i < unit_assumptions.size(); i++){
			// printf(" %s%d", sign(unit_assumptions[i])?"-":"", var(unit_assumptions[i])+1);
			uncheckedEnqueue(unit_assumptions[i]);
		}
		CRef cr = propagate();
		if(cr == CRef_Undef){
			printf("Fail\n");
		}
		else
			printf("Success\n");
		return lit_Undef;
		*/
	}
	else
		return ~a;
}


// Check if 'p' can be removed. 'abstract_levels' is used to abort early if the algorithm is
// visiting literals at levels that cannot be removed later.
bool Solver::litRedundant(Lit p, uint32_t abstract_levels, vec<Lit>& lsr_conflict_side)
{
    analyze_stack.clear(); analyze_stack.push(p);
    int top = analyze_toclear.size();
    while (analyze_stack.size() > 0){
        assert(reason(var(analyze_stack.last())) != CRef_Undef);
        Clause& c = ca[reason(var(analyze_stack.last()))]; analyze_stack.pop();
        // lsr -- grab the vars that any learnt depends upon
		if(c.learnt()){
			for (int i = c.size(); i < c.rsize(); i++){
			  if (!lsr_seen[var(c[i])]){
				lsr_seen[var(c[i])] = 1;
				lsr_conflict_side.push(c[i]);
			  }
			}
		}
        for (int i = 1; i < c.size(); i++){
            Lit p  = c[i];
            if (!seen[var(p)] && level(var(p)) > 0){
                if (reason(var(p)) != CRef_Undef && (abstractLevel(var(p)) & abstract_levels) != 0){
                    seen[var(p)] = 1;
                    analyze_stack.push(p);
                    analyze_toclear.push(p);
                }else{
                    for (int j = top; j < analyze_toclear.size(); j++)
                        seen[var(analyze_toclear[j])] = 0;
                    analyze_toclear.shrink(analyze_toclear.size() - top);
                    return false;
                }
            }
        }
    }

    return true;
}


/*_________________________________________________________________________________________________
|
|  analyzeFinal : (p : Lit)  ->  [void]
|  
|  Description:
|    Specialized analysis procedure to express the final conflict in terms of assumptions.
|    Calculates the (possibly empty) set of assumptions that led to the assignment of 'p', and
|    stores the result in 'out_conflict'.
|________________________________________________________________________________________________@*/
void Solver::analyzeFinal(Lit p, vec<Lit>& out_conflict)
{
    out_conflict.clear();
    out_conflict.push(p);

    if (decisionLevel() == 0)
        return;

    seen[var(p)] = 1;

    for (int i = trail.size()-1; i >= trail_lim[0]; i--){
        Var x = var(trail[i]);
        if (seen[x]){
            if (reason(x) == CRef_Undef){
                assert(level(x) > 0);
                out_conflict.push(~trail[i]);
            }else{
                Clause& c = ca[reason(x)];
                for (int j = 1; j < c.size(); j++)
                    if (level(var(c[j])) > 0)
                        seen[var(c[j])] = 1;
            }
            seen[x] = 0;
        }
    }

    seen[var(p)] = 0;
}


void Solver::uncheckedEnqueue(Lit p, CRef from)
{
    assert(value(p) == l_Undef);
    if(backbone_logging){
    	// if p is a backbone, check if its being flipped
    	Var x = var(p);
    	if(backbone[x] != 0){
    		int new_bb_pol;
    		if(!sign(p))
    			new_bb_pol = 0;
    		else
    			new_bb_pol = 1;
    		if(backbone_prev_polarity[x] != new_bb_pol)
    			num_backbone_flips++;
    		backbone_prev_polarity[x] = new_bb_pol;
    	}
    }
    picked[var(p)] = conflicts;
#if ANTI_EXPLORATION
    uint64_t age = conflicts - canceled[var(p)];
    if (age > 0) {
        double decay = pow(0.95, age);
        activity[var(p)] *= decay;
        if (order_heap.inHeap(var(p))) {
            order_heap.increase(var(p));
        }
    }
#endif
    conflicted[var(p)] = 0;
#if ALMOST_CONFLICT
    almost_conflicted[var(p)] = 0;
#endif
    assigns[var(p)] = lbool(!sign(p));
    vardata[var(p)] = mkVarData(from, decisionLevel(), assumption(var(p)));
    trail.push_(p);
}


/*_________________________________________________________________________________________________
|
|  propagate : [void]  ->  [Clause*]
|  
|  Description:
|    Propagates all enqueued facts. If a conflict arises, the conflicting clause is returned,
|    otherwise CRef_Undef.
|  
|    Post-conditions:
|      * the propagation queue is empty, even if there was a conflict.
|________________________________________________________________________________________________@*/
CRef Solver::propagate()
{
    CRef    confl     = CRef_Undef;
    int     num_props = 0;
    watches.cleanAll();

    while (qhead < trail.size()){
        Lit            p   = trail[qhead++];     // 'p' is enqueued fact to propagate.
        vec<Watcher>&  ws  = watches[p];
        Watcher        *i, *j, *end;
        num_props++;
        vec<Lit> enqueue_lits;
        vec<CRef> enqueue_cr;
        for (i = j = (Watcher*)ws, end = i + ws.size();  i != end;){
            // Try to avoid inspecting the clause:
            Lit blocker = i->blocker;
            if (value(blocker) == l_True){
                *j++ = *i++; continue; }

            // Make sure the false literal is data[1]:
            CRef     cr        = i->cref;
            Clause&  c         = ca[cr];
            Lit      false_lit = ~p;
            if (c[0] == false_lit)
                c[0] = c[1], c[1] = false_lit;
            assert(c[1] == false_lit);
            i++;

            // If 0th watch is true, then clause is already satisfied.
            Lit     first = c[0];
            Watcher w     = Watcher(cr, first);
            if (first != blocker && value(first) == l_True){
                *j++ = w; continue; }

            // Look for new watch:
            for (int k = 2; k < c.size(); k++)
                if (value(c[k]) != l_False){
                    c[1] = c[k]; c[k] = false_lit;
                    watches[~c[1]].push(w);
                    goto NextClause; }

            // Did not find watch -- clause is unit under assignment:
            *j++ = w;
            if (value(first) == l_False){
                confl = cr;
                qhead = trail.size();
                // Copy the remaining watches:
                while (i < end)
                    *j++ = *i++;
            }else{
				uncheckedEnqueue(first, cr);
            }
        NextClause:;
        }


        ws.shrink(i - j);
    }
    propagations += num_props;
    simpDB_props -= num_props;

    return confl;
}

int min(int a, int b) {
    return a < b ? a : b;
}




/*_________________________________________________________________________________________________
|
|  reduceDB : ()  ->  [void]
|  
|  Description:
|    Remove half of the learnt clauses, minus the clauses locked by the current assignment. Locked
|    clauses are clauses that are reason to some assignment. Binary clauses are never removed.
|________________________________________________________________________________________________@*/
struct reduceDB_lt { 
    ClauseAllocator& ca;
#if LBD_BASED_CLAUSE_DELETION
    vec<double>& activity;
    reduceDB_lt(ClauseAllocator& ca_,vec<double>& activity_) : ca(ca_), activity(activity_) {}
#else
    reduceDB_lt(ClauseAllocator& ca_) : ca(ca_) {}
#endif
    bool operator () (CRef x, CRef y) { 
#if LBD_BASED_CLAUSE_DELETION
        return ca[x].activity() > ca[y].activity();
    }
#else
        return ca[x].size() > 2 && (ca[y].size() == 2 || ca[x].activity() < ca[y].activity()); } 
#endif
};

struct reduceDB_lsr_size {
    ClauseAllocator& ca;
    vec<double>& activity;
    reduceDB_lsr_size(ClauseAllocator& ca_,vec<double>& activity_) : ca(ca_), activity(activity_) {}
    bool operator () (CRef x, CRef y) {
        return (ca[x].rsize() - ca[x].size()) > (ca[y].rsize() - ca[y].size());
    }
};

struct reduceDB_lsr_spread {
    ClauseAllocator& ca;
    vec<double>& activity;
    reduceDB_lsr_spread(ClauseAllocator& ca_,vec<double>& activity_) : ca(ca_), activity(activity_) {}
    bool operator () (CRef x, CRef y) {
    	printf("Need to implement LSR spread \n exiting");
    	exit(1);
        return (ca[x].rsize() - ca[x].size()) > (ca[y].rsize() - ca[y].size());
    }
};

void Solver::reduceDB()
{
    int     i, j;
#if LBD_BASED_CLAUSE_DELETION
    if(lsr_deletion_policy == 0)
    	sort(learnts, reduceDB_lt(ca, activity));
    else if(lsr_deletion_policy == 1)
    	sort(learnts, reduceDB_lsr_size(ca, activity));
    else if(lsr_deletion_policy == 2)
    	sort(learnts, reduceDB_lsr_spread(ca, activity));
#else
    double  extra_lim = cla_inc / learnts.size();    // Remove any clause below this activity
    sort(learnts, reduceDB_lt(ca));
#endif

    // Don't delete binary or locked clauses. From the rest, delete clauses from the first half
    // and clauses with activity smaller than 'extra_lim':
#if LBD_BASED_CLAUSE_DELETION
    for (i = j = 0; i < learnts.size(); i++){
        Clause& c = ca[learnts[i]];
        if (c.activity() > 2 && !locked(c) && i < learnts.size() / 2)
#else
    for (i = j = 0; i < learnts.size(); i++){
        Clause& c = ca[learnts[i]];
        if (c.size() > 2 && !locked(c) && (i < learnts.size() / 2 || c.activity() < extra_lim))
#endif
            removeClause(learnts[i]);
        else
            learnts[j++] = learnts[i];
    }
    learnts.shrink(i - j);
    checkGarbage();
}


void Solver::removeSatisfied(vec<CRef>& cs)
{
    int i, j;
    for (i = j = 0; i < cs.size(); i++){
        Clause& c = ca[cs[i]];
        if (satisfied(c))
            removeClause(cs[i]);
        else
            cs[j++] = cs[i];
    }
    cs.shrink(i - j);
}


void Solver::rebuildOrderHeap()
{
    vec<Var> vs;
    for (Var v = 0; v < nVars(); v++)
        if (decision[v] && value(v) == l_Undef)
            vs.push(v);
    order_heap.build(vs);
}


/*_________________________________________________________________________________________________
|
|  simplify : [void]  ->  [bool]
|  
|  Description:
|    Simplify the clause database according to the current top-level assigment. Currently, the only
|    thing done here is the removal of satisfied clauses, but more things can be put here.
|________________________________________________________________________________________________@*/
bool Solver::simplify()
{
    assert(decisionLevel() == 0);

    if (!ok || propagate() != CRef_Undef)
        return ok = false;

    if (nAssigns() == simpDB_assigns || (simpDB_props > 0))
        return true;

    // Remove satisfied clauses:
    removeSatisfied(learnts);
    if (remove_satisfied)        // Can be turned off.
        removeSatisfied(clauses);
    checkGarbage();
    rebuildOrderHeap();

    simpDB_assigns = nAssigns();
    simpDB_props   = clauses_literals + learnts_literals;   // (shouldn't depend on stats really, but it will do for now)

    return true;
}

/*_________________________________________________________________________________________________
|
|  search : (nof_conflicts : int) (params : const SearchParams&)  ->  [lbool]
|  
|  Description:
|    Search for a model the specified number of conflicts. 
|    NOTE! Use negative value for 'nof_conflicts' indicate infinity.
|  
|  Output:
|    'l_True' if a partial assigment that is consistent with respect to the clauseset is found. If
|    all variables are decision variables, this means that the clause set is satisfiable. 'l_False'
|    if the clause set is unsatisfiable. 'l_Undef' if the bound on number of conflicts is reached.
|________________________________________________________________________________________________@*/
lbool Solver::search(int nof_conflicts)
{
    assert(ok);
    int         backtrack_level;
    int         conflictC = 0;
    vec<Lit>    learnt_clause;
#if BRANCHING_HEURISTIC == SGDB
    vec<Var>    sampled;
#endif
    starts++;

    for (;;){
        CRef confl = propagate();
        //printTrail();


#if BRANCHING_HEURISTIC == CHB
        double multiplier = confl == CRef_Undef ? reward_multiplier : 1.0;
        for (int a = action; a < trail.size(); a++) {
            Var v = var(trail[a]);
            uint64_t age = conflicts - last_conflict[v] + 1;
            double reward = multiplier / age ;
            double old_activity = activity[v];
            activity[v] = step_size * reward + ((1 - step_size) * old_activity);
            if (order_heap.inHeap(v)) {
                if (activity[v] > old_activity)
                    order_heap.decrease(v);
                else
                    order_heap.increase(v);
            }
        }
#endif
        if (confl != CRef_Undef){
            // CONFLICT
            conflicts++; conflictC++;
#if BRANCHING_HEURISTIC == CHB || BRANCHING_HEURISTIC == LRB || BRANCHING_HEURISTIC == SGDB
            if (step_size > min_step_size)
                step_size -= step_size_dec;
#endif
            if (decisionLevel() == 0) return l_False;

#if BRANCHING_HEURISTIC == SGDB
            in_conflict.clear();
#endif
            learnt_clause.clear();
            vec<Lit> decision_clause;
            // analyze will now put the dependency lits from the conflict side in decision_clause
            analyze(confl, learnt_clause, decision_clause, backtrack_level);
#if BRANCHING_HEURISTIC == SGDB
            if (trail_lim.size() > 0) {
                for (int k = 0; k < learnt_clause.size(); k++) {
                    Var v = var(learnt_clause[k]);
                    CRef rea = reason(v);
                    if (rea != CRef_Undef) {
                        Clause& reaC = ca[rea];
                        for (int l = 0; l < reaC.size(); l++) {
                            Lit n = reaC[l];
                            Var nV = var(n);
                            if (last_conflict[nV] != conflicts) {
                                last_conflict[nV] = conflicts;
                                in_conflict.push(nV);
                            }
                        }
                    }
                }

                double conflict_class = bias;
                double non_conflict_class = bias;
                for (int i = 0; i < in_conflict.size(); i++) {
                    Var v = in_conflict[i];
                    uint64_t age = conflicts - canceled[v];
                    if (age > 1) {
                        double decay = pow(regularization, age - 1);
                        activity[v] *= decay;
                    }
                    canceled[v] = conflicts;
                    conflict_class += activity[v];
                }
                sampled.clear();
                for (int i = 0; i < trail_lim.size() - 1; i++) {
                    int start = trail_lim[i];
                    int end = trail_lim[i+1];
                    if (level(var(trail[start])) >= level(var(trail[end]))) {
                        printf("fail %d %d\n", level(var(trail[start])), level(var(trail[end])));
                        exit(1);
                    }
                    int t = start + irand(random_seed, end - start);
                    Var v = var(trail[t]);
                    if (canceled[v] != conflicts) {
                        sampled.push(v);
                        uint64_t age = conflicts - canceled[v];
                        if (age > 1) {
                            double decay = pow(regularization, age - 1);
                            activity[v] *= decay;
                        }
                        canceled[v] = conflicts;
                        non_conflict_class += activity[v];
                    }
                }
                in_conflicts += in_conflict.size();
                sampleds += sampled.size();

                double conflict_error = sigmoid(conflict_class) - 1;
                double non_conflict_error = sigmoid(non_conflict_class);

                bias = bias * regularization - step_size * (conflict_error + non_conflict_error);

                for (int i = 0; i < in_conflict.size(); i++) {
                    Var v = in_conflict[i];
                    activity[v] = activity[v] * regularization - step_size * conflict_error;
                    if (order_heap.inHeap(v)) {
                        order_heap.increase(v);
                        order_heap.decrease(v);
                    }
                }
                for (int i = 0; i < sampled.size(); i++) {
                    Var v = sampled[i];
                    activity[v] = activity[v] * regularization - step_size * non_conflict_error;
                    if (order_heap.inHeap(v)) {
                        order_heap.increase(v);
                        order_heap.decrease(v);
                    }
                }
            }
#endif
            if(cmty_logging){
            	int s = learnt_clause.size();
            	for(int i = 0; i < s; i++){
            		int cmty = var_cmty[var(learnt_clause[i])];
            		cmty_clauses[cmty] = cmty_clauses[cmty] + (float(1) / s);
            	}
            }
            if(backbone_logging){
				// if any of the literals in the clause are in the backbone,
            	// then the backbone lit would subsume it
            	for(int i = 0; i < learnt_clause.size(); i++){
            		Lit l = learnt_clause[i];
            		Var x = var(l);
            		if(backbone[x] != 0){
            			if((sign(l) && backbone[x] == -1) || (!sign(l) && backbone[x] == 1)){
            				num_backbone_subsumed_clauses++;
            				break;
            			}
            		}
            	}
			}


            /*
            if(verification_mode){
				printf("Learned: ");
				for(int i = 0 ; i < learnt_clause.size(); i++){
					printf("%s%d ", sign(learnt_clause[i])?"-":"", var(learnt_clause[i]) + 1);
				}
				printf("\n");
            }
            */
            if (learnt_clause.size() == 1){
                unit_assumptions.push(learnt_clause[0]);
                assert (confl != CRef_Undef);
                // if the condition is false, this is the decision level 0 conflict
                if(!assumption(var(learnt_clause[0]))){
                	getDecisions(learnt_clause, decision_clause);
                    assert (decision_clause.size() > 0);
                }

                if (assumption(var(learnt_clause[0]))){
                  assert(unit_lsr[var(learnt_clause[0])]!=CRef_Undef);
                  Clause &c = ca_lsr[unit_lsr[var(learnt_clause[0])]];
                  //printf("Other side: ");
                  for (int i = 0; i < c.size(); i++){
                    Var x = var(c[i]);
                    //printf(" %d", x);

                    seen[x] = 1;
                    lsr_final.push(x);
                  }

                  getDecisionsFinalUnsat(confl, decision_clause);


                  // don't want to add duplicates to lsr_final
                  for (int i = 0; i < lsr_final.size(); i++){
                	  seen[lsr_final[i]] = 1;
                  }
                  for (int i = 0; i < decision_clause.size(); i++){
                    Var x = var(decision_clause[i]);
                    if (!seen[x]){
                      lsr_final.push(x);
                      seen[x] = 1;
                    } 
                  }

                  // clean up
                  for (int i = 0; i < lsr_final.size(); i++)
                    seen[lsr_final[i]] = 0;

                  if(verification_mode){
						printf("Verification succeeded\n");
				  }
                  return l_False;
                }



                setAssumption(var(learnt_clause[0]), true);

                CRef cr = ca_lsr.alloc(decision_clause, false);
                unit_lsr[var(learnt_clause[0])] = cr;
                //unit_lsr.push(cr);

                //printf("variable %d, assumption %d\n",var(learnt_clause[0])+1,assumption(var(learnt_clause[0])));

                cancelUntil(backtrack_level);

#if BRANCHING_HEURISTIC == CHB
            action = trail.size();
#endif

                //uncheckedEnqueue(learnt_clause[0]);
            }else{
                //CRef cr = ca.alloc(learnt_clause, true);

                getDecisions(learnt_clause, decision_clause);
                /*
                if(verification_mode){
					printf("(");
					for (int i = 0; i < decision_clause.size(); i++){
					  printf("%s%d ", sign(decision_clause[i])?"-":"", var(decision_clause[i]) + 1);
					}
					printf(")");
					printf("\n-----------------------\n");
				}*/
                cancelUntil(backtrack_level);

#if BRANCHING_HEURISTIC == CHB
            action = trail.size();
#endif

                //vec<Lit> decision_clause;
                int size = learnt_clause.size();
                for (int i = 0; i < decision_clause.size(); i++)
                  learnt_clause.push(decision_clause[i]);

                CRef cr = ca.alloc(learnt_clause, size, learnt_clause.size(), true);
                learnts.push(cr);
                attachClause(cr);
#if LBD_BASED_CLAUSE_DELETION
                Clause& clause = ca[cr];
                clause.activity() = lbd(clause);
#else
                claBumpActivity(ca[cr]);
#endif
                uncheckedEnqueue(learnt_clause[0], cr);
            }

            if(avg_clause_lsr_out){
				all_learnts++;
				total_clause_lsr_weight += decision_clause.size();
			}


#if BRANCHING_HEURISTIC == VSIDS
            varDecayActivity();
#endif
#if ! LBD_BASED_CLAUSE_DELETION
            claDecayActivity();
#endif

            if (--learntsize_adjust_cnt == 0){
                learntsize_adjust_confl *= learntsize_adjust_inc;
                learntsize_adjust_cnt    = (int)learntsize_adjust_confl;
#if ! RAPID_DELETION
                max_learnts             *= learntsize_inc;
#endif

                if (verbosity >= 1)
                    printf("| %9d | %7d %8d %8d | %8d %8d %6.0f | %6.3f %% |\n", 
                           (int)conflicts, 
                           (int)dec_vars - (trail_lim.size() == 0 ? trail.size() : trail_lim[0]), nClauses(), (int)clauses_literals, 
                           (int)max_learnts, nLearnts(), (double)learnts_literals/nLearnts(), progressEstimate()*100);
            }
            /*
            if (always_restart){
				// Reached bound on number of conflicts:
            	printf("Warning always restart....\n");
				restart_now = false;
				if(generate_certificate && !just_restarted){
					fprintf(certificate_clauses_out, "0 0 0\n");
					just_restarted = true;
				}
				progress_estimate = progressEstimate();
				cancelUntil(0);
				return l_Undef;
            }
            */

        }else{
            // NO CONFLICT
            if (restart_now || (nof_conflicts >= 0 && conflictC >= nof_conflicts || !withinBudget())){
                // Reached bound on number of conflicts:
            	restart_now = false;
                progress_estimate = progressEstimate();
                cancelUntil(0);
                return l_Undef; }

            // Simplify the set of problem clauses:
            if (decisionLevel() == 0 && !simplify())
                return l_False;

            if (learnts.size()-nAssigns() >= max_learnts && !never_gc) {
                // Reduce the set of learnt clauses:
            	//printf("reduce %f\n", max_learnts);
                reduceDB();
#if RAPID_DELETION
                max_learnts += 500;
#endif
            }

            Lit next = lit_Undef;
            while (decisionLevel() < assumptions.size() + unit_assumptions.size()){
                // Perform user provided assumption:
                Lit p = lit_Undef;
                if (decisionLevel() < unit_assumptions.size()){
                  p = unit_assumptions[decisionLevel()];
                  //printf("Asserting unit: %s%d\n", sign(p)?"-":"", var(p));
                }
                else 
                  p = assumptions[decisionLevel()];
                assert (p != lit_Undef);
                if (value(p) == l_True){
                    // Dummy decision level:
                    newDecisionLevel();
                }else if (value(p) == l_False){
                    analyzeFinal(~p, conflict);
                    return l_False;
                }else{
                    next = p;
                    break;
                }
            }

            if (next == lit_Undef){
                // New variable decision:
                decisions++;
                //if(verification_mode)
                //	checkAbsorptionStatus();

                next = pickBranchLit();
                //printf("var %d\n", var(next));
                //if(next == lit_Undef)
                //    printf("picked lit: undef\n");
                //else
		// 	printf("picked lit: %s%d\n", sign(next)?"-":"", var(next) + 1);
                if(structure_logging){
                	if(cmty_logging && var(next) >= 0){
                		int cmty = var_cmty[var(next)];
                		//printf("cmtytest %d %d %d\n", cmty_picks.size(), cmty, var(next));
                		cmty_picks[cmty] = cmty_picks[cmty] + 1;

                	}

                }

                if (restart_immediately){
					// Reached bound on number of conflicts:
                	printf("Preemptive restart\n");
					restart_now = false;
					restart_immediately = false;
					progress_estimate = progressEstimate();
					cancelUntil(0);
					return l_Undef;
                }

                if (next == lit_Undef){
                	if(set_decision_vars || popsim_branching){
                		bool actually_failed = false;
                		if(trail.size() != nVars()){
                			//if(verbosity > 0)
                			//	printf("Failed %d %d\n", trail.size(), nVars());
                			for(int i = 0; i < nVars(); i++){
								if(assigns[i] == l_Undef){
									vec<Watcher>& pos = watches[mkLit(i, true)];
									vec<Watcher>& neg = watches[mkLit(i, false)];
									if(pos.size() != 0 || neg.size() != 0)
										actually_failed = true;
								}
							}
                			if(actually_failed){
								focused_branching_remaining_failures--;
                				if(embed_lsr){
                					if(focused_branching_remaining_failures <= 0){
                						vec<Lit> fake_learnt;
										if(embed_lsr_clause_type == 0) // randomly grab literals from the trail
										{
											// copy trail to auxilliary vec randomized
											vec<Lit> trail_lits;
											trail_lits.growTo(trail.size(), lit_Undef);
											// definitely a better way to do this, but works for now
											int num_inserted = 0;
											for(int i = 0; i < trail.size(); i++){
												Lit l = trail[i];
												int index = irand(random_seed, trail.size() - num_inserted);
												int unused_slots = -1;
												for(int j = 0; j < trail_lits.size(); j++){
													if(trail_lits[j] == lit_Undef) // the spot in the randomized trail is not yet assigned
														unused_slots++;
													if(index == unused_slots){
														trail_lits[j] = l;
														break;
													}
												}
												num_inserted++;
											}
											// this part is the same as case clause_type == 1
											int s = 0;
											while(fake_learnt.size() < embed_lsr_target_clause_size && s < trail_lits.size()){
												if(level(var(trail_lits[s])) != 0)
													fake_learnt.push(~trail_lits[s]);
												s++;
											}
										}
										else if(embed_lsr_clause_type == 1) // start from lowest level literals
										{
											int s = 0;
											while(fake_learnt.size() < embed_lsr_target_clause_size && s < trail.size()){
												if(level(var(trail[s])) != 0)
													fake_learnt.push(~trail[s]);
												s++;
											}
										}
										else if(embed_lsr_clause_type == 2) // start from highest level literals
										{
											int s = trail.size() - 1;
											while(fake_learnt.size() < embed_lsr_target_clause_size && s >= 0){
												fake_learnt.push(~trail[s]);
												s--;
											}
										}

										assert(fake_learnt.size() > 1);
										printClause(fake_learnt, false);
										cancelUntil(0);

										CRef cr = ca.alloc(fake_learnt, fake_learnt.size(), fake_learnt.size(), true);
										learnts.push(cr);
										attachClause(cr);
										Clause& clause = ca[cr];
										// want to keep all these clauses, so just set their lbd to 2
										clause.activity() = 2;
                					}
								}

								cancelUntil(0);
								return l_Undef;
                			}
                		}
                	}
                    vec<Lit> clause;
                    for (int i = 0; i < nVars(); i++){
                    	if(assigns[i] != l_Undef) // some lits not in the trail
                    		clause.push(mkLit(i,false));
                    }
                    vec<Lit> decisions;

                    getDecisions(clause, decisions, false, true);

                    for (int i = 0; i < decisions.size(); i++)
                      lsr_final.push(var(decisions[i]));

                    if(lsr_final_deps_file)
                    	dumpFinalDepsFile();
                    //printLSR();

                    // Model found:
                    return l_True;
                  }
            }

            // Increase decision level and enqueue 'next'
            newDecisionLevel();
#if BRANCHING_HEURISTIC == CHB
            action = trail.size();
#endif

            uncheckedEnqueue(next);
        }
    }
}

void Solver::dumpFinalDepsFile(){
	// dumps all dependencies of each variable individually (can't use getDecisions directly due to 'seen')
	FILE* res_log = NULL;
	if (lsr_final_deps_file != NULL){
	  // Note: variables start at index 0
	  res_log = fopen(lsr_final_deps_file, "wb");
	}
	vec<Lit> decisions;

	for (int i = 0; i < trail.size(); i++){

		Lit p = trail[i];
		Var x = var(p);
		if (assumption(x)){
			decisions.push(p);
		}
		else if (reason(x) == CRef_Undef){
			decisions.push(p);
		}
		else {
			// Reason is a clause
			Clause &c = ca[reason(x)];
			if (c.learnt()){
				decisions.push(p);
			}
			else { // with the absorption result, I don't need to consider original clauses
				continue;
			}
		}
	}

	for(int i = 0; i < decisions.size(); i++){
		fprintf(res_log,"%d\n", var(decisions[i]));
	}
	fclose(res_log);
	for (int i = 0; i < lsr_seen.size(); i++) assert(lsr_seen[i] == 0);
}



double Solver::progressEstimate() const
{
    double  progress = 0;
    double  F = 1.0 / nVars();

    for (int i = 0; i <= decisionLevel(); i++){
        int beg = i == 0 ? 0 : trail_lim[i - 1];
        int end = i == decisionLevel() ? trail.size() : trail_lim[i];
        progress += pow(F, i) * (end - beg);
    }

    return progress / nVars();
}

/*
  Finite subsequences of the Luby-sequence:

  0: 1
  1: 1 1 2
  2: 1 1 2 1 1 2 4
  3: 1 1 2 1 1 2 4 1 1 2 1 1 2 4 8
  ...


 */

static double luby(double y, int x){

    // Find the finite subsequence that contains index 'x', and the
    // size of that subsequence:
    int size, seq;
    for (size = 1, seq = 0; size < x+1; seq++, size = 2*size+1);

    while (size-1 != x){
        size = (size-1)>>1;
        seq--;
        x = x % size;
    }

    return pow(y, seq);
}

// NOTE: assumptions passed in member-variable 'assumptions'.
lbool Solver::solve_()
{
    model.clear();
    conflict.clear();
    if (!ok) return l_False;
    solves++;

#if RAPID_DELETION
    max_learnts               = initial_max_learnts; // 2000
#else
    max_learnts               = nClauses() * learntsize_factor;
#endif
    learntsize_adjust_confl   = learntsize_adjust_start_confl;
    learntsize_adjust_cnt     = (int)learntsize_adjust_confl;
    lbool   status            = l_Undef;

    if (verbosity >= 1){
        printf("LBD Based Clause Deletion : %d\n", LBD_BASED_CLAUSE_DELETION);
        printf("Rapid Deletion : %d\n", RAPID_DELETION);
        printf("Almost Conflict : %d\n", ALMOST_CONFLICT);
        printf("Anti Exploration : %d\n", ANTI_EXPLORATION);
        printf("============================[ Search Statistics ]==============================\n");
        printf("| Conflicts |          ORIGINAL         |          LEARNT          | Progress |\n");
        printf("|           |    Vars  Clauses Literals |    Limit  Clauses Lit/Cl |          |\n");
        printf("===============================================================================\n");
    }

    // Search:
    int curr_restarts = 0;
    while (status == l_Undef){

    	if(popsim_branching || embed_lsr){
    		if(focused_branching_remaining_failures <= 0){
				if(popsim_branching){
					printf("Increasing focused branching limit %d\n", focused_branching_failure_limit);
					setDecisionVar(focused_branching_failure_limit, true);
				}
				focused_branching_remaining_failures = focused_branching_hard_limit; // focused_branching_failure_limit < 50 ? focused_branching_failure_limit : 50;
				focused_branching_failure_limit++;
				//printf("Setting focused branching remaining failures %d\n", focused_branching_remaining_failures);

    		}
    	}


    	double rest_base;
		if(always_restart)
			rest_base = 1;
		else if(never_restart)
			rest_base = -1;
		else if(uniform_restarts > 0)
			rest_base = uniform_restarts;
		else
			rest_base = (luby_restart ? luby(restart_inc, curr_restarts) : pow(restart_inc, curr_restarts)) * restart_first;
		status = search(rest_base);
        if (!withinBudget()) break;
        curr_restarts++;
    }

    if (verbosity >= 1)
        printf("===============================================================================\n");


    if (status == l_True){
        // Extend & copy model:
        model.growTo(nVars());
        for (int i = 0; i < nVars(); i++) model[i] = value(i);
    }else if (status == l_False && conflict.size() == 0)
        ok = false;

    // clean up
    for (int i = 0; i < unit_assumptions.size(); i++){
      CRef cr = unit_lsr[var(unit_assumptions[i])];
      if (cr == CRef_Undef)
        continue;
      assert (cr != CRef_Undef);
      Clause& c = ca_lsr[cr];
      c.mark(1);
      ca_lsr.free(cr);          
      unit_lsr[var(unit_assumptions[i])] = CRef_Undef;
    }
    
    unit_assumptions.clear();

    if (lsr_num){
      printf("LSR Backdoors [%d / %d] : %.2f%%\n",lsr_final.size(),nVars(),(double)lsr_final.size()/nVars()*100);
    }
    if(all_decs_num){
    	int count = 0;
      for(int i = 0; i < all_decisions.size(); i++)
    	  if(all_decisions[i])
    		  count++;

      printf("All Decisions [%d / %d] : %.2f%%\n",count, nVars(),(double)count/nVars()*100);

    }

    if (lsr_filename != NULL){
      // Note: variables start at index 0
      FILE* res_log = fopen(lsr_filename, "wb");
      for (int i = 0; i < lsr_final.size(); i++){
        fprintf(res_log, "%d\n",lsr_final[i]);
      }
    }

    if(all_decisions_filename != NULL){
    	// Note: variables start at index 0
		FILE* res_log = fopen(all_decisions_filename, "wb");
		for (int i = 0; i < all_decisions.size(); i++){
			if(all_decisions[i])
				fprintf(res_log, "%d\n", i);
		}

    }

    lsr_final.clear();

    cancelUntil(0);
    return status;
}

//=================================================================================================
// Writing CNF to DIMACS:
// 
// FIXME: this needs to be rewritten completely.

static Var mapVar(Var x, vec<Var>& map, Var& max)
{
    if (map.size() <= x || map[x] == -1){
        map.growTo(x+1, -1);
        map[x] = max++;
    }
    return map[x];
}


void Solver::toDimacs(FILE* f, Clause& c, vec<Var>& map, Var& max)
{
    if (satisfied(c)) return;

    for (int i = 0; i < c.size(); i++)
        if (value(c[i]) != l_False)
            fprintf(f, "%s%d ", sign(c[i]) ? "-" : "", mapVar(var(c[i]), map, max)+1);
    fprintf(f, "0\n");
}


void Solver::toDimacs(const char *file, const vec<Lit>& assumps)
{
    FILE* f = fopen(file, "wr");
    if (f == NULL)
        fprintf(stderr, "could not open file %s\n", file), exit(1);
    toDimacs(f, assumps);
    fclose(f);
}


void Solver::toDimacs(FILE* f, const vec<Lit>& assumps)
{
    // Handle case when solver is in contradictory state:
    if (!ok){
        fprintf(f, "p cnf 1 2\n1 0\n-1 0\n");
        return; }

    vec<Var> map; Var max = 0;

    // Cannot use removeClauses here because it is not safe
    // to deallocate them at this point. Could be improved.
    int cnt = 0;
    for (int i = 0; i < clauses.size(); i++)
        if (!satisfied(ca[clauses[i]]))
            cnt++;
        
    for (int i = 0; i < clauses.size(); i++)
        if (!satisfied(ca[clauses[i]])){
            Clause& c = ca[clauses[i]];
            for (int j = 0; j < c.size(); j++)
                if (value(c[j]) != l_False)
                    mapVar(var(c[j]), map, max);
        }

    // Assumptions are added as unit clauses:
    cnt += assumptions.size();

    fprintf(f, "p cnf %d %d\n", max, cnt);

    for (int i = 0; i < assumptions.size(); i++){
        assert(value(assumptions[i]) != l_False);
        fprintf(f, "%s%d 0\n", sign(assumptions[i]) ? "-" : "", mapVar(var(assumptions[i]), map, max)+1);
    }

    for (int i = 0; i < clauses.size(); i++)
        toDimacs(f, ca[clauses[i]], map, max);

    if (verbosity > 0)
        printf("Wrote %d clauses with %d variables.\n", cnt, max);
}


//=================================================================================================
// Garbage Collection methods:

void Solver::relocAll(ClauseAllocator& to)
{
    // All watchers:
    //
    // for (int i = 0; i < watches.size(); i++)
    watches.cleanAll();
    for (int v = 0; v < nVars(); v++)
        for (int s = 0; s < 2; s++){
            Lit p = mkLit(v, s);
            // printf(" >>> RELOCING: %s%d\n", sign(p)?"-":"", var(p)+1);
            vec<Watcher>& ws = watches[p];
            for (int j = 0; j < ws.size(); j++)
                ca.reloc(ws[j].cref, to);
        }

    // All reasons:
    //
    for (int i = 0; i < trail.size(); i++){
        Var v = var(trail[i]);

        if (reason(v) != CRef_Undef && (ca[reason(v)].reloced() || locked(ca[reason(v)])))
            ca.reloc(vardata[v].reason, to);
    }

    // All learnt:
    //
    for (int i = 0; i < learnts.size(); i++)
        ca.reloc(learnts[i], to);

    // All original:
    //
    for (int i = 0; i < clauses.size(); i++)
        ca.reloc(clauses[i], to);
}


void Solver::garbageCollect()
{
    // Initialize the next region to a size corresponding to the estimated utilization degree. This
    // is not precise but should avoid some unnecessary reallocations for the new region:
    ClauseAllocator to(ca.size() - ca.wasted()); 

    relocAll(to);
    if (verbosity >= 2)
        printf("|  Garbage collection:   %12d bytes => %12d bytes             |\n", 
               ca.size()*ClauseAllocator::Unit_Size, to.size()*ClauseAllocator::Unit_Size);
    to.moveTo(ca);
}
