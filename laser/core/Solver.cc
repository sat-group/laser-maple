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

#if BRANCHING_HEURISTIC == CHB || BRANCHING_HEURISTIC == LRB
static DoubleOption  opt_step_size         (_cat, "step-size",   "Initial step size",                             0.40,     DoubleRange(0, false, 1, false));
static DoubleOption  opt_step_size_dec     (_cat, "step-size-dec","Step size decrement",                          0.000001, DoubleRange(0, false, 1, false));
static DoubleOption  opt_min_step_size     (_cat, "min-step-size","Minimal step size",                            0.06,     DoubleRange(0, false, 1, false));
#endif
#if BRANCHING_HEURISTIC == VSIDS
static DoubleOption  opt_var_decay         (_cat, "var-decay",   "The variable activity decay factor",            0.95,     DoubleRange(0, false, 1, false));
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

// lsr computation types
static BoolOption    opt_clause_and_conflict_side_lsr      (_cat, "conf-side-lsr",        "Dependencies of a clause are the clause itself and the dependents on the conflict side.", false);



//=================================================================================================
// Constructor/Destructor:


Solver::Solver() :

    // Parameters (user settable):
    //
    verbosity        (0)
#if BRANCHING_HEURISTIC == CHB || BRANCHING_HEURISTIC == LRB
  , step_size        (opt_step_size)
  , step_size_dec    (opt_step_size_dec)
  , min_step_size    (opt_min_step_size)
#endif
#if BRANCHING_HEURISTIC == VSIDS
  , var_decay        (opt_var_decay)
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
{
  //strcpy(lsr_filename,"");
  lsr_filename = NULL;
  all_decisions_filename = NULL;
  lsr_num = false;
} 


Solver::~Solver()
{
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
#if ANTI_EXPLORATION
    canceled.push(0);
#endif
#if BRANCHING_HEURISTIC == CHB
    last_conflict.push(0);
#endif
    total_actual_rewards.push(0);
    total_actual_count.push(0);
    if(set_decision_vars) dvar = decision[v] ? true : false; // EDXXX could probably be done better...
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
            next = order_heap.removeMin();
        }
    /*for(int i = 0; i < decision.size(); i++){
    				printf("%d ", decision[i]);
    			}
    			printf("\n");
	printf("hit %d %d \n", next, decision[next]);
	*/
    //printf("%d\n", next);
    if(next != var_Undef)
    	all_decisions[next] = 1;
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
        // TODO warning -- this now contains the clause itself (i = 0 instead of i = c.size()) and its dependencies, can simplify getDecisions()
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
    	printf("Exiting -- ccmin still broken for lsr\n");
    	exit(1);
        uint32_t abstract_level = 0;
        for (i = 1; i < out_learnt.size(); i++)
            abstract_level |= abstractLevel(var(out_learnt[i])); // (maintain an abstraction of levels involved in conflict)

        for (i = j = 1; i < out_learnt.size(); i++)
            if (reason(var(out_learnt[i])) == CRef_Undef || !litRedundant(out_learnt[i], abstract_level, lsr_conflict_side))
                out_learnt[j++] = out_learnt[i];
        
    }else if (ccmin_mode == 1){
    	printf("Exiting -- did not handle ccmin mode for lsr\n");
    	exit(1);
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



void Solver::getDecisions(vec<Lit>& clause, vec<Lit>& decisions, bool print_flag, bool final_sat_call)
{
  for (int i = 0; i < lsr_seen.size(); i++) assert(lsr_seen[i] == 0);

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
	  }
	  else{ // for the final call to the sat case
		assert(decisions.size() == 0);
		for (int i = 0; i < lsr_toclear.size(); i++) lsr_seen[lsr_toclear[i]] = 0;
		for (int i = 0; i < trail.size(); i++){
			seen[var(decisions[i])] = 0;
			//printf("dec %d\n",var(decisions[i])+1);
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
	  return;
  }

  vec<Lit> workpool; clause.copyTo(workpool);
  while (workpool.size() > 0){
    Var x = var(workpool.last());
    Lit p = workpool.last();
    if (lsr_seen[x]){
    	if(print_flag){
    		printf("Skipping %d\n", x);
    	}
      workpool.pop();
      continue;
    }
    if(print_flag)
    	printf("workpool %s%d\n", sign(p)?"-":"", var(p));

    assert (lsr_seen[x] == 0);
    lsr_seen[x] = 1;
    lsr_toclear.push(x);
    workpool.pop();

    if (assumption(x)){
    	if(print_flag){
    		printf("assumption: %d\n", x);

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
      }
    } else if (reason(x) == CRef_Undef){
    	if(print_flag){
    		printf("decision: %s%d\n", sign(p)?"-":"", var(p));
    	}
      // Decision variable
      if (!seen[x]){
        seen[x] = 1;
        decisions.push(p);
      }
    } else {
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
        for (int i = c.size(); i < c.rsize(); i++){
          if (!seen[var(c[i])]){
            seen[var(c[i])] = 1;
            decisions.push(c[i]);
          }
        }
      } else {
    	  if(print_flag){
    		  printf("clause!\n");
    		  for (int i = 0; i < c.size(); i++){
				  printf("%s%d ", sign(c[i])?"-":"", var(c[i]));
			  }
			  printf("\n");
    	  }
        for (int i = 0; i < c.size(); i++){
          if (!lsr_seen[var(c[i])] && !seen[var(c[i])]){
            workpool.push(c[i]);
          }
        }
      }

    }
  }

  assert (decisions.size() > 0);
  for (int i = 0; i < lsr_toclear.size(); i++) lsr_seen[lsr_toclear[i]] = 0;
  for (int i = 0; i < decisions.size(); i++){
    seen[var(decisions[i])] = 0;
    //printf("dec %d\n",var(decisions[i])+1);
  }

  lsr_toclear.clear();
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
            }else
                uncheckedEnqueue(first, cr);

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
void Solver::reduceDB()
{
    int     i, j;
#if LBD_BASED_CLAUSE_DELETION
    sort(learnts, reduceDB_lt(ca, activity));
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
    starts++;

    for (;;){
        CRef confl = propagate();

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
#if BRANCHING_HEURISTIC == CHB || BRANCHING_HEURISTIC == LRB
            if (step_size > min_step_size)
                step_size -= step_size_dec;
#endif
            if (decisionLevel() == 0) return l_False;

            learnt_clause.clear();
            vec<Lit> decision_clause;
            // analyze will now put the dependency lits from the conflict side in decision_clause
            analyze(confl, learnt_clause, decision_clause, backtrack_level);
            //printf("lsize: %d %d\n", learnt_clause.size(), decision_clause.size());
            //decision_clause.clear();

            if (learnt_clause.size() == 1){
                unit_assumptions.push(learnt_clause[0]);
                // vec<Lit> decision_clause;


                assert (confl != CRef_Undef);
                //printf("------------------------------------\n");
                getDecisions(confl, decision_clause);
                //printf("Deps: ");
                //for(int i = 0; i < decision_clause.size(); i++)
                //	printf("%d ", var(decision_clause[i]));
                //printf("\n");
                //printf("Unit Lit: %s%d\n", sign(learnt_clause[0])?"-":"", var(learnt_clause[0]));

                assert (decision_clause.size() > 0);

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
                  //printf("\n");
                  //printf("Conflict on %d\n", var(learnt_clause[0]));
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

                  //printLSR();
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

                cancelUntil(backtrack_level);

#if BRANCHING_HEURISTIC == CHB
            action = trail.size();
#endif

                //vec<Lit> decision_clause;
                int size = learnt_clause.size();
                getDecisions(learnt_clause, decision_clause);
                //assert (confl != CRef_Undef);
                //getDecisions(confl, decision_clause);
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

        }else{
            // NO CONFLICT
            if (nof_conflicts >= 0 && conflictC >= nof_conflicts || !withinBudget()){
                // Reached bound on number of conflicts:
                progress_estimate = progressEstimate();
                cancelUntil(0);
                return l_Undef; }

            // Simplify the set of problem clauses:
            if (decisionLevel() == 0 && !simplify())
                return l_False;

            if (learnts.size()-nAssigns() >= max_learnts) {
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
                next = pickBranchLit();

                if (next == lit_Undef){
                	if(set_decision_vars){
                		bool actually_failed = false;
                		if(trail.size() != nVars()){
                			if(verbosity > 0)
                				printf("Failed %d %d\n", trail.size(), nVars());
                			for(int i = 0; i < nVars(); i++){
								if(assigns[i] == l_Undef){
									//printf("%d \n", i);
									vec<Watcher>& pos = watches[mkLit(i, true)];
									vec<Watcher>& neg = watches[mkLit(i, false)];
									//printf("w size : %d %d\n", pos.size(), neg.size());
									if(pos.size() != 0 || neg.size() != 0)
										actually_failed = true;
								}
							}
                			if(actually_failed){
								cancelUntil(0);
								return l_Undef;
                			}
                		}
                	}
                    vec<Lit> clause;
                    for (int i = 0; i < nVars(); i++)
                      clause.push(mkLit(i,false));
                    vec<Lit> decisions;

                    getDecisions(clause, decisions, false, true);
                    for (int i = 0; i < decisions.size(); i++)
                      lsr_final.push(var(decisions[i]));

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
    max_learnts               = 2000;
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
    	double rest_base;
		if(always_restart)
			rest_base = 1;
		else if(never_restart)
			rest_base = -1;
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
