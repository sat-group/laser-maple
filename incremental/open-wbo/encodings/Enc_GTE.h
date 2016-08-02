#ifndef Enc_GTE_h
#define Enc_GTE_h

#ifdef SIMP
#include "simp/SimpSolver.h"
#else
#include "core/Solver.h"
#endif

#include "core/SolverTypes.h"
#include "Encodings.h"
#include <map>
#include <vector>
#include <utility>

namespace openwbo
{
struct wlitt
{
	Lit lit;
	uint64_t weight;
};
struct less_than_wlitt
{
	inline bool operator()(const wlitt& wl1, const wlitt& wl2)
	{
		return (wl1.weight < wl2.weight);
	}
};
struct wlit_sumt
{
	inline uint64_t operator()(const uint64_t& wl1, const wlitt& wl2)
	{
		return (wl1 + wl2.weight);
	}
};
typedef std::map<uint64_t,Lit> wlit_mapt;
typedef std::vector<wlitt> weightedlitst;
typedef std::pair<uint64_t,Lit> wlit_pairt;
class GTE : public Encodings
{

public:
  GTE()
  {
    //current_pb_rhs = -1; // -1 corresponds to an unitialized value
    current_pb_rhs = 0;
    nb_clauses = 0;
    nb_variables = 0;
  }
 ~GTE() {}

  // Encode constraint.
  void encode(Solver *S, vec<Lit> &lits, vec<uint64_t> &coeffs, uint64_t rhs);
  
  // Update constraint.
  void update(Solver *S, uint64_t rhs);
  
  // Returns true if the encoding was built, otherwise returns false;
  bool hasCreatedEncoding() { return hasEncoding; }

protected:
  void printLit(Lit l){
    printf("%s%d\n", sign(l) ? "-" : "", var(l)+1);
  }

  bool encodeLeq(uint64_t k, Solver *S, const weightedlitst& iliterals,wlit_mapt& oliterals);
  Lit getNewLit(Solver *S);
  Lit get_var(Solver* S,wlit_mapt& oliterals,uint64_t weight);
  vec<Lit> pb_outlits;    // Stores the outputs of the pseudo-Boolean constraint
                          // encoding for incremental solving.
  uint64_t current_pb_rhs; // Stores the current value of the rhs of the
                          // pseudo-Boolean constraint.

  // Stores unit lits. Used for lits that have a coeff larger than rhs.
  wlit_mapt pb_oliterals;
  vec<Lit> unit_lits;
  vec<uint64_t> unit_coeffs;

  // Number of variables and clauses for statistics.
  int nb_variables;
  int nb_clauses;
};

}

#endif
