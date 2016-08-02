/*****************************************************************************************[Main.cc]
MiniSat  -- Copyright (c) 2003-2006, Niklas Een, Niklas Sorensson
            Copyright (c) 2007-2010, Niklas Sorensson
Open-WBO -- Copyright (c) 2013-2015, Ruben Martins, Vasco Manquinho, Ines Lynce

Permission is hereby granted, free of charge, to any person obtaining a copy of
this software and associated documentation files (the "Software"), to deal in
the Software without restriction, including without limitation the rights to
use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of
the Software, and to permit persons to whom the Software is furnished to do so,
subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS
FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR
COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER
IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
************************************************************************************************/

#include <errno.h>
#include <signal.h>
#include <zlib.h>
#include "utils/System.h"
#include "utils/ParseUtils.h"
#include "utils/Options.h"

#ifdef SIMP
#include "simp/SimpSolver.h"
#else
#include "core/Solver.h"
#endif

#include "MaxTypes.h"
#include "MaxSAT.h"
#include "ParserMaxSAT.h"
#include "ParserPB.h"

//#define LSU 

// Algorithms
//#include "algorithms/Alg_WBO.h"
//#include "algorithms/Alg_incWBO.h"
#include "algorithms/Alg_LinearSU.h"
//#include "algorithms/Alg_LinearUS.h"
//#include "algorithms/Alg_WLinearUS.h"
//#include "algorithms/Alg_MSU3.h"
//#include "algorithms/Alg_WMSU3.h"
//#include "algorithms/Alg_Model.h"
#include "algorithms/Alg_OLL.h"
//#include "algorithms/Alg_PM2.h"

#define VER1_(x) #x
#define VER_(x) VER1_(x)
#define SATVER VER_(SOLVERNAME)
#define VER VER_(VERSION)

using NSPACE::cpuTime;
using NSPACE::OutOfMemoryException;
using NSPACE::IntOption;
using NSPACE::IntRange;
using NSPACE::parseOptions;
using namespace openwbo;

//=================================================================================================

static MaxSAT *mxsolver;

static void SIGINT_exit(int signum)
{
  mxsolver->printAnswer(_UNKNOWN_);
  exit(_UNKNOWN_);
}

//=================================================================================================
// Main:

int main(int argc, char **argv)
{
  printf(
	 "c\nc Open-WBO:\t a Modular MaxSAT Solver -- based on %s (%s version)\n",
	 SATVER, VER);
  //printf("c Version:\t 1.3.2 -- 11 June 2015\n");
#ifdef LSU
  printf("c Version:\t Pseudo-Boolean Evaluation 2015 -- LSU\n");
#else
  printf("c Version:\t Pseudo-Boolean Evaluation 2015\n");
#endif
  printf("c Authors:\t Ruben Martins, Vasco Manquinho, Ines Lynce\n");
  printf("c Contributors:\t Miguel Neves, Saurabh Joshi, Mikolas Janota\n");
  printf("c Contact:\t open-wbo@sat.inesc-id.pt -- "
         "http://sat.inesc-id.pt/open-wbo/\nc\n");
  try {
    NSPACE::setUsageHelp("c USAGE: %s [options] <input-file>\n\n");

#if defined(__linux__)
    fpu_control_t oldcw, newcw;
    _FPU_GETCW(oldcw);
    newcw = (oldcw & ~_FPU_EXTENDED) | _FPU_DOUBLE;
    _FPU_SETCW(newcw);
    printf(
	   "c WARNING: for repeatability, setting FPU to use double precision\n");
#endif

    IntOption algorithm("Open-WBO", "algorithm",
			"Search algorithm (0=oll, 1=linear-su).\n",
			1, IntRange(0, 1));

    IntOption format("Open-WBO","format","Input format (0=MaxSAT, 1=Pseudo-Boolean).\n)",
		     0, IntRange(0, 1));

    parseOptions(argc, argv, true);    

#if 0
    IntOption cpu_lim("Open-WBO", "cpu-lim",
                      "Limit on CPU time allowed in seconds.\n", INT32_MAX,
                      IntRange(0, INT32_MAX));
    IntOption mem_lim("Open-WBO", "mem-lim",
                      "Limit on memory usage in megabytes.\n", INT32_MAX,
                      IntRange(0, INT32_MAX));
    IntOption verbosity("Open-WBO", "verbosity",
                        "Verbosity level (0=minimal, 1=more).\n", 1,
                        IntRange(0, 1));
    IntOption algorithm("Open-WBO", "algorithm",
                        "Search algorithm (0=wbo, 1=linear-su, 2=linear-us, "
                        "3=wlinear-us, "
                        "4=msu3, 5=wmsu3,6=model, 7=oll, 8=pm2, 9=best).\n",
                        8, IntRange(0, 9));
    IntOption incremental("Open-WBO", "incremental",
                          "Incremental level (0=none, 1=blocking, 2=weakening, "
                          "3=iterative) (only for unsat-based algorithms).\n",
                          3, IntRange(0, 3));

    BoolOption bmo("Open-WBO", "bmo", "BMO search.\n", true);

    IntOption cardinality("Encodings", "cardinality",
                          "Cardinality encoding (0=cardinality networks, "
                          "1=totalizer, 2=modulo totalizer).\n",
                          1, IntRange(0, 2));

    IntOption amo("Encodings", "amo", "AMO encoding (0=Ladder).\n", 0,
                  IntRange(0, 0));

    IntOption pb("Encodings", "pb", "PB encoding (0=SWC).\n", 0,
                 IntRange(0, 0));

    IntOption weight(
		     "WBO", "weight-strategy",
		     "Weight strategy (0=none, 1=weight-based, 2=diversity-based).\n", 2,
		     IntRange(0, 2));
    BoolOption symmetry("WBO", "symmetry", "Symmetry breaking.\n", true);
    IntOption symmetry_lim(
			   "WBO", "symmetry-limit",
			   "Limit on the number of symmetry breaking clauses.\n", 500000,
			   IntRange(0, INT32_MAX));

    parseOptions(argc, argv, true);
#else
    int cpu_lim = INT32_MAX;
    int mem_lim = INT32_MAX;
    int verbosity = 0;
    //int algorithm = 1;
    int incremental = 3;
    bool bmo = true;
    int cardinality = 2;
    int amo = 0;
    int pb = 1;
    int weight = 2;
    bool symmetry = true;
    int symmetry_lim = 500000;
#endif

    //IntOption pb("Encodings", "pb", "PB encoding (0=SWC, 1=GTE).\n", 0,
    //             IntRange(0, 1));

    //parseOptions(argc, argv, true);

    double initial_time = cpuTime();
    MaxSAT *S = NULL;

    /*
      switch ((int)algorithm)
      {
      case _ALGORITHM_WBO_:
      if (incremental == _INCREMENTAL_NONE_)
      S = new WBO(verbosity, weight, symmetry, symmetry_lim);
      else if (incremental == _INCREMENTAL_BLOCKING_)
      S = new incWBO(verbosity, weight, symmetry, symmetry_lim, amo);
      else
      {
      printf(
      "c Error: WBO algorithm only supports blocking incrementality.\n");
      printf("s UNKNOWN\n");
      exit(_ERROR_);
      }
      break;

      case _ALGORITHM_LINEAR_SU_:
      S = new LinearSU(verbosity, bmo, cardinality, pb);
      break;

      case _ALGORITHM_LINEAR_US_:
      S = new LinearUS(verbosity, incremental, cardinality);
      break;

      case _ALGORITHM_WLINEAR_US_:
      S = new WLinearUS(verbosity, incremental, cardinality);
      break;

      case _ALGORITHM_MSU3_:
      S = new MSU3(verbosity, incremental, cardinality);
      break;

      case _ALGORITHM_WMSU3_:
      S = new WMSU3(verbosity, incremental, cardinality, pb, bmo);
      break;

      case _ALGORITHM_MODEL_:
      S = new Model(verbosity, incremental, cardinality);
      break;

      case _ALGORITHM_OLL_:
      S = new OLL(verbosity, incremental, cardinality);
      break;

      case _ALGORITHM_PM2_:
      S = new PM2(verbosity, incremental, cardinality);
      break;

      case _ALGORITHM_BEST_:
      S = new MSU3(_VERBOSITY_SOME_, _INCREMENTAL_ITERATIVE_, _CARD_TOTALIZER_);
      break;

      default:
      printf("c Error: Invalid MaxSAT algorithm.\n");
      printf("s UNKNOWN\n");
      exit(_ERROR_);
      }
    */

    if (algorithm == 1) 
      S = new LinearSU(verbosity, bmo, cardinality, pb);
    else {
      algorithm = _ALGORITHM_OLL_;
      cardinality = 1;
      S = new OLL(verbosity, incremental, cardinality);
    }
      
      mxsolver = S;
      mxsolver->setInitialTime(initial_time);

      signal(SIGXCPU, SIGINT_exit);
      signal(SIGTERM, SIGINT_exit);

      // Set limit on CPU-time:
      if (cpu_lim != INT32_MAX) {
	rlimit rl;
	getrlimit(RLIMIT_CPU, &rl);
	if (rl.rlim_max == RLIM_INFINITY || (rlim_t)cpu_lim < rl.rlim_max) {
	  rl.rlim_cur = cpu_lim;
	  if (setrlimit(RLIMIT_CPU, &rl) == -1)
	    printf("c WARNING! Could not set resource limit: CPU-time.\n");
	}
      }

      // Set limit on virtual memory:
      if (mem_lim != INT32_MAX) {
	rlim_t new_mem_lim = (rlim_t)mem_lim * 1024 * 1024;
	rlimit rl;
	getrlimit(RLIMIT_AS, &rl);
	if (rl.rlim_max == RLIM_INFINITY || new_mem_lim < rl.rlim_max) {
	  rl.rlim_cur = new_mem_lim;
	  if (setrlimit(RLIMIT_AS, &rl) == -1)
	    printf("c WARNING! Could not set resource limit: Virtual memory.\n");
	}
      }

      if (argc == 1){
	printf("c Error: no filename.\n");
	printf("c UNKNOWN\n");
	exit(_ERROR_);
      }
      
      gzFile in = (argc == 1) ? gzdopen(0, "rb") : gzopen(argv[1], "rb");
      if (in == NULL)
	printf("c ERROR! Could not open file: %s\n",
	       argc == 1 ? "<stdin>" : argv[1]),
	  printf("c UNKNOWN\n"), exit(_ERROR_);

      MaxSATFormula *maxsat_formula = new MaxSATFormula();
      ParserPB *parser_pb = new ParserPB();

      if (format == 0){
	parseMaxSATFormula(in, maxsat_formula);
	maxsat_formula->setFormat(_FORMAT_MAXSAT_);
      } else {
	parser_pb->parsePBFormula(argv[1], maxsat_formula);
	maxsat_formula->setFormat(_FORMAT_PB_);
      }

      gzclose(in);
      mxsolver->loadFormula(maxsat_formula);
  

      printf("c |                                                                "
	     "                                       |\n");
      printf("c ========================================[ Problem Statistics "
	     "]===========================================\n");
      printf("c |                                                                "
	     "                                       |\n");

      if (maxsat_formula->getProblemType() == _UNWEIGHTED_)
	printf("c |  Problem Type:  %19s                                         "
	       "                          |\n",
	       "Unweighted");
      else
	printf("c |  Problem Type:  %19s                                         "
	       "                          |\n",
	       "Weighted");

      printf("c |  Number of variables:  %12d                                    "
	     "                               |\n",
	     maxsat_formula->nVars());
      printf("c |  Number of hard clauses:    %7d                                "
	     "                                   |\n",
	     maxsat_formula->nHard());
      printf("c |  Number of soft clauses:    %7d                                "
	     "                                   |\n",
	     maxsat_formula->nSoft());
      printf("c |  Number of cardinality:     %7d                                "
	     "                                   |\n",
	     maxsat_formula->nCard());
      printf("c |  Number of PB :             %7d                                "
	     "                                   |\n",
	     maxsat_formula->nPB());
      double parsed_time = cpuTime();

      printf("c |  Parse time:           %12.2f s                                "
	     "                                 |\n",
	     parsed_time - initial_time);
      printf("c |                                                                "
	     "                                       |\n");

      /*
      //algorithm = _ALGORITHM_MSE15_;
      if (algorithm == _ALGORITHM_MSE15_)
      {
      MaxSAT *solver;
      if (S->nHard() == 0)
      {
      //MSU3
      solver = new MSU3(_VERBOSITY_MINIMAL_, _INCREMENTAL_ITERATIVE_, _CARD_TOTALIZER_);
      S->copySolver(solver);
      delete S;
      mxsolver = solver;
      mxsolver->setInitialTime(initial_time);
      mxsolver->search();
      } else {
      if (S->getProblemType() == _WEIGHTED_)
      {
      bool bmo_strategy = S->isBMO(false);
      if (bmo_strategy) {
      solver = new WMSU3(_VERBOSITY_MINIMAL_, _INCREMENTAL_ITERATIVE_,
      _CARD_TOTALIZER_, _PB_SWC_, true);

      S->copySolver(solver);
      delete S;
      mxsolver = solver;
      mxsolver->setInitialTime(initial_time);
      }
      mxsolver->search();
      } else {
      //default OLL
      mxsolver->search();
      }
      }
      }
      else if (algorithm == _ALGORITHM_BEST_ && S->getProblemType() == _WEIGHTED_)
      {
      // Check if the formula is BMO without caching the bmo functions.
      bool bmo_strategy = S->isBMO(false);
      MaxSAT *solver;
      if (bmo_strategy)
      solver = new WMSU3(_VERBOSITY_SOME_, _INCREMENTAL_ITERATIVE_,
      _CARD_TOTALIZER_, _PB_SWC_, true);
      else
      solver = new WBO(_VERBOSITY_SOME_, _WEIGHT_DIVERSIFY_, true, 500000);

      // HACK: Copy the contents of a solver to a fresh solver.
      // This could be avoided if the parsing was done to a DB that is
      // independent from the solver.
      S->copySolver(solver);
      delete S;
      mxsolver = solver;
      mxsolver->setInitialTime(initial_time);
      mxsolver->search();
      }
      else

      */
      mxsolver->search();


    } catch (OutOfMemoryException &) {
      sleep(1);
      printf("c Error: Out of memory.\n");
      printf("s UNKNOWN\n");
      exit(_ERROR_);
    }
  }

