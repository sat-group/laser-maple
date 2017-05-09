# Laser

Suite of tools to compute LSR-backdoors, built on top of the MapleSat SAT solver.

The laser directory contains the main tool, which heuristically computes upper bounds on LSR-backdoors, along with computing many other metrics.

The minlaser directory computes exact minimum LS and LSR backdoors, effectively by performing BFS over all search-tree explorations.


# Typical Usage:

Assume throughout that:
```
name=`basename my_sat_instance.cnf .cnf`
```

We also typically use -no-pre in favor of pre-simplifying most instances. 
For the lens experiments, it is also important that the metrics (e.g. community structure) are computed over the simplified instance.


For computing LSR-backdoor upper bounds:
```
./laser/simp/laser
   -verb=0                                 // lower output verbosity
   -cpu-lim=18000                          // cpu limit in seconds 
   -mem-lim=6600                           // memory limit in MB 
   -no-pre                                 // turn off preprocessing 
   -lsr-num                                // output the size of the found LSR-backdoor to stdout 
   -lsr-out=${name}.lsr                    // out file for LSR-backdoor variables (not minimized) 
   ${name}.cnf                             // the cnf file
```

For computing LSR-backdoor upper bounds, and further minimizing backdoors for satisfiable instances:
```
./laser/simp/laser
   -verb=0                                 // lower output verbosity
   -cpu-lim=18000                          // cpu limit in seconds
   -mem-lim=6600                           // memory limit in MB
   -no-pre                                 // turn off preprocessing
   -lsr-num                                // output the size of the found LSR-backdoor to stdout
   -lsr-out=${name}.lsr                    // out file for LSR-backdoor variables (not minimized)
   -lsr-final-deps=${name}.min_sat_lsr     // out file for LSR-backdoor variables (minimized, only works in SAT case)
   ${name}.cnf                             // the cnf file
```

For computing lens experiments:
```
./laser/simp/laser
   -verb=0                                 // lower output verbosity
   -cpu-lim=18000                          // cpu limit in seconds
   -mem-lim=6600                           // memory limit in MB
   -no-pre                                 // turn off preprocessing
   -lsr-num                                // output the size of the found LSR-backdoor to stdout
   -lsr-out=${name}.lsr                    // out file for LSR-backdoor variables
   -all-dec-out=${name}.all_decs           // out file for all variables that were branched on
   -avg-clause-lsr-out=${name}.avg_lsr     // out file for average dependency set size of all learnt clauses
   -cmty-loc-out=${name}.cmty_loc          // out file for community based locality
   -cmty-file=${name}.cmty                 // in file for community map from variables to communities
   -backbone-file=${name}.backbone         // in file for backbone
   -bb-metrics-out=${name}.bb_metrics      // out file for backbone-based measures
   -always-restart                         // restart after every conflict
   ${name}.cnf                             // the cnf file
```

For computing LSR-backdoors and further verifying their correctness (see laser/simp/run_and_verify.sh):
```
./laser/simp/laser 
   -lsr-out=${name}.lsr                    // out file for LSR-backdoor variables
   ${name}.cnf                             // the cnf file
   
./laser/simp/laser  
   -lsr-in=${name}.lsr                     // in file of LSR-backdoor variables (from previous call to laser)
   -lsr-cert-cls-out=${name}.clauses       // out file of all clauses dependent only on LSR-backdoor variables
   ${name}.cnf                             // the cnf file

./laser/simp/laser  
   -lsr-in=${name}.lsr                     // in file of LSR-backdoor variables (from previous call to laser)
   -lsr-cert-cls-in=${name}.clauses        // in file of learnt clauses (from previous call to laser)
   -never-gc                               // disable clause database reductions
   ${name}.cnf                             // the cnf file
```

For computing minimum LSR-backdoors:
```
./minlaser/simp/minlaser
   -lsr-mode=true                          // allow restarts in the extended trail (turn off to compute LS-backdoors)
   -min-lsr-size=2                         // only consider sets of [potential] backdoors of size at least 2
   -max-lsr-size=4                         // only consider sets of [potential] backdoors of size at least 4
   -decision-vars=${name}.vars             // if given, only look for backdoors over the listed variables
   -replay-bd-out=${name}.decisions        // out file of the sequence of decisions that witness the backdoor
   ${name}.cnf                             // the cnf file
```

If a minimum backdoor is found, output will contain a sequence of decisions that witness the backdoor. To run the solver
against the sequence of decisions (effectively overriding the solver's branching heuristic):
```
./minlaser/simp/minlaser
   -replay-bd-file=${name}.decisions       // in file of sequence of decisions (from previous call to minlaser)
   ${name}.cnf                             // the cnf file
```



# Options for laser (many are inherited from MapleSat):
```
USAGE: ./laser [options] <input-file> <result-output-file>

  where input may be either in plain or gzipped DIMACS.

CORE OPTIONS:

  -conf-side-lsr, -no-conf-side-lsr       (default: on)

        Dependencies of a clause are the clause itself and the dependents on the conflict side.

  -never-gc, -no-never-gc                 (default: off)

        Never remove clauses.

  -never-restart, -no-never-restart       (default: off)

        Restart never.

  -always-restart, -no-always-restart     (default: off)

        Restart after every conflict.

  -luby, -no-luby                         (default: on)

        Use the Luby restart sequence

  -rnd-pol, -no-rnd-pol                   (default: off)

        Randomize the polarity selection

  -rnd-init, -no-rnd-init                 (default: off)

        Randomize the initial activity


  -step-size    = <double> (   0 ..    1) (default: 0.4)

        Initial step size

  -step-size-dec = <double> (   0 ..    1) (default: 1e-06)

        Step size decrement

  -min-step-size = <double> (   0 ..    1) (default: 0.06)

        Minimal step size

  -rnd-freq     = <double> [   0 ..    1] (default: 0)

        The frequency with which the decision heuristic tries to choose a random variable

  -rnd-seed     = <double> (   0 ..  inf) (default: 9.16483e+07)

        Used by the random variable selection

  -gc-frac      = <double> (   0 ..  inf) (default: 0.2)

        The fraction of wasted memory allowed before a garbage collection is triggered

  -rinc         = <double> (   1 ..  inf) (default: 2)

        Restart interval increase factor


  -rfirst       = <int32>  [   1 .. imax] (default: 100)

        The base restart interval

  -phase-saving = <int32>  [   0 ..    2] (default: 2)

        Controls the level of phase saving (0=none, 1=limited, 2=full)

  -ccmin-mode   = <int32>  [   0 ..    2] (default: 2)

        Controls conflict clause minimization (0=none, 1=basic, 2=deep)


LASER OPTIONS:

  -all-decs-num, -no-all-decs-num         (default: off)

        Number of unique decision variables.


  -lsr-num, -no-lsr-num                   (default: off)

        Number of LSR backdoor variables.



  -backbone-file = <string>

        backbone literals (one-based)

  -cmty-file  = <string>

        var+cmty pairs, where vars should be 0 based

  -lsr-cert-cls-in = <string>

        File containing the clause sequence witnessing a backdoor (one-based).


  -lsr-cert-cls-out = <string>

        File to output the clauses sequence witnessing a backdoor (one-based).


  -lsr-in     = <string>

        Used to create a certificate for an lsr, generate with -lsr-out (zero-based).


  -bb-metrics-out = <string>

        Output backbone-based metrics to given file

  -cmty-loc-out = <string>

        Output locality measures to given file

  -all-dec-out = <string>

        Write all unique decision vars to a file (same as LS paper) (zero-based).


  -lsr-out    = <string>

        Write LSR backdoor to a file (zero-based).


  -avg-clause-lsr-out = <string>

        For each learnt, record its lsr size, compute the average. Dump to given file.

  -lsr-frequency-out = <string>

        Record how many times each variable is a dependent of a clause.

  -lsr-final-deps = <string>

        Smaller LSR for SAT case. Don't look at deps of propagated vars, just add it.


MAIN OPTIONS:

  -pre, -no-pre                           (default: on)

        Completely turn on/off any preprocessing.


  -mem-lim      = <int32>  [   0 .. imax] (default: 2147483647)

        Limit on memory usage in megabytes.


  -verb         = <int32>  [  -1 ..    2] (default: 1)

        Verbosity level (0=silent, 1=some, 2=more).

  -cpu-lim      = <int32>  [   0 .. imax] (default: 2147483647)

        Limit on CPU time allowed in seconds.



  -assumptions = <string>

        If given, use the assumptions in the file.

  -dimacs     = <string>

        If given, stop after preprocessing and write the result to this file.

  -decision-vars = <string>

        Only branch on the listed vars (zero-based).


SIMP OPTIONS:

  -asymm, -no-asymm                       (default: off)

        Shrink clauses by asymmetric branching.

  -rcheck, -no-rcheck                     (default: off)

        Check if a clause is already implied. (costly)

  -elim, -no-elim                         (default: on)

        Perform variable elimination.


  -simp-gc-frac = <double> (   0 ..  inf) (default: 0.5)

        The fraction of wasted memory allowed before a garbage collection is triggered during simplification.


  -cl-lim       = <int32>  [  -1 .. imax] (default: 20)

        Variables are not eliminated if it produces a resolvent with a length above this limit. -1 means no limit

  -sub-lim      = <int32>  [  -1 .. imax] (default: 1000)

        Do not check if subsumption against a clause larger than this. -1 means no limit.

  -grow         = <int32>  [imin .. imax] (default: 0)

        Allow a variable elimination step to grow by a number of clauses.


HELP OPTIONS:

  --help        Print help message.
  --help-verb   Print verbose help message.
```

# Options for minlaser:
```
USAGE: ./minlaser [options] <input-file> <result-output-file>

  where input may be either in plain or gzipped DIMACS.

CORE OPTIONS:

  -lsr-mode, -no-lsr-mode                 (default: off)

        Whether restarts are allowed in trail.

  -clause-deletion, -no-clause-deletion   (default: off)

        Whether clause deletion happens.

  -never-restart, -no-never-restart       (default: on)

        Never restart.

  -always-restart, -no-always-restart     (default: off)

        Restart after every conflict.

  -luby, -no-luby                         (default: on)

        Use the Luby restart sequence

  -rnd-pol, -no-rnd-pol                   (default: off)

        Randomize the polarity selection

  -rnd-init, -no-rnd-init                 (default: off)

        Randomize the initial activity


  -step-size    = <double> (   0 ..    1) (default: 0.4)

        Initial step size

  -step-size-dec = <double> (   0 ..    1) (default: 1e-06)

        Step size decrement

  -min-step-size = <double> (   0 ..    1) (default: 0.06)

        Minimal step size

  -rnd-freq     = <double> [   0 ..    1] (default: 0)

        The frequency with which the decision heuristic tries to choose a random variable

  -rnd-seed     = <double> (   0 ..  inf) (default: 9.16483e+07)

        Used by the random variable selection

  -rinc         = <double> (   1 ..  inf) (default: 2)

        Restart interval increase factor

  -gc-frac      = <double> (   0 ..  inf) (default: 0.2)

        The fraction of wasted memory allowed before a garbage collection is triggered


  -rfirst       = <int32>  [   1 .. imax] (default: 100)

        The base restart interval

  -ccmin-mode   = <int32>  [   0 ..    2] (default: 2)

        Controls conflict clause minimization (0=none, 1=basic, 2=deep)

  -phase-saving = <int32>  [   0 ..    2] (default: 2)

        Controls the level of phase saving (0=none, 1=limited, 2=full)


LASER OPTIONS:

  -lsr-num, -no-lsr-num                   (default: off)

        Number of LSR backdoor variables.


  -lsr-verb, -no-lsr-verb                 (default: off)

        Debugging for MinLS tool.


  -replay-bd, -no-replay-bd               (default: off)

        Replay the decisions that witness the backdoor.



  -max-lsr-size = <int32>  [   1 .. imax] (default: 5)

        Maximum size LSR set to try.


  -min-lsr-size = <int32>  [   1 .. imax] (default: 1)

        Minimum size LSR set to try.



  -lsr-out    = <string>

        Write LSR backdoor to a file (zero-based).



MAIN OPTIONS:

  -pre, -no-pre                           (default: on)

        Completely turn on/off any preprocessing.


  -verb         = <int32>  [  -1 ..    2] (default: 1)

        Verbosity level (0=silent, 1=some, 2=more).

  -cpu-lim      = <int32>  [   0 .. imax] (default: 2147483647)

        Limit on CPU time allowed in seconds.


  -mem-lim      = <int32>  [   0 .. imax] (default: 2147483647)

        Limit on memory usage in megabytes.



  -assumptions = <string>

        If given, use the assumptions in the file.

  -var-sets-file = <string>

        List of potential backdoor sets to try (dimacs format without header).

  -dimacs     = <string>

        If given, stop after preprocessing and write the result to this file.

  -decision-vars = <string>

        Only branch on the listed vars (zero-based).

  -replay-bd-file = <string>

        Infile of decisions that witness the backdoor (one-based).

  -replay-bd-out = <string>

        Outfile of decisions that witness the backdoor (one-based).


SIMP OPTIONS:

  -asymm, -no-asymm                       (default: off)

        Shrink clauses by asymmetric branching.

  -rcheck, -no-rcheck                     (default: off)

        Check if a clause is already implied. (costly)

  -elim, -no-elim                         (default: on)

        Perform variable elimination.


  -simp-gc-frac = <double> (   0 ..  inf) (default: 0.5)

        The fraction of wasted memory allowed before a garbage collection is triggered during simplification.


  -cl-lim       = <int32>  [  -1 .. imax] (default: 20)

        Variables are not eliminated if it produces a resolvent with a length above this limit. -1 means no limit

  -sub-lim      = <int32>  [  -1 .. imax] (default: 1000)

        Do not check if subsumption against a clause larger than this. -1 means no limit.

  -grow         = <int32>  [imin .. imax] (default: 0)

        Allow a variable elimination step to grow by a number of clauses.


HELP OPTIONS:

  --help        Print help message.
  --help-verb   Print verbose help message.
```






