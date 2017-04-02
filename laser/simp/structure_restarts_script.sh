#!/bin/bash

echo "-------------"
echo "Never restart"
echo "-------------"
./maplesat -lsr-out=lsr -all-dec-out=alldecs -lsr-num -ccmin-mode=2 -verb=0  -conf-side-lsr -avg-clause-lsr -never-restart -cmty-file=bench_2390.smt2.cmty bench_2390.smt2.cnf
echo "--------------"
echo "Always restart"
echo "--------------"
./maplesat -lsr-out=lsr -all-dec-out=alldecs -lsr-num -ccmin-mode=2 -verb=0  -conf-side-lsr -avg-clause-lsr -always-restart -cmty-file=bench_2390.smt2.cmty bench_2390.smt2.cnf
echo "-----------------"
echo "Standard restarts"
echo "-----------------"
./maplesat -lsr-out=lsr -all-dec-out=alldecs -lsr-num -ccmin-mode=2 -verb=0  -conf-side-lsr -avg-clause-lsr  -cmty-file=bench_2390.smt2.cmty bench_2390.smt2.cnf
