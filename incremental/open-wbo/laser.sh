#!/bin/bash

rm -rf solvers/laser
mkdir solvers/laser
cp -r ../../laser/* solvers/laser

sed -i '/VERSION    =/c\VERSION    = core' Makefile
sed -i '/SOLVERNAME =/c\SOLVERNAME = "laser"' Makefile
sed -i '/SOLVERDIR  =/c\SOLVERDIR  = laser' Makefile
sed -i '/NSPACE     =/c\NSPACE     = Minisat' Makefile

make clean && make
