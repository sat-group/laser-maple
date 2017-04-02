#!/bin/bash

maplesat="./maplesat"
name=$1

bb="-backbone-file=${1}.backbone"

echo "-------------"
echo "Never restart"
echo "-------------"
$maplesat -lsr-out=lsr -all-dec-out=alldecs -lsr-num -ccmin-mode=2 -verb=0  -conf-side-lsr -avg-clause-lsr -never-restart -cmty-file=${name}.cmty ${name}.cnf  $bb
echo "--------------"
echo "Always restart"
echo "--------------"
$maplesat -lsr-out=lsr -all-dec-out=alldecs -lsr-num -ccmin-mode=2 -verb=0  -conf-side-lsr -avg-clause-lsr -always-restart -cmty-file=${name}.cmty ${name}.cnf $bb
echo "-----------------"
echo "Standard restarts"
echo "-----------------"
$maplesat -lsr-out=lsr -all-dec-out=alldecs -lsr-num -ccmin-mode=2 -verb=0  -conf-side-lsr -avg-clause-lsr  -cmty-file=${name}.cmty ${name}.cnf $bb
