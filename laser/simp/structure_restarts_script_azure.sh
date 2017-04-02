#!/bin/bash

maplesat="/home/ezulkosk/laser-maple-structure/laser/simp/maplesat"

base_dir=$1
name=$2

if [ -s ${base_dir}/backbones/${name}.backbone ]
then
    echo "Backbone exists"
    bb="-backbone-file=${base_dir}/backbones/${name}.backbone"
else
    echo "Backbone unfound"
    bb=""
fi



echo "-------------"
echo "Never restart"
echo "-------------"
timeout 300 $maplesat -lsr-out=lsr -all-dec-out=alldecs -lsr-num -ccmin-mode=2 -verb=0  -conf-side-lsr -avg-clause-lsr -never-restart -cmty-file=${base_dir}/cmty/${name}.cmty ${base_dir}/cnf/${name}.cnf  $bb
echo "--------------"
echo "Always restart"
echo "--------------"
timeout 300 $maplesat -lsr-out=lsr -all-dec-out=alldecs -lsr-num -ccmin-mode=2 -verb=0  -conf-side-lsr -avg-clause-lsr -always-restart -cmty-file=${base_dir}/cmty/${name}.cmty ${base_dir}/cnf/${name}.cnf $bb
echo "-----------------"
echo "Standard restarts"
echo "-----------------"
timeout 300 $maplesat -lsr-out=lsr -all-dec-out=alldecs -lsr-num -ccmin-mode=2 -verb=0  -conf-side-lsr -avg-clause-lsr  -cmty-file=${base_dir}/cmty/${name}.cmty ${base_dir}/cnf/${name}.cnf $bb
