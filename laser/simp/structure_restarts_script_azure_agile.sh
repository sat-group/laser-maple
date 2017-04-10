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
timeout 306 $maplesat -no-pre  -lsr-num -ccmin-mode=2 -verb=0 -lsr-out=${base_dir}/ccmin2_nr_confl_side_lsr/${name}.lsr -all-dec-out=${base_dir}/ccmin2_nr_all_decs/${name}.all_decs  -conf-side-lsr -avg-clause-lsr -never-restart -cmty-file=${base_dir}/cmty/${name}.cmty ${base_dir}/cnf/${name}.cnf  $bb
echo "--------------"
echo "Always restart"
echo "--------------"
timeout 306 $maplesat -no-pre  -lsr-num -ccmin-mode=2 -verb=0 -lsr-out=${base_dir}/ccmin2_ar_confl_side_lsr/${name}.lsr -all-dec-out=${base_dir}/ccmin2_ar_all_decs/${name}.all_decs -conf-side-lsr -avg-clause-lsr -always-restart -cmty-file=${base_dir}/cmty/${name}.cmty ${base_dir}/cnf/${name}.cnf $bb
echo "-----------------"
echo "Standard restarts"
echo "-----------------"
timeout 306 $maplesat -no-pre  -lsr-num -ccmin-mode=2 -verb=0 -lsr-out=${base_dir}/ccmin2_confl_side_lsr/${name}.lsr -all-dec-out=${base_dir}/ccmin2_all_decs/${name}.all_decs -conf-side-lsr -avg-clause-lsr  -cmty-file=${base_dir}/cmty/${name}.cmty ${base_dir}/cnf/${name}.cnf $bb
