#!/bin/bash

maplesat="/home/ezulkosk/laser-maple-structure/laser/simp/maplesat"

base_dir=$1
name=$2
ccmin=$3 # ccmin2 or ccmin2_ar or ccmin2_nr
res=$4 # "" or -never-restart or -always-restart

if [ -s ${base_dir}/backbones/${name}.backbone ]
then
    echo "Backbone exists"
    bb="-backbone-file=${base_dir}/backbones/${name}.backbone -bb-metrics-out=${base_dir}/${ccmin}_bb_metrics/${name}.bb_metrics "
else
    echo "Backbone unfound"
    bb=""
fi

timeout 18020  $maplesat -cpu-lim=18000 -mem-lim=6600 -no-pre  -lsr-num -ccmin-mode=2 -verb=0 -lsr-out=${base_dir}/${ccmin}_confl_side_lsr/${name}.lsr -all-dec-out=${base_dir}/${ccmin}_all_decs/${name}.all_decs -conf-side-lsr -avg-clause-lsr-out=${base_dir}/${ccmin}_avg_lsr/${name}.avg_lsr -cmty-loc-out=${base_dir}/${ccmin}_cmty_loc/${name}.cmty_loc $res  -cmty-file=${base_dir}/cmty/${name}.cmty ${base_dir}/cnf/${name}.cnf $bb

