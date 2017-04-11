#!/bin/bash

maplesat="/home/ezulkosk/laser-maple-structure/laser/simp/maplesat"

base_dir=$1
name=$2

if [ -s ${base_dir}/backbones/${name}.backbone ]
then
    echo "Backbone exists"
    bb="-backbone-file=${base_dir}/backbones/${name}.backbone -bb-metrics-out=${base_dir}/ccmin2_bb_metrics/${name}.bb_metrics "
else
    echo "Backbone unfound"
    bb=""
fi



echo "-------------"
echo "Never restart"
echo "-------------"
#timeout 5100 $maplesat -no-pre  -lsr-num -ccmin-mode=2 -verb=0 -lsr-out=${base_dir}/ccmin2_nr_confl_side_lsr/${name}.lsr -all-dec-out=${base_dir}/ccmin2_nr_all_decs/${name}.all_decs  -conf-side-lsr -avg-clause-lsr -never-restart -cmty-file=${base_dir}/cmty/${name}.cmty ${base_dir}/cnf/${name}.cnf  $bb
echo "--------------"
echo "Always restart"
echo "--------------"
#timeout 5100 $maplesat -no-pre  -lsr-num -ccmin-mode=2 -verb=0 -lsr-out=${base_dir}/ccmin2_ar_confl_side_lsr/${name}.lsr -all-dec-out=${base_dir}/ccmin2_ar_all_decs/${name}.all_decs -conf-side-lsr -avg-clause-lsr -always-restart -cmty-file=${base_dir}/cmty/${name}.cmty ${base_dir}/cnf/${name}.cnf $bb
echo "-----------------"
echo "Standard restarts"
echo "-----------------"
timeout 18020  $maplesat -cpu-lim=18000 -mem-lim=13312 -no-pre  -lsr-num -ccmin-mode=2 -verb=0 -lsr-out=${base_dir}/ccmin2_confl_side_lsr/${name}.lsr -all-dec-out=${base_dir}/ccmin2_all_decs/${name}.all_decs -conf-side-lsr -avg-clause-lsr-out=${base_dir}/ccmin2_avg_lsr/${name}.avg_lsr -cmty-loc-out=${base_dir}/ccmin2_cmty_loc/${name}.cmty_loc  -cmty-file=${base_dir}/cmty/${name}.cmty ${base_dir}/cnf/${name}.cnf $bb
#timeout ${maplesat -cpu-lim=18000 -mem-lim=13312 -no-pre  -lsr-num -ccmin-mode=2  -lsr-out=./ccmin2_confl_side_lsr/uf100-0999.lsr -all-dec-out=./ccmin2_all_decs/uf100-0999.all_decs -conf-side-lsr -avg-clause-lsr-out=avglsrout  -cmty-file=./cmty/uf100-0999.cmty ./cnf/uf100-0999.cnf -backbone-file=backbones/uf100-0999.backbone -bb-metrics-out=bbout -cmty-loc-out=cmtyout
