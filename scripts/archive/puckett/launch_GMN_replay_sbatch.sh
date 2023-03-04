#!/bin/bash

runnum=$1
maxsegments=$2

for ((i=0; i<=$2; i++))
do
    fnameout_pattern='/farm_out/puckett/gmn_replayed_'$runnum'_segment'$i'.out'
    sbatch --output=$fnameout_pattern run_GMN_swif2.sh $runnum -1 0 e1209019 $i 1
done

