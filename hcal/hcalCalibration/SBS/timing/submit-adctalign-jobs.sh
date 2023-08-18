#!/bin/bash

# swif2 job submit for adct_align.C for hcal calibration
# note that diagnostics should be done first in order to have good elastic cuts by kine/target/field-setting located in configuration files
#Usage
#./submit-adctalign-jobs.sh <string: experiment> <int: configuration> <bool: quasi-replay> <int: replay-pass> <bool: hydrogen-only calibration>

exper=$1
kine=$2
iter=$3
pass=$4
lh2only=$5

cp $SBS/run_replay_here/.rootrc $SWIF_JOB_WORK_DIR

work_flow='adctalign_'$exper'_sbs'$kine'_i'$iter'_pass'$pass'_lh2o'$lh2only

swif2 create $work_flow

jobname=$kine-$iter

script='/w/halla-scshelf2102/sbs/seeds/HCal_replay/hcal/hcalCalibration/SBS/timing/run_adct_align.sh'

swif2 add-job -workflow $work_flow -partition production -name adctalign-$jobname -cores 1 -disk 25GB -ram 2000MB $script $exper $kine $iter $pass $lh2only

# run the workflow and then print status
swif2 run $work_flow
echo -e "\n Getting workflow status.. \n"
swif2 status $work_flow
