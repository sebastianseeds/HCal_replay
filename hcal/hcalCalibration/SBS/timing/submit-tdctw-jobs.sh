#!/bin/bash

# swif2 job submit for tdc_tw.C for hcal calibration
# note that diagnostics should be done first in order to have good elastic cuts by kine/target/field-setting located in configuration files
#Usage
#./submit-tdctw-jobs.sh <string: experiment> <int: configuration> <bool: quasi-replay> <int: replay pass> <int: number of energy calibrations> <bool: parameter override option> <string: parameter override timestamp>

exper=$1
kine=$2
iter=$3
pass=$4
nset=$5
pover=$6
pts=$7

cp $SBS/run_replay_here/.rootrc $SWIF_JOB_WORK_DIR

work_flow='tdctw_'$exper'_sbs'$kine'_i'$iter'_pass'$pass'_nset'$nset

swif2 create $work_flow

jobname=$kine-$iter

script='/w/halla-scshelf2102/sbs/seeds/HCal_replay/hcal/hcalCalibration/SBS/timing/run_tdc_tw.sh'

swif2 add-job -workflow $work_flow -partition production -name tdctw-$jobname -cores 1 -disk 25GB -ram 2000MB $script $exper $kine $iter $pass $nset $pover $pts

# run the workflow and then print status
swif2 run $work_flow
echo -e "\n Getting workflow status.. \n"
swif2 status $work_flow
