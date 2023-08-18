#!/bin/bash

# swif2 job submit for ecal.C for hcal calibration
# note that diagnostics should be done first in order to have good elastic cuts by kine/target/field-setting located in configuration files
#Usage
#./submit-qreplay-jobs.sh <experiment> <configuration> <replay-pass> <target-option> <selected-energy-timestamp> <replacement-experiment> <replacement-configuration> <replacement-pass> <replacement-database-timestamp> <verbose-option>

exper=$1
config=$2
pass=$3
sts=$4
rexper=$5
rconfig=$6
rpass=$7
qts=$8
h2only=$9
verbose=false

cp $SBS/run_replay_here/.rootrc $SWIF_JOB_WORK_DIR

work_flow='qreplay_'$exper'_from'$rconfig'_to'$config
swif2 create $work_flow

jobname=$rconfig-$config

script='/w/halla-scshelf2102/sbs/seeds/HCal_replay/hcal/hcalCalibration/SBS/qreplay/run_qreplay.sh'

swif2 add-job -workflow $work_flow -partition production -name qreplay-$jobname -cores 1 -disk 25GB -ram 2000MB $script $exper $config $pass $sts $rexper $rconfig $rpass $qts $h2only $verbose

# run the workflow and then print status
swif2 run $work_flow
echo -e "\n Getting workflow status.. \n"
swif2 status $work_flow
