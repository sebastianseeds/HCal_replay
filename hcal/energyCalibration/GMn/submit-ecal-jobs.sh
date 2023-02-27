#!/bin/bash

# swif2 job submit for ecal.C for hcal calibration
# note that diagnostics should be done first in order to have good elastic cuts by kine/target/field-setting located in configuration files
#Usage
#./submit-ecal-jobs.sh <kinematic> <iteration>

kine=$1
iter=$4

cp $SBS/run_replay_here/.rootrc $SWIF_JOB_WORK_DIR

work_flow='seeds_GMn_ecal_sbs'$kine
swif2 create $work_flow

jobname=$kine-$iter

script='/w/halla-scshelf2102/sbs/seeds/HCal_replay/hcal/energyCalibration/GMn/run_ecal.sh'

swif2 add-job -workflow $work_flow -partition production -name ecal-$jobname -cores 1 -disk 25GB -ram 1500MB $script $kine $iter

# run the workflow and then print status
swif2 run $work_flow
echo -e "\n Getting workflow status.. \n"
swif2 status $work_flow
