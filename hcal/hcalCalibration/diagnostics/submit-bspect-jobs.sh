#!/bin/bash

# swif2 job submit for ecal_blockE_spectra.C for determination of sig_E to E ratio
#Usage
#./submit_bspect_jobs.sh <kinematic>

kine=$1

cp $SBS/run_replay_here/.rootrc $SWIF_JOB_WORK_DIR

work_flow='seeds_bspect_sbs'$kine
swif2 create $work_flow

jobname=$kine

script='/w/halla-scshelf2102/sbs/seeds/HCal_replay/hcal/hcalCalibration/diagnostics/run_bspect.sh'

swif2 add-job -workflow $work_flow -partition production -name ecalbspect-$jobname -cores 1 -disk 25GB -ram 2000MB $script $kine

# run the workflow and then print status
swif2 run $work_flow
echo -e "\n Getting workflow status.. \n"
swif2 status $work_flow
