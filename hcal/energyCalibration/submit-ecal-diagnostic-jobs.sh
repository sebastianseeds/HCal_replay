#!/bin/bash

# swif2 job submit for ecal_diagnostic.C on road to hcal calibration
#Usage
#./submit-ecal-diagnostic-jobs.sh <kinematic> <target (lh2/ld2)> <field> <iteration>

kine=$1
tar=$2
field=$3
iter=$4

cp $SBS/run_replay_here/.rootrc $SWIF_JOB_WORK_DIR

work_flow='seeds_ecalDiagnostic_'$tar'_sbs'$kine'_f'$field
swif2 create $work_flow

jobname=$kine-$tar-$field-$iter

script='/w/halla-scshelf2102/sbs/seeds/HCal_replay/hcal/energyCalibration/run_diagnostic.sh'

swif2 add-job -workflow $work_flow -partition production -name ecaldiag-$jobname -cores 1 -disk 25GB -ram 1500MB $script $kine $tar $field $iter

# run the workflow and then print status
swif2 run $work_flow
echo -e "\n Getting workflow status.. \n"
swif2 status $work_flow
