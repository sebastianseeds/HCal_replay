#!/bin/bash

# swif2 job submit for simTOF.C
#Usage
#./submit-mcTOF-jobs.sh <kinematic> <field> <iteration>

kine=$1
field=$2
iter=$3

cp $SBS/run_replay_here/.rootrc $SWIF_JOB_WORK_DIR

work_flow='seeds_MCTOF_sbs'$kine'_field'$field
swif2 create $work_flow

jobname=$kine-$field-$iter

script='/w/halla-scshelf2102/sbs/seeds/HCal_replay/hcal/analysis/TOF/simulation/run_mcTOF.sh'

swif2 add-job -workflow $work_flow -partition production -name MCTOF-$jobname -cores 1 -disk 25GB -ram 2000MB $script $kine $field $iter

# run the workflow and then print status
swif2 run $work_flow
echo -e "\n Getting workflow status.. \n"
swif2 status $work_flow
