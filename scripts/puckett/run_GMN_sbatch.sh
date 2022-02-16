#!/bin/bash

#SBATCH --partition=production
#SBATCH --account=halla
#SBATCH --chdir=/work/halla/sbs/puckett/GMN_ANALYSIS
#SBATCH --mem-per-cpu=2000
#SBATCH --constraint=farm19

source /site/12gev_phys/softenv.sh 2.5
#source /work/halla/sbs/puckett/install/bin/g4sbs.sh 
cd /work/halla/sbs/puckett/GMN_ANALYSIS

#root -b -q 'GEP_FOM_quick_and_dirty.C+("setup_GEP_FOM_option1_67cm.txt","GEP_FOM_option1_67cm_take2.root")' >output_GEP_FOM_option1_67cm_take2.txt

#root -b -q 'GEM_reconstruct_standalone_consolidated.C+("../UVA_EEL_DATA/gem_hit_2811*.root","config_UVA_EEL_5layer_Jan2021.txt", "temp2811_gainmatch.root")' >output2811_gainmatch.txt

# setup environment for ANALYZER and SBS-offline:

export ANALYZER=/work/halla/sbs/ANALYZER/install
source $ANALYZER/bin/setup.sh
source /work/halla/sbs/SBS_OFFLINE/install/bin/sbsenv.sh

export DB_DIR=$SBS_REPLAY/DB
export DATA_DIR=/cache/mss/halla/sbs/raw
export OUT_DIR=/volatile/halla/sbs/puckett/GMN_REPLAYS/SBS1_OPTICS/rootfiles
export LOG_DIR=/volatile/halla/sbs/puckett/GMN_REPLAYS/SBS1_OPTICS/logs

runnum=$1
maxevents=$2
firstevent=$3

prefix=$4
firstsegment=$5
maxsegments=$6




#cmd="analyzer -b -q 'replay/replay_BBGEM.C+("$runnum","$firstsegment","$maxsegments")'"

#echo $cmd


analyzer -b -q 'replay_gmn.C+('$runnum','$maxevents','$firstevent','\"$prefix\"','$firstsegment','$maxsegments')'


