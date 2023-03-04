#!/bin/csh 

#SBATCH --partition=production
#SBATCH --account=halla
#SBATCH --mem-per-cpu=2000
#SBATCH --constraint=farm19

#cd /work/halla/sbs/puckett/GMN_ANALYSIS

echo 'working directory ='
echo $PWD

# source login stuff since swif2 completely nukes any sensible default software environment
#source /home/puckett/.cshrc
#source /home/puckett/.login

#echo 'before sourcing environment, env = '
#env 

source /site/12gev_phys/softenv.csh 2.5

#echo 'after attempting to source environment, env = '
#env
#source /work/halla/sbs/puckett/install/bin/g4sbs.sh 
#cd /work/halla/sbs/puckett/GMN_ANALYSIS

#cd /scratch/slurm/puckett

#root -b -q 'GEP_FOM_quick_and_dirty.C+("setup_GEP_FOM_option1_67cm.txt","GEP_FOM_option1_67cm_take2.root")' >output_GEP_FOM_option1_67cm_take2.txt

#root -b -q 'GEM_reconstruct_standalone_consolidated.C+("../UVA_EEL_DATA/gem_hit_2811*.root","config_UVA_EEL_5layer_Jan2021.txt", "temp2811_gainmatch.root")' >output2811_gainmatch.txt

# setup environment for ANALYZER and SBS-offline:

echo 'working directory = '$PWD

setenv ANALYZER /work/halla/sbs/ANALYZER/install
source $ANALYZER/bin/setup.csh
source /work/halla/sbs/SBS_OFFLINE/install/bin/sbsenv.csh

#cp $SBS/run_replay_here/.rootrc $PWD

#mkdir -p $PWD/in
#mkdir -p $PWD/rootfiles
#mkdir -p $PWD/logs

setenv SBS_REPLAY /work/halla/sbs/SBS_REPLAY/SBS-replay
setenv DB_DIR $SBS_REPLAY/DB
setenv DATA_DIR /cache/mss/halla/sbs/raw
# try this under swif2:
#export DATA_DIR=$PWD
#export OUT_DIR=/volatile/halla/sbs/puckett/GMN_REPLAYS/SBS1_OPTICS/rootfiles
#export LOG_DIR=/volatile/halla/sbs/puckett/GMN_REPLAYS/SBS1_OPTICS/logs
mkdir -p /volatile/halla/sbs/puckett/GMN_REPLAYS/rootfiles
mkdir -p /volatile/halla/sbs/puckett/GMN_REPLAYS/logs

setenv OUT_DIR /volatile/halla/sbs/puckett/GMN_REPLAYS/rootfiles
setenv LOG_DIR /volatile/halla/sbs/puckett/GMN_REPLAYS/logs
setenv ANALYZER_CONFIGPATH $SBS_REPLAY/replay

set runnum=$1
set maxevents=$2
set firstevent=$3

set prefix=$4
set firstsegment=$5
set maxsegments=$6

#cmd="analyzer -b -q 'replay/replay_BBGEM.C+("$runnum","$firstsegment","$maxsegments")'"

#echo $cmd


analyzer -b -q 'replay_gmn.C+('$runnum','$maxevents','$firstevent','\"$prefix\"','$firstsegment','$maxsegments')'


