#!/bin/bash

#SBATCH --partition=production
#SBATCH --account=halla
#SBATCH --mem-per-cpu=1500

#cd /work/halla/sbs/puckett/GMN_ANALYSIS

echo 'working directory ='
echo $PWD

echo 'swif_job_work_dir='$SWIF_JOB_WORK_DIR

# source login stuff since swif2 completely nukes any sensible default software environment
#source /home/puckett/.cshrc
#source /home/puckett/.login

#echo 'before sourcing environment, env = '
#env 

#module load gcc/9.2.0

#!/bin/bash 

MODULES=/etc/profile.d/modules.sh 

if [[ $(type -t module) != function && -r ${MODULES} ]]; then 
source ${MODULES} 
fi 

if [ -d /apps/modulefiles ]; then 
module use /apps/modulefiles 
fi 

module load gcc/9.2.0 

source /site/12gev_phys/softenv.sh 2.5

ldd /work/halla/sbs/ANALYZER/install/bin/analyzer |& grep not

#source /etc/profile.d/modules.sh

#module load gcc/9.2.0


#export ROOTSYS=/site/12gev_phys/2.5/Linux_CentOS7.7.1908-gcc9.2.0/root/6.24.06

#source $ROOTSYS/bin/thisroot.sh

#echo 'after attempting to source environment, env = '
#env
#source /work/halla/sbs/puckett/install/bin/g4sbs.sh 
#cd /work/halla/sbs/puckett/GMN_ANALYSIS

#cd /scratch/slurm/puckett

#root -b -q 'GEP_FOM_quick_and_dirty.C+("setup_GEP_FOM_option1_67cm.txt","GEP_FOM_option1_67cm_take2.root")' >output_GEP_FOM_option1_67cm_take2.txt

#root -b -q 'GEM_reconstruct_standalone_consolidated.C+("../UVA_EEL_DATA/gem_hit_2811*.root","config_UVA_EEL_5layer_Jan2021.txt", "temp2811_gainmatch.root")' >output2811_gainmatch.txt

# setup environment for ANALYZER and SBS-offline:

echo 'working directory = '$PWD

export ANALYZER=/work/halla/sbs/ANALYZER/install
source $ANALYZER/bin/setup.sh
source /work/halla/sbs/SBS_OFFLINE/install/bin/sbsenv.sh

#cp $SBS/run_replay_here/.rootrc $PWD

#mkdir -p $PWD/in
#mkdir -p $PWD/rootfiles
#mkdir -p $PWD/logs

export SBS_REPLAY=/work/halla/sbs/SBS_REPLAY/SBS-replay
export DB_DIR=$SBS_REPLAY/DB
export DATA_DIR=/cache/mss/halla/sbs/raw

export OUT_DIR=$SWIF_JOB_WORK_DIR
export LOG_DIR=$SWIF_JOB_WORK_DIR

echo 'OUT_DIR='$OUT_DIR
echo 'LOG_DIR='$LOG_DIR
# try this under swif2:
#export DATA_DIR=$PWD
#export OUT_DIR=/volatile/halla/sbs/puckett/GMN_REPLAYS/SBS1_OPTICS/rootfiles
#export LOG_DIR=/volatile/halla/sbs/puckett/GMN_REPLAYS/SBS1_OPTICS/logs
#mkdir -p /volatile/halla/sbs/puckett/GMN_REPLAYS/rootfiles
#mkdir -p /volatile/halla/sbs/puckett/GMN_REPLAYS/logs

#export OUT_DIR=/volatile/halla/sbs/puckett/GMN_REPLAYS/SBS4/rootfiles
#export LOG_DIR=/volatile/halla/sbs/puckett/GMN_REPLAYS/SBS4/rootfiles
export ANALYZER_CONFIGPATH=$SBS_REPLAY/replay

runnum=$1
maxevents=$2
firstevent=$3

prefix=$4
firstsegment=$5
maxsegments=$6

cp $SBS/run_replay_here/.rootrc $SWIF_JOB_WORK_DIR

#cmd="analyzer -b -q 'replay/replay_BBGEM.C+("$runnum","$firstsegment","$maxsegments")'"

#echo $cmd


analyzer -b -q 'replay_gmn.C+('$runnum','$maxevents','$firstevent','\"$prefix\"','$firstsegment','$maxsegments')'

outfilename=$OUT_DIR'/e1209019_*'$runnum'*.root'
logfilename=$LOG_DIR'/replay_gmn_'$runnum'*.log' 

mv $outfilename /volatile/halla/sbs/puckett/GMN_REPLAYS/rootfiles
mv $logfilename /volatile/halla/sbs/puckett/GMN_REPLAYS/logs
