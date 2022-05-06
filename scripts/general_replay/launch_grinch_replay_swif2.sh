#!/bin/bash

# single-run version for swif2:
# usage:
# verify that a workflow exists with name msatnik_GMN_analysis
## if not: swif2 create -workflow msatnik_GMN_analysis
### to check on workflow status: swif2 status msatnik_GMN_analysis
# from working directory with script:
# ./launch_HCal_replay_swif2.sh <run number> <last segment>
# after script completes:
# swif2 run msatnik_GMN_analysis
# see https://scicomp.jlab.org/scicomp/swif/active for more information

runnum=$1
maxsegments=$2

#optional 3rd argument for output directory:
out_dir='/volatile/halla/sbs/msatnik/GMN_Replays'

if [ $# -eq 3 ];
then
    out_dir=$3
fi

for ((i=0; i<=$2; i++))
do
    fnameout_pattern='/farm_out/msatnik/replay_grinch_'$runnum'_segment'$i'.out'
    #    sbatch --output=$fnameout_pattern run_GMN_sbatch_nohodo.sh $runnum -1 0 e1209019 $i 1
    jobname='replay_grinch_'$runnum'_segment'$i
    
    # look for first segment on cache disk:
    firstsegname='e1209019_'$runnum'.evio.0.0'
    mssfirst='mss:/mss/halla/sbs/raw/'$firstsegname
    cachefirst='/cache/mss/halla/sbs/raw/'$firstsegname
    
    eviofilename='e1209019_'$runnum'.evio.0.'$i
    mssfilename='mss:/mss/halla/sbs/raw/'$eviofilename
    cachefile='/cache/mss/halla/sbs/raw/'$eviofilename
    
    script='/w/halla-scshelf2102/sbs/msatnik/run_grinch_swif2.sh'
    testfilename='/mss/halla/sbs/raw/'$eviofilename
    
    outfilename='match:e1209019_gmn_'$runnum'*seg'$i'*.root'
    logfilename='match:replay_gmn_'$runnum'*seg'$i'*.log'

    outcopydir=$out_dir
    logcopydir=$out_dir'/logs'

    if [ -f "$testfilename" ]; 
    then
	echo 'Adding new swif2 job, runnum='$runnum', segment='$i 
    
	if [ $i -gt 0 ]
	then
	    echo 'segment '$i' also requires first segment'
	    swif2 add-job -workflow msatnik_GMN_analysis -partition production -name $jobname -cores 1 -disk 25GB -ram 1500MB -input $cachefile $mssfilename -input $cachefirst $mssfirst $script $runnum -1 0 e1209019 $i 1
	else
	    echo 'segment '$i' IS first segment'
	    swif2 add-job -workflow msatnik_GMN_analysis -partition production -name $jobname -cores 1 -disk 25GB -ram 1500MB -input $cachefile $mssfilename $script $runnum -1 0 e1209019 $i 1
	fi
    fi
done
