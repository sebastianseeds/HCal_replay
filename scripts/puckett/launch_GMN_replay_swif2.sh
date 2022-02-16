#!/bin/bash

# single-run version for swif2:

runnum=$1
maxsegments=$2

#optional 3rd argument for output directory:
out_dir='/volatile/halla/sbs/puckett/GMN_REPLAYS'

if [ $# -eq 3 ];
then
    out_dir=$3
fi

for ((i=0; i<=$2; i++))
do
    fnameout_pattern='/farm_out/puckett/gmn_replayed_'$runnum'_segment'$i'.out'
    #    sbatch --output=$fnameout_pattern run_GMN_sbatch_nohodo.sh $runnum -1 0 e1209019 $i 1
    jobname='gmn_replay_'$runnum'_segment'$i
    
    # look for first segment on cache disk:
    firstsegname='e1209019_'$runnum'.evio.0.0'
    mssfirst='mss:/mss/halla/sbs/raw/'$firstsegname
    cachefirst='/cache/mss/halla/sbs/raw/'$firstsegname
    
    eviofilename='e1209019_'$runnum'.evio.0.'$i
    mssfilename='mss:/mss/halla/sbs/raw/'$eviofilename
    cachefile='/cache/mss/halla/sbs/raw/'$eviofilename
    
    script='/work/halla/sbs/puckett/GMN_ANALYSIS/run_GMN_swif2.sh'

    testfilename='/mss/halla/sbs/raw/'$eviofilename
    
    outfilename='match:e1209019_fullreplay_'$runnum'*seg'$i'*.root'
    logfilename='match:replay_gmn_'$runnum'*seg'$i'*.log'

    outcopydir=$out_dir'/rootfiles'
    logcopydir=$out_dir'/logs'

    if [ -f "$testfilename" ]; 
    then
	echo 'Adding new swif2 job, runnum='$runnum', segment='$i 
    
	if [ $i -gt 0 ]
	then
	    echo 'segment '$i' also requires first segment'
	    swif2 add-job -workflow puckett_GMN_analysis -partition production -name $jobname -cores 1 -disk 25GB -ram 1500MB -input $cachefile $mssfilename -input $cachefirst $mssfirst $script $runnum -1 0 e1209019 $i 1
	else
	    echo 'segment '$i' IS first segment'
	    swif2 add-job -workflow puckett_GMN_analysis -partition production -name $jobname -cores 1 -disk 25GB -ram 1500MB -input $cachefile $mssfilename $script $runnum -1 0 e1209019 $i 1
	fi
    fi
done
