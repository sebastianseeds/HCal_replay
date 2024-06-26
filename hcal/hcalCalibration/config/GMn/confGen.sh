#!/bin/sh

# SSeeds - 2.27.23 - produce configuration files by pass/kinematic/target/field for both ecal.C and tcal.C

## Usage
# Simple script takes data from mss and places it into the cache 
#./confGen.sh <pass> <kinematic>

echo 'This script produces .cfg files for calibration scripts.'
echo 'Which pass would you like to work with?'
read pass
echo 'Which kinematic?'
read kine

#Set input limits
if [ $pass -ne 0 ] && [ $pass -ne 1 ]
then
echo 'This script only configured for pass0 and pass1 as of 2.23.23. Please enter 0 or 1.'
exit
fi

if [ $kine -ne 4 ] && [ $kine -ne 7 ] && [ $kine -ne 11 ] && [ $kine -ne 14 ] && [ $kine -ne 8 ] && [ $kine -ne 9 ]
then
echo 'GMn kinematics available: 4, 7, 11, 14, 8, and 9. Please enter a valid kinematic.'
exit
fi

#Declare general cut parameters over all kinematics
#Order (by idx) here: kine = (4 7 11 14 8 9)
globalcut=('bb.tr.n==1&&bb.ps.e>0.2&&abs(bb.tr.vz[0])<0.08&&bb.gem.track.nhits>3&&abs(bb.etot_over_p-0.92)<0.2&&sbs.hcal.e>0.01&&bb.ps.e+bb.sh.e>1.7' 'bb.tr.n>0&&abs(bb.tr.vz[0])<0.08&&bb.gem.track.nhits>2&&bb.tr.p[0]>2.0&&sbs.hcal.e>0.01&&bb.ps.e>0.2' 'bb.tr.n>0&&abs(bb.tr.vz[0])<0.08&&bb.gem.track.nhits>2&&bb.tr.p[0]>2.0&&sbs.hcal.e>0.01&&bb.ps.e>0.2' 'bb.tr.n>0&&abs(bb.tr.vz[0])<0.08&&bb.gem.track.nhits>2&&bb.tr.p[0]>1.6&&sbs.hcal.e>0.01&&bb.ps.e>0.2' 'bb.tr.n==1&&bb.ps.e>0.2&&abs(bb.tr.vz[0])<0.08&&bb.gem.track.nhits>3&&sbs.hcal.e>0.01&&abs(bb.tr.tg_th[0])<0.15&&abs(bb.tr.tg_ph[0])<0.3' 'bb.tr.n==1&&bb.ps.e>0.2&&abs(bb.tr.vz[0])<0.08&&bb.gem.track.nhits>3&&sbs.hcal.e>0.01&&abs(bb.tr.tg_th[0])<0.15&&abs(bb.tr.tg_ph[0])<0.3')
E_e=(3.7278 7.9072 9.8594 5.9649 5.9648 4.0148)
HCal_d=(11. 14. 14.5 14. 11. 11.)
HCal_th=(31.9 16.1 13.3 17.3 29.4 22.)
BB_th=(36. 40. 42. 46.5 26.5 49.)
W2_mean=(0.917994 0.88 0.925 0.870819 0.91 0.91)
W2_sig=(0.167922 0.5 0.325 0.19443 0.21 0.17)
#dy0=(-0.0270143 0.0158817 -0.0119 0.00220912 -0.0119 -0.0175041)
#dy_sig=(0.0869443 0.155968 0.08675 0.0824202 0.08675 0.134702)
atime0=(51.466 63.1685 50.36 58.7825 50.36 50.6158)
atime_sig=(3.67744 3.69617 3.73 3.68101 3.73 3.52261)
useAlshield=(0 0 0 0 0 0)

#Get the path to files
######################
if [ $pass -eq 0 ] || [ $pass -eq 1 ]
then
    if [ $kine -eq 4 ]
    then
    lh2path='/w/halla-scshelf2102/sbs/sbs-gmn/pass0/SBS4/LH2/rootfiles/'
    ld2path='/w/halla-scshelf2102/sbs/sbs-gmn/pass0/SBS4/LD2/rootfiles/'
    fi

    if [ $kine -eq 7 ]
    then
    lh2path='/w/halla-scshelf2102/sbs/sbs-gmn/pass0/SBS7/LH2/rootfiles/'
    ld2path='/w/halla-scshelf2102/sbs/sbs-gmn/pass0/SBS7/LD2/rootfiles/'
    fi

    if [ $kine -eq 11 ]
    then
    lh2path='/w/halla-scshelf2102/sbs/sbs-gmn/pass1/SBS11/LH2/rootfiles/'
    ld2path='/w/halla-scshelf2102/sbs/sbs-gmn/pass1/SBS11/LD2/rootfiles/'
    fi

    if [ $kine -eq 14 ]
    then
    lh2path='/w/halla-scshelf2102/sbs/sbs-gmn/pass1/SBS14/LH2/rootfiles/'
    ld2path='/w/halla-scshelf2102/sbs/sbs-gmn/pass1/SBS14/LD2/rootfiles/'
    fi

    if [ $kine -eq 8 ]
    then
    lh2path='/w/halla-scshelf2102/sbs/sbs-gmn/pass1/SBS8/LH2/rootfiles/'
    ld2path='/w/halla-scshelf2102/sbs/sbs-gmn/pass1/SBS8/LD2/rootfiles/'
    fi

    if [ $kine -eq 9 ]
    then
    lh2path='/w/halla-scshelf2102/sbs/sbs-gmn/pass1/SBS9/LH2/rootfiles/'
    ld2path='/w/halla-scshelf2102/sbs/sbs-gmn/pass1/SBS9/LD2/rootfiles/'
    fi
fi

#Generate full file extensions for desired kinematic
####################################################
#SBS4
if [ $kine -eq 4 ]
then

    #Index for general cut params
    idx=0
    #LH2 runs
    tar='lh2'
    #SBS field 0%
    field=0
    dx0_n=0
    dx0_p=0
    dx_sig_n=1
    dx_sig_p=1
    dy0=-0.0270143
    dy_sig=0.0869443
    configname='secal_'$tar'_sbs'$kine'_f'$field'.cfg'
    echo '#sseeds 2.23.23' >> $configname
    echo '#Config file for sbs'$kine' generated from pass'$pass >> $configname
    for a in {11573,11587,11588};
    do
    ls $lh2path*$a* | sort -V >> $configname
    done
    echo 'endlist' >> $configname
    echo ${globalcut[$idx]} >> $configname
    echo 'endcut' >> $configname
    echo 'E_e '${E_e[$idx]} >> $configname
    echo 'HCal_d '${HCal_d[$idx]} >> $configname
    echo 'HCal_th '${HCal_th[$idx]} >> $configname
    echo 'BB_th '${BB_th[$idx]} >> $configname
    echo 'W2_mean '${W2_mean[$idx]} >> $configname
    echo 'W2_sig '${W2_sig[$idx]} >> $configname
    echo 'dx0_n '$dx0_n >> $configname
    echo 'dx0_p '$dx0_p >> $configname
    echo 'dy0 '$dy0 >> $configname
    echo 'dx_sig_n '$dx_sig_n >> $configname
    echo 'dx_sig_p '$dx_sig_p >> $configname
    echo 'dy_sig '$dy_sig >> $configname
    echo 'atime0 '${atime0[$idx]} >> $configname
    echo 'atime_sig '${atime_sig[$idx]} >> $configname
    echo 'useAlshield '${useAlshield[$idx]} >> $configname
    mv $configname 'SBS'$kine'/'

    #SBS field 30%
    field=30
    dx0_n=0
    dx0_p=0
    dx_sig_n=1
    dx_sig_p=1
    dy0=-0.0270143
    dy_sig=0.0869443
    configname='secal_'$tar'_sbs'$kine'_f'$field'.cfg'
    echo '#sseeds 2.23.23' >> $configname
    echo '#Config file for sbs'$kine' generated from pass'$pass >> $configname
    for a in {11436,11500,11547,11548,11616};
    do
    ls $lh2path*$a* | sort -V >> $configname
    done
    echo 'endlist' >> $configname
    echo ${globalcut[$idx]} >> $configname
    echo 'endcut' >> $configname
    echo 'E_e '${E_e[$idx]} >> $configname
    echo 'HCal_d '${HCal_d[$idx]} >> $configname
    echo 'HCal_th '${HCal_th[$idx]} >> $configname
    echo 'BB_th '${BB_th[$idx]} >> $configname
    echo 'W2_mean '${W2_mean[$idx]} >> $configname
    echo 'W2_sig '${W2_sig[$idx]} >> $configname
    echo 'dx0_n '$dx0_n >> $configname
    echo 'dx0_p '$dx0_p >> $configname
    echo 'dy0 '$dy0 >> $configname
    echo 'dx_sig_n '$dx_sig_n >> $configname
    echo 'dx_sig_p '$dx_sig_p >> $configname
    echo 'dy_sig '$dy_sig >> $configname
    echo 'atime0 '${atime0[$idx]} >> $configname
    echo 'atime_sig '${atime_sig[$idx]} >> $configname
    echo 'useAlshield '${useAlshield[$idx]} >> $configname
    mv $configname 'SBS'$kine'/'

    #SBS field 50%
    field=50
    dx0_n=0
    dx0_p=0
    dx_sig_n=1
    dx_sig_p=1
    dy0=-0.0270143
    dy_sig=0.0869443
    configname='secal_'$tar'_sbs'$kine'_f'$field'.cfg'
    echo '#sseeds 2.23.23' >> $configname
    echo '#Config file for sbs'$kine' generated from pass'$pass >> $configname
    for a in {11589,11590,11592};
    do
    ls $lh2path*$a* | sort -V >> $configname
    done
    echo 'endlist' >> $configname
    echo ${globalcut[$idx]} >> $configname
    echo 'endcut' >> $configname
    echo 'E_e '${E_e[$idx]} >> $configname
    echo 'HCal_d '${HCal_d[$idx]} >> $configname
    echo 'HCal_th '${HCal_th[$idx]} >> $configname
    echo 'BB_th '${BB_th[$idx]} >> $configname
    echo 'W2_mean '${W2_mean[$idx]} >> $configname
    echo 'W2_sig '${W2_sig[$idx]} >> $configname
    echo 'dx0_n '$dx0_n >> $configname
    echo 'dx0_p '$dx0_p >> $configname
    echo 'dy0 '$dy0 >> $configname
    echo 'dx_sig_n '$dx_sig_n >> $configname
    echo 'dx_sig_p '$dx_sig_p >> $configname
    echo 'dy_sig '$dy_sig >> $configname
    echo 'atime0 '${atime0[$idx]} >> $configname
    echo 'atime_sig '${atime_sig[$idx]} >> $configname
    echo 'useAlshield '${useAlshield[$idx]} >> $configname
    mv $configname 'SBS'$kine'/'

    #LD2 runs
    tar='ld2'
    #SBS field 0%
    field=0
    dx0_n=0
    dx0_p=0
    dx_sig_n=1
    dx_sig_p=1
    dy0=-0.0270143
    dy_sig=0.0869443
    configname='secal_'$tar'_sbs'$kine'_f'$field'.cfg'
    echo '#sseeds 2.23.23' >> $configname
    echo '#Config file for sbs'$kine' generated from pass'$pass >> $configname
    for a in {11579,11580,11581,11582,11583,11586};
    do
    ls $ld2path*$a* | sort -V >> $configname
    done
    echo 'endlist' >> $configname
    echo ${globalcut[$idx]} >> $configname
    echo 'endcut' >> $configname
    echo 'E_e '${E_e[$idx]} >> $configname
    echo 'HCal_d '${HCal_d[$idx]} >> $configname
    echo 'HCal_th '${HCal_th[$idx]} >> $configname
    echo 'BB_th '${BB_th[$idx]} >> $configname
    echo 'W2_mean '${W2_mean[$idx]} >> $configname
    echo 'W2_sig '${W2_sig[$idx]} >> $configname
    echo 'dx0_n '$dx0_n >> $configname
    echo 'dx0_p '$dx0_p >> $configname
    echo 'dy0 '$dy0 >> $configname
    echo 'dx_sig_n '$dx_sig_n >> $configname
    echo 'dx_sig_p '$dx_sig_p >> $configname
    echo 'dy_sig '$dy_sig >> $configname
    echo 'atime0 '${atime0[$idx]} >> $configname
    echo 'atime_sig '${atime_sig[$idx]} >> $configname
    echo 'useAlshield '${useAlshield[$idx]} >> $configname
    mv $configname 'SBS'$kine'/'

    #SBS field 30%
    field=30
    dx0_n=0
    dx0_p=0
    dx_sig_n=1
    dx_sig_p=1
    dy0=-0.0270143
    dy_sig=0.0869443
    configname='secal_'$tar'_sbs'$kine'_f'$field'.cfg'
    echo '#sseeds 2.23.23' >> $configname
    echo '#Config file for sbs'$kine' generated from pass'$pass >> $configname
    for a in {11449,11451,11452,11456,11493,11494,11495,11496,11551,11554,11562,11563,11564,11565,11568,11570,11571};
    do
    ls $ld2path*$a* | sort -V >> $configname
    done
    echo 'endlist' >> $configname
    echo ${globalcut[$idx]} >> $configname
    echo 'endcut' >> $configname
    echo 'E_e '${E_e[$idx]} >> $configname
    echo 'HCal_d '${HCal_d[$idx]} >> $configname
    echo 'HCal_th '${HCal_th[$idx]} >> $configname
    echo 'BB_th '${BB_th[$idx]} >> $configname
    echo 'W2_mean '${W2_mean[$idx]} >> $configname
    echo 'W2_sig '${W2_sig[$idx]} >> $configname
    echo 'dx0_n '$dx0_n >> $configname
    echo 'dx0_p '$dx0_p >> $configname
    echo 'dy0 '$dy0 >> $configname
    echo 'dx_sig_n '$dx_sig_n >> $configname
    echo 'dx_sig_p '$dx_sig_p >> $configname
    echo 'dy_sig '$dy_sig >> $configname
    echo 'atime0 '${atime0[$idx]} >> $configname
    echo 'atime_sig '${atime_sig[$idx]} >> $configname
    echo 'useAlshield '${useAlshield[$idx]} >> $configname
    mv $configname 'SBS'$kine'/'

    #SBS field 50%
    field=50
    dx0_n=0
    dx0_p=0
    dx_sig_n=1
    dx_sig_p=1
    dy0=-0.0270143
    dy_sig=0.0869443
    configname='secal_'$tar'_sbs'$kine'_f'$field'.cfg'
    echo '#sseeds 2.23.23' >> $configname
    echo '#Config file for sbs'$kine' generated from pass'$pass >> $configname
    for a in {11593,11595};
    do
    ls $ld2path*$a* | sort -V >> $configname
    done
    echo 'endlist' >> $configname
    echo ${globalcut[$idx]} >> $configname
    echo 'endcut' >> $configname
    echo 'E_e '${E_e[$idx]} >> $configname
    echo 'HCal_d '${HCal_d[$idx]} >> $configname
    echo 'HCal_th '${HCal_th[$idx]} >> $configname
    echo 'BB_th '${BB_th[$idx]} >> $configname
    echo 'W2_mean '${W2_mean[$idx]} >> $configname
    echo 'W2_sig '${W2_sig[$idx]} >> $configname
    echo 'dx0_n '$dx0_n >> $configname
    echo 'dx0_p '$dx0_p >> $configname
    echo 'dy0 '$dy0 >> $configname
    echo 'dx_sig_n '$dx_sig_n >> $configname
    echo 'dx_sig_p '$dx_sig_p >> $configname
    echo 'dy_sig '$dy_sig >> $configname
    echo 'atime0 '${atime0[$idx]} >> $configname
    echo 'atime_sig '${atime_sig[$idx]} >> $configname
    echo 'useAlshield '${useAlshield[$idx]} >> $configname
    mv $configname 'SBS'$kine'/'

fi

#SBS7
if [ $kine -eq 7 ]
then

    #Index for general cut params
    idx=1
    #LH2 runs
    tar='lh2'
    #SBS field 85%
    field=85
    dx0_n=0
    dx0_p=0
    dx_sig_n=1
    dx_sig_p=1
    dy0=0.0158817
    dy_sig=0.155968
    configname='secal_'$tar'_sbs'$kine'_f'$field'.cfg'
    echo '#sseeds 2.23.23' >> $configname
    echo '#Config file for sbs'$kine' generated from pass'$pass >> $configname
    for a in {11989,11990,11991,11992,11993,11994,12000,12008,12022,12035,12049,12052,12053,12058,12064,12072};
    do
    ls $lh2path*$a* | sort -V >> $configname
    done
    echo 'endlist' >> $configname
    echo ${globalcut[$idx]} >> $configname
    echo 'endcut' >> $configname
    echo 'E_e '${E_e[$idx]} >> $configname
    echo 'HCal_d '${HCal_d[$idx]} >> $configname
    echo 'HCal_th '${HCal_th[$idx]} >> $configname
    echo 'BB_th '${BB_th[$idx]} >> $configname
    echo 'W2_mean '${W2_mean[$idx]} >> $configname
    echo 'W2_sig '${W2_sig[$idx]} >> $configname
    echo 'dx0_n '$dx0_n >> $configname
    echo 'dx0_p '$dx0_p >> $configname
    echo 'dy0 '$dy0 >> $configname
    echo 'dx_sig_n '$dx_sig_n >> $configname
    echo 'dx_sig_p '$dx_sig_p >> $configname
    echo 'dy_sig '$dy_sig >> $configname
    echo 'atime0 '${atime0[$idx]} >> $configname
    echo 'atime_sig '${atime_sig[$idx]} >> $configname
    echo 'useAlshield '${useAlshield[$idx]} >> $configname
    mv $configname 'SBS'$kine'/'

    #LD2 runs
    tar='ld2'
    #SBS field 85%
    field=85
    dx0_n=0
    dx0_p=0
    dx_sig_n=1
    dx_sig_p=1
    dy0=0.0158817
    dy_sig=0.155968
    configname='secal_'$tar'_sbs'$kine'_f'$field'.cfg'
    echo '#sseeds 2.23.23' >> $configname
    echo '#Config file for sbs'$kine' generated from pass'$pass >> $configname
    for a in {11996,11997,11998,11999,12001,12002,12004,12006,12013,12014,12017,12019,12021,12029,12030,12038,12039,12040,12041,12042,12043,12044,12045,12046,12047,12048,12050,12051,12055,12056,12057,12059,12060,12062,12063,12065,12066,12068,12069,12070,12071,12073};
    do
    ls $ld2path*$a* | sort -V >> $configname
    done
    echo 'endlist' >> $configname
    echo ${globalcut[$idx]} >> $configname
    echo 'endcut' >> $configname
    echo 'E_e '${E_e[$idx]} >> $configname
    echo 'HCal_d '${HCal_d[$idx]} >> $configname
    echo 'HCal_th '${HCal_th[$idx]} >> $configname
    echo 'BB_th '${BB_th[$idx]} >> $configname
    echo 'W2_mean '${W2_mean[$idx]} >> $configname
    echo 'W2_sig '${W2_sig[$idx]} >> $configname
    echo 'dx0_n '$dx0_n >> $configname
    echo 'dx0_p '$dx0_p >> $configname
    echo 'dy0 '$dy0 >> $configname
    echo 'dx_sig_n '$dx_sig_n >> $configname
    echo 'dx_sig_p '$dx_sig_p >> $configname
    echo 'dy_sig '$dy_sig >> $configname
    echo 'atime0 '${atime0[$idx]} >> $configname
    echo 'atime_sig '${atime_sig[$idx]} >> $configname
    echo 'useAlshield '${useAlshield[$idx]} >> $configname
    mv $configname 'SBS'$kine'/'

fi

#SBS11
if [ $kine -eq 11 ]
then

    #Index for general cut params
    idx=2
    #LH2 runs,    
    tar='lh2'
    #SBS field 0%
    field=0
    dx0_n=0
    dx0_p=0
    dx_sig_n=1
    dx_sig_p=1
    dy0=-0.0119
    dy_sig=0.08675
    configname='secal_'$tar'_sbs'$kine'_f'$field'.cfg'
    echo '#sseeds 2.23.23' >> $configname
    echo '#Config file for sbs'$kine' generated from pass'$pass >> $configname
    for a in {12730,12731};
    do
    ls $lh2path*$a* | sort -V >> $configname
    done
    echo 'endlist' >> $configname
    echo ${globalcut[$idx]} >> $configname
    echo 'endcut' >> $configname
    echo 'E_e '${E_e[$idx]} >> $configname
    echo 'HCal_d '${HCal_d[$idx]} >> $configname
    echo 'HCal_th '${HCal_th[$idx]} >> $configname
    echo 'BB_th '${BB_th[$idx]} >> $configname
    echo 'W2_mean '${W2_mean[$idx]} >> $configname
    echo 'W2_sig '${W2_sig[$idx]} >> $configname
    echo 'dx0_n '$dx0_n >> $configname
    echo 'dx0_p '$dx0_p >> $configname
    echo 'dy0 '$dy0 >> $configname
    echo 'dx_sig_n '$dx_sig_n >> $configname
    echo 'dx_sig_p '$dx_sig_p >> $configname
    echo 'dy_sig '$dy_sig >> $configname
    echo 'atime0 '${atime0[$idx]} >> $configname
    echo 'atime_sig '${atime_sig[$idx]} >> $configname
    echo 'useAlshield '${useAlshield[$idx]} >> $configname
    mv $configname 'SBS'$kine'/'

    #SBS field 100%
    field=100
    dx0_n=0
    dx0_p=0
    dx_sig_n=1
    dx_sig_p=1
    dy0=-0.0119
    dy_sig=0.08675
    configname='secal_'$tar'_sbs'$kine'_f'$field'.cfg'
    echo '#sseeds 2.23.23' >> $configname
    echo '#Config file for sbs'$kine' generated from pass'$pass >> $configname
    for a in {12313,12320,12335,12336,12340,12345,12355,12358,12363,12367,12368,12369,12370,12380,12382,12400,12401,12414,12415,12427,12428,12429,12471,12472,12491,12496,12497,12524,12525,12526,12548,12561,12563,12575,12576,12587,12588,12625,12626,12657,12658,12675,12676,12677,12691,12692,12704,12705,12716,12747,12758,12765,12776,12777,12809,12810,12812,12822,12823,12895,12910,12923,12930,12931,12959,12960,12973,12974,13027,13028,13041,13042,13043,13056,13057};
    do
    ls $lh2path*$a* | sort -V >> $configname
    done
    echo 'endlist' >> $configname
    echo ${globalcut[$idx]} >> $configname
    echo 'endcut' >> $configname
    echo 'E_e '${E_e[$idx]} >> $configname
    echo 'HCal_d '${HCal_d[$idx]} >> $configname
    echo 'HCal_th '${HCal_th[$idx]} >> $configname
    echo 'BB_th '${BB_th[$idx]} >> $configname
    echo 'W2_mean '${W2_mean[$idx]} >> $configname
    echo 'W2_sig '${W2_sig[$idx]} >> $configname
    echo 'dx0_n '$dx0_n >> $configname
    echo 'dx0_p '$dx0_p >> $configname
    echo 'dy0 '$dy0 >> $configname
    echo 'dx_sig_n '$dx_sig_n >> $configname
    echo 'dx_sig_p '$dx_sig_p >> $configname
    echo 'dy_sig '$dy_sig >> $configname
    echo 'atime0 '${atime0[$idx]} >> $configname
    echo 'atime_sig '${atime_sig[$idx]} >> $configname
    echo 'useAlshield '${useAlshield[$idx]} >> $configname
    mv $configname 'SBS'$kine'/'

    #LD2 runs
    tar='ld2'
    #SBS field 0%
    field=0
    dx0_n=0
    dx0_p=0
    dx_sig_n=1
    dx_sig_p=1
    dy0=-0.0119
    dy_sig=0.08675
    configname='secal_'$tar'_sbs'$kine'_f'$field'.cfg'
    echo '#sseeds 2.23.23' >> $configname
    echo '#Config file for sbs'$kine' generated from pass'$pass >> $configname
    for a in {12726,12728,12729};
    do
    ls $ld2path*$a* | sort -V >> $configname
    done
    echo 'endlist' >> $configname
    echo ${globalcut[$idx]} >> $configname
    echo 'endcut' >> $configname
    echo 'E_e '${E_e[$idx]} >> $configname
    echo 'HCal_d '${HCal_d[$idx]} >> $configname
    echo 'HCal_th '${HCal_th[$idx]} >> $configname
    echo 'BB_th '${BB_th[$idx]} >> $configname
    echo 'W2_mean '${W2_mean[$idx]} >> $configname
    echo 'W2_sig '${W2_sig[$idx]} >> $configname
    echo 'dx0_n '$dx0_n >> $configname
    echo 'dx0_p '$dx0_p >> $configname
    echo 'dy0 '$dy0 >> $configname
    echo 'dx_sig_n '$dx_sig_n >> $configname
    echo 'dx_sig_p '$dx_sig_p >> $configname
    echo 'dy_sig '$dy_sig >> $configname
    echo 'atime0 '${atime0[$idx]} >> $configname
    echo 'atime_sig '${atime_sig[$idx]} >> $configname
    echo 'useAlshield '${useAlshield[$idx]} >> $configname
    mv $configname 'SBS'$kine'/'

    #SBS field 100%
    field=100
    dx0_n=0
    dx0_p=0
    dx_sig_n=1
    dx_sig_p=1
    dy0=-0.0119
    dy_sig=0.08675
    configname='secal_'$tar'_sbs'$kine'_f'$field'.cfg'
    echo '#sseeds 2.23.23' >> $configname
    echo '#Config file for sbs'$kine' generated from pass'$pass >> $configname
    for a in {12314,12315,12316,12318,12319,12321,12322,12323,12333,12337,12338,12339,12341,12342,12343,12344,12346,12347,12348,12349,12356,12357,12359,12360,12362,12364,12366,12371,12372,12373,12374,12375,12376,12377,12378,12379,12383,12384,12385,12387,12388,12394,12395,12396,12403,12404,12405,12406,12407,12408,12409,12410,12411,12412,12416,12423,12424,12425,12473,12474,12478,12479,12480,12481,12482,12485,12486,12498,12499,12500,12501,12502,12512,12513,12514,12515,12516,12517,12521,12523,12527,12532,12534,12535,12536,12537,12538,12539,12541,12542,12543,12544,12547,12549,12550,12551,12552,12553,12554,12556,12558,12559,12565,12566,12567,12568,12569,12570,12571,12572,12574,12577,12580,12581,12582,12583,12584,12585,12586,12589,12590,12592,12593,12597,12598,12599,12614,12618,12619,12620,12623,12627,12628,12638,12639,12641,12642,12643,12644,12645,12647,12653,12654,12655,12656,12659,12660,12661,12662,12664,12666,12671,12673,12674,12678,12679,12680,12681,12682,12683,12684,12686,12687,12688,12693,12694,12695,12697,12699,12700,12701,12702,12706,12707,12708,12709,12710,12711,12713,12714,12715,12721,12722,12732,12733,12734,12735,12737,12738,12739,12740,12742,12743,12767,12768,12769,12770,12771,12772,12773,12778,12779,12781,12782,12796,12797,12798,12799,12801,12802,12813,12814,12815,12816,12817,12819,12820,12824,12825,12826,12827,12828,12829,12830,12894,12896,12897,12902,12903,12904,12906,12907,12908,12909,12911,12912,12913,12914,12915,12916,12917,12918,12919,12920,12921,12932,12933,12934,12935,12936,12939,12940,12946,12957,12961,12962,12963,12965,12968,12969,12970,12971,12972,12975,12976,12977,12979,12980,12982,12983,12984,12986,12987,12988,12989,13000,13029,13030,13034,13035,13036,13038,13039,13040,13044,13045,13048,13049,13050,13051,13052,13053,13054,13055,13058,13059,13060,13061,13062,13063};
    do
    ls $ld2path*$a* | sort -V >> $configname
    done
    echo 'endlist' >> $configname
    echo ${globalcut[$idx]} >> $configname
    echo 'endcut' >> $configname
    echo 'E_e '${E_e[$idx]} >> $configname
    echo 'HCal_d '${HCal_d[$idx]} >> $configname
    echo 'HCal_th '${HCal_th[$idx]} >> $configname
    echo 'BB_th '${BB_th[$idx]} >> $configname
    echo 'W2_mean '${W2_mean[$idx]} >> $configname
    echo 'W2_sig '${W2_sig[$idx]} >> $configname
    echo 'dx0_n '$dx0_n >> $configname
    echo 'dx0_p '$dx0_p >> $configname
    echo 'dy0 '$dy0 >> $configname
    echo 'dx_sig_n '$dx_sig_n >> $configname
    echo 'dx_sig_p '$dx_sig_p >> $configname
    echo 'dy_sig '$dy_sig >> $configname
    echo 'atime0 '${atime0[$idx]} >> $configname
    echo 'atime_sig '${atime_sig[$idx]} >> $configname
    echo 'useAlshield '${useAlshield[$idx]} >> $configname
    mv $configname 'SBS'$kine'/'

fi

#SBS14
if [ $kine -eq 14 ]
then

    #Index for general cut params
    idx=3
    #LH2 runs
    tar='lh2'
    #SBS field 0%
    field=0
    dx0_n=0
    dx0_p=0
    dx_sig_n=1
    dx_sig_p=1
    dy0=0.00220912
    dy_sig=0.0824202
    configname='secal_'$tar'_sbs'$kine'_f'$field'.cfg'
    echo '#sseeds 2.23.23' >> $configname
    echo '#Config file for sbs'$kine' generated from pass'$pass >> $configname
    for a in {13375,13376,13377,13378};
    do
    ls $lh2path*$a* | sort -V >> $configname
    done
    echo 'endlist' >> $configname
    echo ${globalcut[$idx]} >> $configname
    echo 'endcut' >> $configname
    echo 'E_e '${E_e[$idx]} >> $configname
    echo 'HCal_d '${HCal_d[$idx]} >> $configname
    echo 'HCal_th '${HCal_th[$idx]} >> $configname
    echo 'BB_th '${BB_th[$idx]} >> $configname
    echo 'W2_mean '${W2_mean[$idx]} >> $configname
    echo 'W2_sig '${W2_sig[$idx]} >> $configname
    echo 'dx0_n '$dx0_n >> $configname
    echo 'dx0_p '$dx0_p >> $configname
    echo 'dy0 '$dy0 >> $configname
    echo 'dx_sig_n '$dx_sig_n >> $configname
    echo 'dx_sig_p '$dx_sig_p >> $configname
    echo 'dy_sig '$dy_sig >> $configname
    echo 'atime0 '${atime0[$idx]} >> $configname
    echo 'atime_sig '${atime_sig[$idx]} >> $configname
    echo 'useAlshield '${useAlshield[$idx]} >> $configname
    mv $configname 'SBS'$kine'/'

    #SBS field 70%
    field=70
    dx0_n=0
    dx0_p=0
    dx_sig_n=1
    dx_sig_p=1
    dy0=0.00220912
    dy_sig=0.0824202
    configname='secal_'$tar'_sbs'$kine'_f'$field'.cfg'
    echo '#sseeds 2.23.23' >> $configname
    echo '#Config file for sbs'$kine' generated from pass'$pass >> $configname
    for a in {13239,13240,13241,13242,13243,13244,13312,13313,13320,13321,13345,13346,13348,13349,13351,13352,13379,13396,13397,13405};
    do
    ls $lh2path*$a* | sort -V >> $configname
    done
    echo 'endlist' >> $configname
    echo ${globalcut[$idx]} >> $configname
    echo 'endcut' >> $configname
    echo 'E_e '${E_e[$idx]} >> $configname
    echo 'HCal_d '${HCal_d[$idx]} >> $configname
    echo 'HCal_th '${HCal_th[$idx]} >> $configname
    echo 'BB_th '${BB_th[$idx]} >> $configname
    echo 'W2_mean '${W2_mean[$idx]} >> $configname
    echo 'W2_sig '${W2_sig[$idx]} >> $configname
    echo 'dx0_n '$dx0_n >> $configname
    echo 'dx0_p '$dx0_p >> $configname
    echo 'dy0 '$dy0 >> $configname
    echo 'dx_sig_n '$dx_sig_n >> $configname
    echo 'dx_sig_p '$dx_sig_p >> $configname
    echo 'dy_sig '$dy_sig >> $configname
    echo 'atime0 '${atime0[$idx]} >> $configname
    echo 'atime_sig '${atime_sig[$idx]} >> $configname
    echo 'useAlshield '${useAlshield[$idx]} >> $configname
    mv $configname 'SBS'$kine'/'

    #LD2 runs
    tar='ld2'
    #SBS field 70%
    field=70
    dx0_n=0
    dx0_p=0
    dx_sig_n=1
    dx_sig_p=1
    dy0=0.00220912
    dy_sig=0.0824202
    configname='secal_'$tar'_sbs'$kine'_f'$field'.cfg'
    echo '#sseeds 2.23.23' >> $configname
    echo '#Config file for sbs'$kine' generated from pass'$pass >> $configname
    for a in {13305,13306,13307,13308,13309,13314,13315,13316,13317,13318,13319,13322,13323,13324,13325,13342,13344,13353,13357,13358,13359,13360,13361,13362,13363,13364,13368,13369,13370,13371,13372,13373,13381,13382,13384,13385,13387,13388,13389,13390,13391,13392,13393,13394,13395,13398,13399,13400,13402,13403,13406,13407};
    do
    ls $ld2path*$a* | sort -V >> $configname
    done
    echo 'endlist' >> $configname
    echo ${globalcut[$idx]} >> $configname
    echo 'endcut' >> $configname
    echo 'E_e '${E_e[$idx]} >> $configname
    echo 'HCal_d '${HCal_d[$idx]} >> $configname
    echo 'HCal_th '${HCal_th[$idx]} >> $configname
    echo 'BB_th '${BB_th[$idx]} >> $configname
    echo 'W2_mean '${W2_mean[$idx]} >> $configname
    echo 'W2_sig '${W2_sig[$idx]} >> $configname
    echo 'dx0_n '$dx0_n >> $configname
    echo 'dx0_p '$dx0_p >> $configname
    echo 'dy0 '$dy0 >> $configname
    echo 'dx_sig_n '$dx_sig_n >> $configname
    echo 'dx_sig_p '$dx_sig_p >> $configname
    echo 'dy_sig '$dy_sig >> $configname
    echo 'atime0 '${atime0[$idx]} >> $configname
    echo 'atime_sig '${atime_sig[$idx]} >> $configname
    echo 'useAlshield '${useAlshield[$idx]} >> $configname
    mv $configname 'SBS'$kine'/'

fi

#SBS8
if [ $kine -eq 8 ]
then

    #Index for general cut params
    idx=4
    #LH2 runs
    tar='lh2'
    #SBS field 0%
    field=0
    dx0_n=0
    dx0_p=0
    dx_sig_n=1
    dx_sig_p=1
    dy0=-0.0119
    dy_sig=0.08675
    configname='secal_'$tar'_sbs'$kine'_f'$field'.cfg'
    echo '#sseeds 2.23.23' >> $configname
    echo '#Config file for sbs'$kine' generated from pass'$pass >> $configname
    for a in {13459,13460,13461,13463,13464,13465,13466};
    do
    ls $lh2path*$a* | sort -V >> $configname
    done
    echo 'endlist' >> $configname
    echo ${globalcut[$idx]} >> $configname
    echo 'endcut' >> $configname
    echo 'E_e '${E_e[$idx]} >> $configname
    echo 'HCal_d '${HCal_d[$idx]} >> $configname
    echo 'HCal_th '${HCal_th[$idx]} >> $configname
    echo 'BB_th '${BB_th[$idx]} >> $configname
    echo 'W2_mean '${W2_mean[$idx]} >> $configname
    echo 'W2_sig '${W2_sig[$idx]} >> $configname
    echo 'dx0_n '$dx0_n >> $configname
    echo 'dx0_p '$dx0_p >> $configname
    echo 'dy0 '$dy0 >> $configname
    echo 'dx_sig_n '$dx_sig_n >> $configname
    echo 'dx_sig_p '$dx_sig_p >> $configname
    echo 'dy_sig '$dy_sig >> $configname
    echo 'atime0 '${atime0[$idx]} >> $configname
    echo 'atime_sig '${atime_sig[$idx]} >> $configname
    echo 'useAlshield '${useAlshield[$idx]} >> $configname
    mv $configname 'SBS'$kine'/'

    #SBS field 50%
    field=50
    dx0_n=0
    dx0_p=0
    dx_sig_n=1
    dx_sig_p=1
    dy0=-0.0119
    dy_sig=0.08675
    configname='secal_'$tar'_sbs'$kine'_f'$field'.cfg'
    echo '#sseeds 2.23.23' >> $configname
    echo '#Config file for sbs'$kine' generated from pass'$pass >> $configname
    #for a in {13580};
    #do
    #ls $lh2path*$a* | sort -V >> $configname
    ls $lh2path*13580* | sort -V >> $configname
    #done
    echo 'endlist' >> $configname
    echo ${globalcut[$idx]} >> $configname
    echo 'endcut' >> $configname
    echo 'E_e '${E_e[$idx]} >> $configname
    echo 'HCal_d '${HCal_d[$idx]} >> $configname
    echo 'HCal_th '${HCal_th[$idx]} >> $configname
    echo 'BB_th '${BB_th[$idx]} >> $configname
    echo 'W2_mean '${W2_mean[$idx]} >> $configname
    echo 'W2_sig '${W2_sig[$idx]} >> $configname
    echo 'dx0_n '$dx0_n >> $configname
    echo 'dx0_p '$dx0_p >> $configname
    echo 'dy0 '$dy0 >> $configname
    echo 'dx_sig_n '$dx_sig_n >> $configname
    echo 'dx_sig_p '$dx_sig_p >> $configname
    echo 'dy_sig '$dy_sig >> $configname
    echo 'atime0 '${atime0[$idx]} >> $configname
    echo 'atime_sig '${atime_sig[$idx]} >> $configname
    echo 'useAlshield '${useAlshield[$idx]} >> $configname
    mv $configname 'SBS'$kine'/'

    #SBS field 70%
    field=70
    dx0_n=0
    dx0_p=0
    dx_sig_n=1
    dx_sig_p=1
    dy0=-0.0119
    dy_sig=0.08675
    configname='secal_'$tar'_sbs'$kine'_f'$field'.cfg'
    echo '#sseeds 2.23.23' >> $configname
    echo '#Config file for sbs'$kine' generated from pass'$pass >> $configname
    for a in {13450,13451,13452,13482,13483,13484,13485,13486,13487,13488,13489,13490,13573,13574,13575,13576,13577,13578};
    do
    ls $lh2path*$a* | sort -V >> $configname
    done
    echo 'endlist' >> $configname
    echo ${globalcut[$idx]} >> $configname
    echo 'endcut' >> $configname
    echo 'E_e '${E_e[$idx]} >> $configname
    echo 'HCal_d '${HCal_d[$idx]} >> $configname
    echo 'HCal_th '${HCal_th[$idx]} >> $configname
    echo 'BB_th '${BB_th[$idx]} >> $configname
    echo 'W2_mean '${W2_mean[$idx]} >> $configname
    echo 'W2_sig '${W2_sig[$idx]} >> $configname
    echo 'dx0_n '$dx0_n >> $configname
    echo 'dx0_p '$dx0_p >> $configname
    echo 'dy0 '$dy0 >> $configname
    echo 'dx_sig_n '$dx_sig_n >> $configname
    echo 'dx_sig_p '$dx_sig_p >> $configname
    echo 'dy_sig '$dy_sig >> $configname
    echo 'atime0 '${atime0[$idx]} >> $configname
    echo 'atime_sig '${atime_sig[$idx]} >> $configname
    echo 'useAlshield '${useAlshield[$idx]} >> $configname
    mv $configname 'SBS'$kine'/'

    #SBS field 100%
    field=100
    dx0_n=0
    dx0_p=0
    dx_sig_n=1
    dx_sig_p=1
    dy0=-0.0119
    dy_sig=0.08675
    configname='secal_'$tar'_sbs'$kine'_f'$field'.cfg'
    echo '#sseeds 2.23.23' >> $configname
    echo '#Config file for sbs'$kine' generated from pass'$pass >> $configname
    for a in {13537,13538,13539,13540,13541,13542,13543};
    do
    ls $lh2path*$a* | sort -V >> $configname
    done
    echo 'endlist' >> $configname
    echo ${globalcut[$idx]} >> $configname
    echo 'endcut' >> $configname
    echo 'E_e '${E_e[$idx]} >> $configname
    echo 'HCal_d '${HCal_d[$idx]} >> $configname
    echo 'HCal_th '${HCal_th[$idx]} >> $configname
    echo 'BB_th '${BB_th[$idx]} >> $configname
    echo 'W2_mean '${W2_mean[$idx]} >> $configname
    echo 'W2_sig '${W2_sig[$idx]} >> $configname
    echo 'dx0_n '$dx0_n >> $configname
    echo 'dx0_p '$dx0_p >> $configname
    echo 'dy0 '$dy0 >> $configname
    echo 'dx_sig_n '$dx_sig_n >> $configname
    echo 'dx_sig_p '$dx_sig_p >> $configname
    echo 'dy_sig '$dy_sig >> $configname
    echo 'atime0 '${atime0[$idx]} >> $configname
    echo 'atime_sig '${atime_sig[$idx]} >> $configname
    echo 'useAlshield '${useAlshield[$idx]} >> $configname
    mv $configname 'SBS'$kine'/'

    #LD2 runs
    tar='ld2'
    #SBS field 0%
    field=0
    dx0_n=0
    dx0_p=0
    dx_sig_n=1
    dx_sig_p=1
    dy0=-0.0119
    dy_sig=0.08675
    configname='secal_'$tar'_sbs'$kine'_f'$field'.cfg'
    echo '#sseeds 2.23.23' >> $configname
    echo '#Config file for sbs'$kine' generated from pass'$pass >> $configname
    for a in {13468,13470,13472,13473,13474,13475,13476,13477,13478,13479};
    do
    ls $ld2path*$a* | sort -V >> $configname
    done
    echo 'endlist' >> $configname
    echo ${globalcut[$idx]} >> $configname
    echo 'endcut' >> $configname
    echo 'E_e '${E_e[$idx]} >> $configname
    echo 'HCal_d '${HCal_d[$idx]} >> $configname
    echo 'HCal_th '${HCal_th[$idx]} >> $configname
    echo 'BB_th '${BB_th[$idx]} >> $configname
    echo 'W2_mean '${W2_mean[$idx]} >> $configname
    echo 'W2_sig '${W2_sig[$idx]} >> $configname
    echo 'dx0_n '$dx0_n >> $configname
    echo 'dx0_p '$dx0_p >> $configname
    echo 'dy0 '$dy0 >> $configname
    echo 'dx_sig_n '$dx_sig_n >> $configname
    echo 'dx_sig_p '$dx_sig_p >> $configname
    echo 'dy_sig '$dy_sig >> $configname
    echo 'atime0 '${atime0[$idx]} >> $configname
    echo 'atime_sig '${atime_sig[$idx]} >> $configname
    echo 'useAlshield '${useAlshield[$idx]} >> $configname
    mv $configname 'SBS'$kine'/'

    #SBS field 50%
    field=50
    dx0_n=0
    dx0_p=0
    dx_sig_n=1
    dx_sig_p=1
    dy0=-0.0119
    dy_sig=0.08675
    configname='secal_'$tar'_sbs'$kine'_f'$field'.cfg'
    echo '#sseeds 2.23.23' >> $configname
    echo '#Config file for sbs'$kine' generated from pass'$pass >> $configname
    for a in {13581,13582,13583,13584,13585,13586};
    do
    ls $ld2path*$a* | sort -V >> $configname
    done
    echo 'endlist' >> $configname
    echo ${globalcut[$idx]} >> $configname
    echo 'endcut' >> $configname
    echo 'E_e '${E_e[$idx]} >> $configname
    echo 'HCal_d '${HCal_d[$idx]} >> $configname
    echo 'HCal_th '${HCal_th[$idx]} >> $configname
    echo 'BB_th '${BB_th[$idx]} >> $configname
    echo 'W2_mean '${W2_mean[$idx]} >> $configname
    echo 'W2_sig '${W2_sig[$idx]} >> $configname
    echo 'dx0_n '$dx0_n >> $configname
    echo 'dx0_p '$dx0_p >> $configname
    echo 'dy0 '$dy0 >> $configname
    echo 'dx_sig_n '$dx_sig_n >> $configname
    echo 'dx_sig_p '$dx_sig_p >> $configname
    echo 'dy_sig '$dy_sig >> $configname
    echo 'atime0 '${atime0[$idx]} >> $configname
    echo 'atime_sig '${atime_sig[$idx]} >> $configname
    echo 'useAlshield '${useAlshield[$idx]} >> $configname
    mv $configname 'SBS'$kine'/'

    #SBS field 70%
    field=70
    dx0_n=0
    dx0_p=0
    dx_sig_n=1
    dx_sig_p=1
    dy0=-0.0119
    dy_sig=0.08675
    configname='secal_'$tar'_sbs'$kine'_f'$field'.cfg'
    echo '#sseeds 2.23.23' >> $configname
    echo '#Config file for sbs'$kine' generated from pass'$pass >> $configname
    for a in {13453,13454,13455,13491,13492,13493,13494,13495,13496,13497,13502,13503,13504,13505,13558,13559,13560,13561,13562,13563,13564,13565,13566,13567,13568,13569,13570,13571,13587,13588,13589,13590,13591,13592,13593,13596,13597,13608,13609,13610,13612,13613,13614,13615,13616,13617,13618,13619,13620};
    do
    ls $ld2path*$a* | sort -V >> $configname
    done
    echo 'endlist' >> $configname
    echo ${globalcut[$idx]} >> $configname
    echo 'endcut' >> $configname
    echo 'E_e '${E_e[$idx]} >> $configname
    echo 'HCal_d '${HCal_d[$idx]} >> $configname
    echo 'HCal_th '${HCal_th[$idx]} >> $configname
    echo 'BB_th '${BB_th[$idx]} >> $configname
    echo 'W2_mean '${W2_mean[$idx]} >> $configname
    echo 'W2_sig '${W2_sig[$idx]} >> $configname
    echo 'dx0_n '$dx0_n >> $configname
    echo 'dx0_p '$dx0_p >> $configname
    echo 'dy0 '$dy0 >> $configname
    echo 'dx_sig_n '$dx_sig_n >> $configname
    echo 'dx_sig_p '$dx_sig_p >> $configname
    echo 'dy_sig '$dy_sig >> $configname
    echo 'atime0 '${atime0[$idx]} >> $configname
    echo 'atime_sig '${atime_sig[$idx]} >> $configname
    echo 'useAlshield '${useAlshield[$idx]} >> $configname
    mv $configname 'SBS'$kine'/'

    #SBS field 100%
    field=100
    dx0_n=0
    dx0_p=0
    dx_sig_n=1
    dx_sig_p=1
    dy0=-0.0119
    dy_sig=0.08675
    configname='secal_'$tar'_sbs'$kine'_f'$field'.cfg'
    echo '#sseeds 2.23.23' >> $configname
    echo '#Config file for sbs'$kine' generated from pass'$pass >> $configname
    for a in {13544,13545,13546,13547,13548,13549,13550,13551,13552,13554,13556,13557};
    do
    ls $ld2path*$a* | sort -V >> $configname
    done
    echo 'endlist' >> $configname
    echo ${globalcut[$idx]} >> $configname
    echo 'endcut' >> $configname
    echo 'E_e '${E_e[$idx]} >> $configname
    echo 'HCal_d '${HCal_d[$idx]} >> $configname
    echo 'HCal_th '${HCal_th[$idx]} >> $configname
    echo 'BB_th '${BB_th[$idx]} >> $configname
    echo 'W2_mean '${W2_mean[$idx]} >> $configname
    echo 'W2_sig '${W2_sig[$idx]} >> $configname
    echo 'dx0_n '$dx0_n >> $configname
    echo 'dx0_p '$dx0_p >> $configname
    echo 'dy0 '$dy0 >> $configname
    echo 'dx_sig_n '$dx_sig_n >> $configname
    echo 'dx_sig_p '$dx_sig_p >> $configname
    echo 'dy_sig '$dy_sig >> $configname
    echo 'atime0 '${atime0[$idx]} >> $configname
    echo 'atime_sig '${atime_sig[$idx]} >> $configname
    echo 'useAlshield '${useAlshield[$idx]} >> $configname
    mv $configname 'SBS'$kine'/'

fi

#SBS9
if [ $kine -eq 9 ]
then

    #Index for general cut params
    idx=5
    #LH2 runs
    tar='lh2'
    #SBS field 70%
    field=70
    dx0_n=0
    dx0_p=0
    dx_sig_n=1
    dx_sig_p=1
    dy0=-0.0175041
    dy_sig=0.134702
    configname='secal_'$tar'_sbs'$kine'_f'$field'.cfg'
    echo '#sseeds 2.23.23' >> $configname
    echo '#Config file for sbs'$kine' generated from pass'$pass >> $configname
    for a in {13656,13657,13663,13676,13683,13696,13697,13719,13720,13730,13747,13768,13769,13795,13796};
    do
    ls $lh2path*$a* | sort -V >> $configname
    done
    echo 'endlist' >> $configname
    echo ${globalcut[$idx]} >> $configname
    echo 'endcut' >> $configname
    echo 'E_e '${E_e[$idx]} >> $configname
    echo 'HCal_d '${HCal_d[$idx]} >> $configname
    echo 'HCal_th '${HCal_th[$idx]} >> $configname
    echo 'BB_th '${BB_th[$idx]} >> $configname
    echo 'W2_mean '${W2_mean[$idx]} >> $configname
    echo 'W2_sig '${W2_sig[$idx]} >> $configname
    echo 'dx0_n '$dx0_n >> $configname
    echo 'dx0_p '$dx0_p >> $configname
    echo 'dy0 '$dy0 >> $configname
    echo 'dx_sig_n '$dx_sig_n >> $configname
    echo 'dx_sig_p '$dx_sig_p >> $configname
    echo 'dy_sig '$dy_sig >> $configname
    echo 'atime0 '${atime0[$idx]} >> $configname
    echo 'atime_sig '${atime_sig[$idx]} >> $configname
    echo 'useAlshield '${useAlshield[$idx]} >> $configname
    mv $configname 'SBS'$kine'/'

    #LD2 runs
    tar='ld2'
    #SBS field 70%
    field=70
    dx0_n=0
    dx0_p=0
    dx_sig_n=1
    dx_sig_p=1
    dy0=-0.0175041
    dy_sig=0.134702
    configname='secal_'$tar'_sbs'$kine'_f'$field'.cfg'
    echo '#sseeds 2.23.23' >> $configname
    echo '#Config file for sbs'$kine' generated from pass'$pass >> $configname
    for a in {13660,13661,13662,13664,13665,13666,13677,13678,13679,13680,13681,13682,13684,13685,13686,13687,13688,13689,13694,13695,13698,13699,13700,13710,13711,13712,13714,13715,13716,13717,13721,13723,13724,13727,13728,13729,13731,13732,13734,13736,13737,13746,13748,13749,13753,13754,13755,13756,13757,13758,13760,13761,13764,13765,13766,13767,13770,13771,13773,13775,13776,13777,13778,13779,13793,13797,13798,13799};
    do
    ls $ld2path*$a* | sort -V >> $configname
    done
    echo 'endlist' >> $configname
    echo ${globalcut[$idx]} >> $configname
    echo 'endcut' >> $configname
    echo 'E_e '${E_e[$idx]} >> $configname
    echo 'HCal_d '${HCal_d[$idx]} >> $configname
    echo 'HCal_th '${HCal_th[$idx]} >> $configname
    echo 'BB_th '${BB_th[$idx]} >> $configname
    echo 'W2_mean '${W2_mean[$idx]} >> $configname
    echo 'W2_sig '${W2_sig[$idx]} >> $configname
    echo 'dx0_n '$dx0_n >> $configname
    echo 'dx0_p '$dx0_p >> $configname
    echo 'dy0 '$dy0 >> $configname
    echo 'dx_sig_n '$dx_sig_n >> $configname
    echo 'dx_sig_p '$dx_sig_p >> $configname
    echo 'dy_sig '$dy_sig >> $configname
    echo 'atime0 '${atime0[$idx]} >> $configname
    echo 'atime_sig '${atime_sig[$idx]} >> $configname
    echo 'useAlshield '${useAlshield[$idx]} >> $configname
    mv $configname 'SBS'$kine'/'

fi
