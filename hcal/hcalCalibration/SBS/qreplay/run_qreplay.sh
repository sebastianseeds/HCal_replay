#!/bin/sh

# SSeeds - 8.5.23 - simple shell script to run ecal.C

## Usage
#./run_qreplay.sh <experiment> <configuration> <replay-pass> <target-option> <selected-energy-timestamp> <replacement-experiment> <replacement-configuration> <replacement-pass> <replacement-database-timestamp> <verbose-option>

exper=$1
config=$2
pass=$3
sts=$4
rexper=$5
rconfig=$6
rpass=$7
qts=$8
h2only=$9
verbose=true

cd /w/halla-scshelf2102/sbs/seeds/HCal_replay/hcal/hcalCalibration/SBS/qreplay

root -l -b -q 'qreplay_standalone.C("'$exper'",'$config','$pass',"'$sts'","'$rexper'",'$rconfig','$rpass',"'$qts'",'$h2only')'
