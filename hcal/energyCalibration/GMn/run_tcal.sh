#!/bin/sh

# SSeeds - 2.27.23 - simple shell script to run tcal.C for use with batch farm scripts

## Usage
#./run_tcal.sh <kine> <iter>

kine=$1
iter=$2

cd /w/halla-scshelf2102/sbs/seeds/HCal_replay/hcal/energyCalibration/GMn

root -l -b -q 'tcal.C('$kine','$iter')'

