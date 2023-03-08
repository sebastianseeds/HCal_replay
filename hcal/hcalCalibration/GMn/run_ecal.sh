#!/bin/sh

# SSeeds - 2.24.23 - simple shell script to run ecal.C for use with batch farm scripts

## Usage
#./run_ecal.sh <kine> <iter>

kine=$1
iter=$2

cd /w/halla-scshelf2102/sbs/seeds/HCal_replay/hcal/hcalCalibration/GMn

root -l -b -q 'ecal.C('$kine','$iter')'

