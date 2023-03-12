#!/bin/sh

# SSeeds - 3.9.23 - simple shell script to run simTOF.C for use with batch farm scripts

## Usage
#./run_mcTOF.sh <kine> <iter>

kine=$1
iter=$2

cd /w/halla-scshelf2102/sbs/seeds/HCal_replay/hcal/analysis/TOF/simulation

root -l -b -q 'simTOF.C('$kine','$iter')'

