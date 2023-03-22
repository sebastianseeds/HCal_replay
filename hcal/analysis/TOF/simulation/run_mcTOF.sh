#!/bin/sh

# SSeeds - 3.9.23 - simple shell script to run simTOF.C for use with batch farm scripts

## Usage
#./run_mcTOF.sh <kine> <field> <iter>

kine=$1
field=$2
iter=$3

cd /w/halla-scshelf2102/sbs/seeds/HCal_replay/hcal/analysis/TOF/simulation

root -l -b -q 'simTOF.C('$kine','$field','$iter')'

