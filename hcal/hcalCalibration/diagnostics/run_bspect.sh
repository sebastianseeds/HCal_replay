#!/bin/sh

# SSeeds - 3.4.23 - simple shell script to run ecal_blockE_spectra.C for use with batch farm scripts

## Usage
#./run_bspect.sh <kine>

kine=$1

cd /w/halla-scshelf2102/sbs/seeds/HCal_replay/hcal/hcalCalibration/diagnostics

root -l -b -q 'ecal_blockE_spectra.C('$kine')'

