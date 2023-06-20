#!/bin/sh

# SSeeds - 5.28.23 - simple shell script to run quality_check.C (and quality_canv.C, for now) intended to save plots only

## Usage
#./run_quality_check.sh <string: experiment> <int: configuration> <int: replay-pass> <bool: hydrogen only calibration>

exp=$1
config=$2
pass=$3
h2only=$4

root -l -b -q 'quality_check.C("'$exp'",'$config','$pass','$h2only')'

#Can safely comment lines below unless adct/tdc fit plots do not exist before quality check

wait

root -l -b -q 'quality_canv.C("'$exp'",'$config','$pass','$h2only')'
