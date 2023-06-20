#!/bin/sh

# SSeeds - 5.28.23 - simple shell script to run ecal.C

## Usage
#./run_ecal.sh <string: experiment> <int: configuration> <bool: quasi-replay> <int: replay-pass> <bool: hydrogen-only calibration>

exp=$1
config=$2
qreplay=$3
pass=$4
h2only=$5

root -l 'ecal.C("'$exp'",'$config','$qreplay','$pass','$h2only')'
