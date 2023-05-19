#!/bin/sh

# SSeeds - 5.28.23 - simple shell script to run tcal.C

## Usage
#./run_tdc_align.sh <string: experiment> <int: configuration> <bool: quasi-replay> <int: replay-pass>

exp=$1
config=$2
qreplay=$3
pass=$4

root -l 'tdc_align.C("'$exp'",'$config','$qreplay','$pass')'
