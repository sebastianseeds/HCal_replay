#!/bin/sh

# SSeeds - 5.28.23 - simple shell script to run adct_align.C

## Usage
#./run_adct_align.sh <string: experiment> <int: configuration> <bool: quasi-replay> <int: replay-pass>

exp=$1
config=$2
qreplay=$3
pass=$4

root -l 'adct_align.C("'$exp'",'$config','$qreplay','$pass')'
