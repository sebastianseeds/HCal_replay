#!/bin/sh

# SSeeds - 5.28.23 - simple shell script to run tcal.C

## Usage
#./run_tdc_align.sh <string: experiment> <int: configuration> <bool: quasi-replay> <int: replay pass> <int: number of energy calibrations> <bool: parameter override option> <string: parameter override timestamp>

exp=$1
config=$2
qreplay=$3
pass=$4
nset=$5
pover=$6
pts=$7

root -l 'tdc_tw.C("'$exp'",'$config','$qreplay','$pass','$nset',"'$pover'","'$pts'")'
