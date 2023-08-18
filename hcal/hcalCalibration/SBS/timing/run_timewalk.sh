#!/bin/sh

# SSeeds - 5.28.23 - simple shell script to run tcal.C

## Usage
#./run_timewalk.sh <string: experiment> <int: configuration> <bool: quasi-replay> <int: replay pass> <int: number of energy calibrations>
# Note: remaining parameters should not be changed on this version of the program unless changes to code to enable their use is performed

exp=$1
config=$2
qreplay=$3
pass=$4

root -l --web=off 'timewalk.C("'$exp'",'$config','$qreplay','$pass')'
