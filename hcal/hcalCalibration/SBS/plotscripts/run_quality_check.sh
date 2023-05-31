#!/bin/sh

# SSeeds - 5.28.23 - simple shell script to run quality_check.C

## Usage
#./run_quality_check.sh <string: experiment> <int: configuration> <int: replay-pass>

exp=$1
config=$2
pass=$3

root -l 'quality_check.C("'$exp'",'$config','$pass')'
