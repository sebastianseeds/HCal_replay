#!/bin/sh

# SSeeds - 5.28.23 - simple shell script to run hcalE_mc.C

## Usage
#./run_hcalE_mc.sh <string: experiment> <int: configuration>

exp=$1
config=$2

root -l 'hcalE_digmc.C("'$exp'",'$config')'
