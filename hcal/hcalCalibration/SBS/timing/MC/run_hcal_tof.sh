#!/bin/sh

# SSeeds - 5.23.23 - simple shell script to run hcal_tof.C

## Usage
#./run_hcal_tof.sh <kine> <field> <iter>

kine=$1
field=$2
iter=$3

root -l 'hcal_tof.C('$kine','$field','$iter')'

