#!/bin/sh

# SSeeds - 5.15.23: Script to display .csv in human readable format

## Usage
#./run_mcreplayed_hde.sh <path-to-csv>

path=$1

column -s, -t < $1 | less -#2 -N -S
