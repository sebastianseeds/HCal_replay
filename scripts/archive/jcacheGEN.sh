#!/bin/sh

# SSeeds - 7.22.22

## Usage
# Simple script takes GEN data from mss and places it into the cache 
#./jcache.sh <run> <nsegments>

echo -e "\n Requesting cache up to segment $2 from zero of run $1... \n"

for (( i=0; i<=$2; i++ ))
do
    #echo "$i"
    jcache get /mss/halla/sbs/raw/e1209016_$1.evio.0.$i
    sleep 5s
done
