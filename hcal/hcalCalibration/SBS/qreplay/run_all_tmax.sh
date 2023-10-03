#!/bin/sh

# SSeeds - 9.23.23 - run all timing calibration sets tmax.C

echo "Setting up tmax run all with a 5 run maximum"

echo "Running tmax.C on SBS4 set 1"

root -l -b -q 'tmax.C("gmn",4,0,1,5,0,0)'

wait

echo "Running tmax.C on SBS7 set 1"

root -l -b -q 'tmax.C("gmn",7,0,1,5,0,0)'

wait

echo "Running tmax.C on SBS14 set 1"

root -l -b -q 'tmax.C("gmn",14,1,1,5,0,0)'

wait

echo "Running tmax.C on SBS14 set 2"

root -l -b -q 'tmax.C("gmn",14,1,2,5,0,0)'

wait

echo "Running tmax.C on SBS11 set 1"

root -l -b -q 'tmax.C("gmn",11,1,1,5,0,0)'

wait

echo "Running tmax.C on SBS11 set 2"

root -l -b -q 'tmax.C("gmn",11,1,2,5,0,0)'

wait

echo "Running tmax.C on SBS8 set 1"

root -l -b -q 'tmax.C("gmn",8,1,1,5,0,0)'

wait

echo "Running tmax.C on SBS9 set 1"

root -l -b -q 'tmax.C("gmn",9,1,1,5,0,0)'

echo "Finished!"
