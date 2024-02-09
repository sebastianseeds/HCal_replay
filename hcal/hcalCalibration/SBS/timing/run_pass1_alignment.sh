#!/bin/sh

# SSeeds - 12.12.23
# Contains all GMn timing sets accounting for timing shifts where they are manifest

# root --web=off -l 'simple_setalign.C("gmn",4,0,11589,11595,0,0,false,0,0,-249.448)'

# # Wait for user confirmation
# echo "Please check the output plots. Press any key to continue..."
# read -n 1 -s -r

root --web=off -l 'simple_setalign.C("gmn",4,0,0,0,0,0,false,0,0,-249.948)'

# Wait for user confirmation
echo "Please check the output plots. Press any key to continue..."
read -n 1 -s -r

root --web=off -l 'simple_setalign.C("gmn",7,0,0,0,0,0,false,0,0,-239.234)'

# Wait for user confirmation
echo "Please check the output plots. Press any key to continue..."
read -n 1 -s -r

root --web=off -l 'simple_setalign.C("gmn",11,1,12450,12860,0,0,false,0,0,-271.841)'

# Wait for user confirmation
echo "Please check the output plots. Press any key to continue..."
read -n 1 -s -r

root --web=off -l 'simple_setalign.C("gmn",11,1,0,0,12450,12860,false,0,0,-267.471)'

# Wait for user confirmation
echo "Please check the output plots. Press any key to continue..."
read -n 1 -s -r

root --web=off -l 'simple_setalign.C("gmn",14,1,13239,13260,0,0,false,0,0,-243.765)'

# Wait for user confirmation
echo "Please check the output plots. Press any key to continue..."
read -n 1 -s -r

root --web=off -l 'simple_setalign.C("gmn",14,1,13261,13407,0,0,false,0,0,-243.765)'

# Wait for user confirmation
echo "Please check the output plots. Press any key to continue..."
read -n 1 -s -r

root --web=off -l 'simple_setalign.C("gmn",8,1,0,0,0,0,false,0,0,-253.499)'

# Wait for user confirmation
echo "Please check the output plots. Press any key to continue..."
read -n 1 -s -r

root --web=off -l 'simple_setalign.C("gmn",9,1,0,0,0,0,false,0,0,-250.951)'

# Wait for user confirmation
echo "Please check the output plots. Press any key to continue..."
read -n 1 -s -r

