#!/bin/sh

# SSeeds - 02.03.24
# Contains all GMn timing sets accounting for timing shifts where they are manifest

#Visually verified exclusion range SBS11

root --web=off -l 'simple_setalign.C("gmn",11,1,12450,12860,0,0,false,0,0,-271.841)'

# Wait for user confirmation
echo "Please check the output plots. Press any key to continue..."
read -n 1 -s -r

root --web=off -l 'simple_setalign.C("gmn",11,1,0,0,12450,12860,false,0,0,-267.471)'

# Wait for user confirmation
echo "Please check the output plots. Press any key to continue..."
read -n 1 -s -r

#Visually verified exclusion range SBS14

root --web=off -l 'simple_setalign.C("gmn",14,1,13239,13260,0,0,false,0,0,-243.765)'

# Wait for user confirmation
echo "Please check the output plots. Press any key to continue..."
read -n 1 -s -r

root --web=off -l 'simple_setalign.C("gmn",14,1,13261,13407,0,0,false,0,0,-243.765)'

# Wait for user confirmation
echo "Please check the output plots. Press any key to continue..."
read -n 1 -s -r


