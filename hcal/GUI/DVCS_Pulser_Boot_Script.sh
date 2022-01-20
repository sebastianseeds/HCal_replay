#!/bin/bash
input="/home/daq/test_fadc/DVCS_Pulser_Settings.json"
#input="DVCS_Pulser_Settings_Log_File.txt"

#Load the library.
remex -d hcalROC17 linuxvme/dvcsPulser/libdvcsPulser.so

#Initialize the DVCS pulser.
remex hcalROC17 'dvcsPulserInit(0xEF0,0)'

#Enable the output channels.
remex hcalROC17 'dvcsPulserEnable(0x1d)'

#Read in last frequency setting.
frequency=$(awk '{print $1}' $input)
#frequency=$20
#Set the frequency to previous setting.
#remex hcalROC17 'dvcsPulserSetFrequency(50)'
remex hcalROC17 'dvcsPulserSetFrequency('$frequency')'

echo "$frequency"
