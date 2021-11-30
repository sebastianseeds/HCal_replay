#!/bin/sh

## Edits
# P. Datta <pdbforce@jlab.org> Created 11-14-2021

## Usage
# This script generates calibrated HVs for BB Shower and
# PreShower using cosmic data. Example execution: 
#./get_calibrated_hv.sh <nrun> <desired_trigger_amplitude>

cd macros
echo -e "\n"

# Calculating preshower HVs
echo -e "Getting PreShwoer HVs \n"
sleep 2
root -l -b -q 'PreShower_macros/ps_HVUpdate_cosmic.C('$1','$2')'
evince plots/ps_hv_calib_run_$1_$2mV.pdf
echo -e "\n"

sleep 2
# Calculating shower HVs
echo -e "Getting Shwoer HVs \n"
sleep 2
root -l -b -q 'Shower_macros/sh_HVUpdate_cosmic.C('$1','$2')'
evince plots/sh_hv_calib_run_$1_$2mV.pdf
echo -e "\n"

sleep 2
# Combine the HVs to get new settings
echo -e "Combining Shower and PreShower HVs \n"
sleep 2
root -l -b -q 'Combined_macros/Combine_HV.C(0,"sh_hv_calib_run_'$1'_'$2'mV.set","ps_hv_calib_run_'$1'_'$2'mV.set","hv_calibrated_run_'$1'_'$2'mV")'
evince plots/hv_calibrated_run_$1_$2mV.pdf &
