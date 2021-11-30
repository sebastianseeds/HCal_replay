#!/bin/sh

## Edits
# P. Datta <pdbforce@jlab.org> Created 11-14-2021

## Usage
# This macro performs the analysis of cosmic data 
# for BB Shower and PreShower. Example execution: 
#./run_cosmic_analysis.sh <nrun> <nevent>

cd macros
root -l 'PreShower_macros/bbps_cos_cal.C('$1','$2',0)'
root -l 'Shower_macros/bbsh_cos_cal.C('$1','$2',0)'

