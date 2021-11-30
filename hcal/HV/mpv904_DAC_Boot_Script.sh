#!/bin/bash
input="/home/daq/test_fadc/mpv904_DAC_Voltage_Set_Values.txt"
#input="mpv904_DAC_Voltage_Settings_Log_File.txt"

let "min_volt=-10"
let "max_volt=10"
let "max_value=4095"
let "volt_range=($max_volt-$min_volt)"
#echo $max_value
#echo $volt_range

#Load the library.
remex -d hcalROC17 linuxvme/mpv904/libmpv904.so

#Initialize the mpv904 DAC.
remex hcalROC17 'mpv904Init(0x700000,0,1,0x6)'

#Read in the previous set values from the txt file (no floating point in bash).
#awk '{ print $1}' $input
f0=$(awk '{ print $1}' $input)
f1=$(awk '{ print $2}' $input)
f2=$(awk '{ print $3}' $input)
f3=$(awk '{ print $4}' $input)
f4=$(awk '{ print $5}' $input)
f5=$(awk '{ print $6}' $input)
f6=$(awk '{ print $7}' $input)
f7=$(awk '{ print $8}' $input)
f8=$(awk '{ print $9}' $input)
f9=$(awk '{ print $10}' $input)
f10=$(awk '{ print $11}' $input)
f11=$(awk '{ print $12}' $input)
f12=$(awk '{ print $13}' $input)
f13=$(awk '{ print $14}' $input)
f14=$(awk '{ print $15}' $input)
f15=$(awk '{ print $16}' $input)

#Set the DAC voltage outputs for each channel.
remex hcalROC17 'mpv904SetDAC(0,0,'$f0')'
remex hcalROC17 'mpv904SetDAC(0,1,'$f1')'
remex hcalROC17 'mpv904SetDAC(0,2,'$f2')'
remex hcalROC17 'mpv904SetDAC(0,3,'$f3')'
remex hcalROC17 'mpv904SetDAC(0,4,'$f4')'
remex hcalROC17 'mpv904SetDAC(0,5,'$f5')'
remex hcalROC17 'mpv904SetDAC(0,6,'$f6')'
remex hcalROC17 'mpv904SetDAC(0,7,'$f7')'
remex hcalROC17 'mpv904SetDAC(0,8,'$f8')'
remex hcalROC17 'mpv904SetDAC(0,9,'$f9')'
remex hcalROC17 'mpv904SetDAC(0,10,'$f10')'
remex hcalROC17 'mpv904SetDAC(0,11,'$f11')'
remex hcalROC17 'mpv904SetDAC(0,12,'$f12')'
remex hcalROC17 'mpv904SetDAC(0,13,'$f13')'
remex hcalROC17 'mpv904SetDAC(0,14,'$f14')'
remex hcalROC17 'mpv904SetDAC(0,15,'$f15')'

#echo "$f0"
#echo "$f1"
