#!bin/bash

# sseeds 10.15.22

## Usage
# ./kin_gen.sh <kinematic>

echo -e "\n Generating expected elastic parameters for GEn kinematic $1.. \n"

if [ $1 == 1 ]
then
root -l -q -b 'GEnElasPeak.C(2.1,47.5,1.63,1.7,34.7,17)'

elif [ $1 == 2 ]
then
root -l -q -b 'GEnElasPeak.C(4.291,29.5,1.63,2.9,34.7,17)'

elif [ $1 == 3 ]
then
root -l -q -b 'GEnElasPeak.C(6.379,35.9,1.63,6.6,22.1,17)'

elif [ $1 == 4 ]
then
root -l -q -b 'GEnElasPeak.C(8.4,35.0,1.63,9.7,18.0,17)'

else
    echo -e "\n Invalid GEn kinematic entered. Enter a number between 1 and 4 and try again."
fi
