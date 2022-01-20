#!bin/bash

echo -e "\n"
echo "You've arrived at the center for HCal expert analysis. Setting path variables for local scripts.." 
echo -e "\n"

source $HOME/sbs/sbs_devel/install/bin/sbsenv.sh

export L_DIR=/adaqfs/home/a-onl/sbs/HCal_replay/replay
export DATA_DIR=/adaqeb1/data1
export OUT_DIR=/adaqfs/home/a-onl/sbs/Rootfiles
#export DB_DIR=/adaqfs/home/a-onl/sbs/HCal_replay/replay/DB
export DB_DIR=/adaqfs/home/a-onl/sbs/sbs_devel/SBS-replay/DB
