#!bin/sh

# sseeds 10.13.22

## Usage
# ./run_hcal_gen_replay.sh <run>

echo -e "\n"
echo "Replaying all streams from run $1.." 
echo -e "\n"

#For three streams
for (( i=0; i<3; i++ ))
do
    analyzer -l -b -q 'replay_hcal_GEn.C('$1',-1,0,0,'$i')'
done

echo -e "\n"
echo "Replay over run $1 complete." 
echo -e "\n"
