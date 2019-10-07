#!/bin/sh
#SBATCH -p general                          # partition to run on
#SBATCH -n 1                                # number of cores
#SBATCH -t 3-00:00                         # time (D-HH:MM)
#SBATCH -o placeFrames.out

frame=$1
script_path=$2
sedtype='sed'

# Name of wall composite
wallFrame="wall_${frame}.png"
# First make the frame
python ${script_path}/makeWallFrame.py ${wallFrame}

# Find positions by variable
pe=()
ep=()
phi=()

# Read filenames for parameters
for sim in $(ls pe*.gsd)
do

    pe+=($(echo $sim | $sedtype 's/^.*pe\([0-9]*\)_.*/\1/'))
    ep+=($(echo $sim | $sedtype 's/^.*ep\([0-9]*\)_.*/\1/'))
    phi+=($(echo $sim | $sedtype 's/^.*phi\([0-9]*\)..*/\1/'))

done

# That IFS line came from the below discussion, thanks antak!
# https://stackoverflow.com/questions/7442417/how-to-sort-an-array-in-bash

if [ ${pe[0]} -ne ${pe[1]} ]; then
    var="pe"
    IFS=$'\n' sortVar=($(sort -g <<<"${pe[*]}"))
elif [ ${ep[0]} -ne ${ep[1]} ]; then
    var="ep"
    IFS=$'\n' sortVar=($(sort -g <<<"${ep[*]}"))
elif [ ${phi[0]} -ne ${phi[1]} ]; then
    var="phi"
    IFS=$'\n' sortVar=($(sort -g <<<"${phi[*]}"))
fi

## You'll assign x,y position using this sorted array
#echo "Variable is: ${var}"
#echo "${sortVar[@]}"
#len=${#sortVar[@]}
#echo "Consists of ${len} simulations"

for sim in $(ls pe*.gsd)
do

    # Get everything before the file extension
    pe=$(echo $sim | $sedtype 's/^.*pe\([0-9]*\)_.*/\1/')
    ep=$(echo $sim | $sedtype 's/^.*ep\([0-9]*\)_.*/\1/')
    phi=$(echo $sim | $sedtype 's/^.*phi\([0-9]*\)..*/\1/')
    
    # Get the png's index
    if [ ${var} == "pe" ]; then
        for iii in ${!sortVar[@]}; do
            if [ "${sortVar[$iii]}" == "${pe}" ]; then
                pos=$iii
            fi
        done

    elif [ ${var} == "ep" ]; then
        for iii in ${!sortVar[@]}; do
            if [ "${sortVar[$iii]}" == "${ep}" ]; then
                pos=$iii
            fi
        done

    elif [ ${var} == "phi" ]; then
        for iii in ${!sortVar[@]}; do
            if [ "${sortVar[$iii]}" == "${phi}" ]; then
                pos=$iii
            fi
        done
        
    fi

    # Name of component simulation png
    image="pe${pe}_ep${ep}_phi${phi}_frame_${frame}.png"
    # Add all pngs for individual frame
    python ${script_path}/placeInFrame.py ${image} ${wallFrame} ${pos}

done
