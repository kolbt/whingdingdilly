#!/bin/sh
#SBATCH -p general                          # partition to run on
#SBATCH -n 1                                # number of cores
#SBATCH -t 3-00:00                         # time (D-HH:MM)
#SBATCH -o placeFrames.out

# Command to increase memory allocated --mem=100g

pe=$1
ep=$2
phi=$3
pos=$4
script_path=$5

sedtype='sed'
#sedtype='gsed'

# Place image on wall composite
for image in $(ls pe${pe}_ep${ep}_phi${phi}*frame_*.png)
do

    # Get frame number
    frame=$(echo $image | $sedtype 's/^.*frame_\([0-9]*\)..*/\1/')
    # Check if wall frame exists, if not, create it
    wallFrame="wall_${frame}.png"
    if [ ! -f ${wallFrame} ]; then
        python ${script_path}/makeWallFrame.py ${wallFrame}
    fi
    # Add each png to corresponding wall frame
    python ${script_path}/placeMovieFrame.py ${image} ${wallFrame} ${pos}
    # Remove pngs after adding
#    rm ${image}

done


