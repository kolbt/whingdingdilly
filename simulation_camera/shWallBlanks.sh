#!/bin/sh

script_path="/nas/longleaf/home/kolbt/whingdingdilly/simulation_camera"
sedtype='sed'

# Place image on wall composite
for image in $(ls pe*frame_*.png)
do

    # Get frame number
    frame=$(echo $image | $sedtype 's/^.*frame_\([0-9]*\)..*/\1/')
    # Check if wall frame exists, if not, create it
    wallFrame="wall_${frame}.png"
    if [ ! -f ${wallFrame} ]; then
        sbatch ${script_path}/sbatchWallBlanks.sh ${wallFrame} ${script_path}
    fi

done
