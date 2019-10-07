#!/bin/sh

# Longleaf
script_path="/nas/longleaf/home/kolbt/whingdingdilly/simulation_camera"
sedtype='sed'
submit='sbatch'

## Desktop
#script_path="/Users/kolbt/Desktop/compiled/whingdingdilly/simulation_camera"
#sedtype='gsed'
#submit='sh'

# Get a specific name structure to use for frames
for sim in $(ls pe*gsd)
do

    # Get everything before the file extension
    pe=$(echo $sim | $sedtype 's/^.*pe\([0-9]*\)_.*/\1/')
    ep=$(echo $sim | $sedtype 's/^.*ep\([0-9]*\)_.*/\1/')
    phi=$(echo $sim | $sedtype 's/^.*phi\([0-9]*\)..*/\1/')

done

search="pe${pe}_ep${ep}_phi${phi}_frame_*.png"
for number in $(ls ${search})
do

    # Grab the frame number
    frame=$(echo $number | $sedtype 's/^.*frame_\([0-9]*\)..*/\1/')
    # Submit by frame number
    ${submit} ${script_path}/sbByFrame.sh ${frame} ${script_path}

done
