#!/bin/sh

script_path="/nas/longleaf/home/kolbt/whingdingdilly/simulation_camera"
sedtype='sed'

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

# See what parameter is changing

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

# You'll assign x,y position using this sorted array
echo "Variable is: ${var}"
echo "${sortVar[@]}"
len=${#sortVar[@]}
echo "Consists of ${len} simulations"

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
        rm ${image}

    done

done

# This will make a movie using the wall pngs, via ffmpeg
ffmpeg -framerate 10 -i wall_%04d.png\
 -vcodec libx264 -s 6000x4500 -pix_fmt yuv420p\
 -threads 1 wallOutLarge.mp4

ffmpeg -framerate 10 -i wall_%04d.png\
 -vcodec libx264 -s 2000x1500 -pix_fmt yuv420p\
 -threads 1 wallOut.mp4

ffmpeg -framerate 10 -i wall_%04d.png\
 -vcodec libx264 -s 1600x1200 -pix_fmt yuv420p\
 -threads 1 wallOutmedium.mp4

ffmpeg -framerate 10 -i wall_%04d.png\
 -vcodec libx264 -s 1000x750 -pix_fmt yuv420p\
 -threads 1 wallOutsmall.mp4
