#!/bin/sh

script_path="/Users/kolbt/Desktop/compiled/whingdingdilly/simulation_camera"
sedtype='gsed'


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

    # Make png series for a simulation
    ovitos ${script_path}/pngSeries.py ${sim}
    # Get everything before the file extension
    pe=$(echo $sim | $sedtype 's/^.*pa\([0-9]*\)_.*/\1/')
    ep=$(echo $sim | $sedtype 's/^.*pb\([0-9]*\)_.*/\1/')
    phi=$(echo $sim | $sedtype 's/^.*xa\([0-9]*\)..*/\1/')
    # Make an individual ffmpeg movie
    ffmpeg -framerate 10 -i pa${pa}_pb${pb}_xa${xa}_ep${ep}_frame_%04d.png\
     -vcodec libx264 -s 2000x2000 -pix_fmt yuv420p\
     -threads 1 pa${pa}_pb${pb}_xa${xa}_ep${ep}.mp4

    # Get the png's index
    if [ ${var} == "pa" ]; then
        for iii in ${!sortVar[@]}; do
            if [ "${sortVar[$iii]}" == "${pa}" ]; then
                pos=$iii
            fi
        done

    elif [ ${var} == "pb" ]; then
        for iii in ${!sortVar[@]}; do
            if [ "${sortVar[$iii]}" == "${pb}" ]; then
                pos=$iii
            fi
        done

    elif [ ${var} == "xa" ]; then
        for iii in ${!sortVar[@]}; do
            if [ "${sortVar[$iii]}" == "${xa}" ]; then
                pos=$iii
            fi
        done

    elif [ ${var} == "ep" ]; then
        for iii in ${!sortVar[@]}; do
            if [ "${sortVar[$iii]}" == "${ep}" ]; then
                pos=$iii
        fi
    done

    fi

    # Place image on wall composite
    for image in $(ls pa${pa}_pb${pb}_xa${xa}*frame_*.png)
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

# Get rid of the source pngs
#rm wall_*.png
