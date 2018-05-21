#!/bin/sh

script_path="/Users/kolbt/Desktop/compiled/whingdingdilly/simulation_camera"
sedtype='gsed'

pa=()
pb=()
xa=()

# Read filenames for parameters
for sim in $(ls pa*.gsd)
do

    pa+=($(echo $sim | $sedtype 's/^.*pa\([0-9]*\)_.*/\1/'))
    pb+=($(echo $sim | $sedtype 's/^.*pb\([0-9]*\)_.*/\1/'))
    xa+=($(echo $sim | $sedtype 's/^.*xa\([0-9]*\)..*/\1/'))

done

# See what parameter is changing

# That IFS line came from the below discussion, thanks antak!
# https://stackoverflow.com/questions/7442417/how-to-sort-an-array-in-bash

if [ ${pa[0]} -ne ${pa[1]} ]; then
    var="pa"
    IFS=$'\n' sortVar=($(sort -g <<<"${pa[*]}"))
elif [ ${pb[0]} -ne ${pb[1]} ]; then
    var="pb"
    IFS=$'\n' sortVar=($(sort -g <<<"${pb[*]}"))
elif [ ${xa[0]} -ne ${xa[1]} ]; then
    var="xa"
    IFS=$'\n' sortVar=($(sort -g <<<"${xa[*]}"))
fi

# You'll assign x,y position using this sorted array
echo "Variable is: ${var}"
echo "${sortVar[@]}"
len=${#sortVar[@]}
echo "Consists of ${len} simulations"

for sim in $(ls pa*.gsd)
do

    # Make png series for a simulation
    ovitos ${script_path}/pngSeries.py ${sim}
    # Get everything before the file extension
    pa=$(echo $sim | $sedtype 's/^.*pa\([0-9]*\)_.*/\1/')
    pb=$(echo $sim | $sedtype 's/^.*pb\([0-9]*\)_.*/\1/')
    xa=$(echo $sim | $sedtype 's/^.*xa\([0-9]*\)..*/\1/')

    # Get the png's index
    if [ ${var} == "pa" ]; then
        for iii in ${!sortVar[@]}; do
            if [ "${sortVar[$iii]}" == "${pa}" ]; then
                pos=$iii
            fi
        done
    elif [ ${var} == "pb" ]; then
        for iii in {sortVar[@]}; do
            if [ ${sortVar[$iii]} == ${pb} ]; then
                pos=$iii
            fi
        done

    elif [ ${var} == "xa" ]; then
        for iii in {sortVar[@]}; do
            if [ ${sortVar[$iii]} == ${xa} ]; then
                pos=$iii
            fi
        done

    fi

    # Place image on wall composite
    for image in $(ls pa${pa}_pb${pb}_xa${xa}_frame_*.png)
    do

        # Get frame number
        frame=$(echo $image | $sedtype 's/^.*frame_\([0-9]*\)..*/\1/')
        # Check if wall frame exists, if not, create it
        wallFrame="wall_${frame}_pa${pa}_pb${pb}_xa${xa}.png"
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
