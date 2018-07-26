#!/bin/sh

script_path="/Users/kolbt/Desktop/compiled/whingdingdilly/simulation_camera"
sedtype='gsed'

echo "Varying epsilon directly? (Y/n)"
read answer

pa=()
pb=()
xa=()

if [ ${answer} == "n" ]; then
    # Read filenames for parameters
    for sim in $(ls pa*.gsd)
    do

        pa+=($(echo $sim | $sedtype 's/^.*pa\([0-9]*\)_.*/\1/'))
        pb+=($(echo $sim | $sedtype 's/^.*pb\([0-9]*\)_.*/\1/'))
        xa+=($(echo $sim | $sedtype 's/^.*xa\([0-9]*\)..*/\1/'))

    done

elif [ ${answer} == "Y" ]; then
    ep=()
    # Read filenames for parameters
    for sim in $(ls pa*.gsd)
    do

        pa+=($(echo $sim | $sedtype 's/^.*pa\([0-9]*\)_.*/\1/'))
        pb+=($(echo $sim | $sedtype 's/^.*pb\([0-9]*\)_.*/\1/'))
        xa+=($(echo $sim | $sedtype 's/^.*xa\([0-9]*\)_.*/\1/'))
        ep+=($(echo $sim | $sedtype 's/^.*ep\([0-9]*\)..*/\1/'))

    done
fi

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
elif [ ${ep[0]} -ne ${ep[1]} ]; then
    var="ep"
    IFS=$'\n' sortVar=($(sort -g <<<"${ep[*]}"))
fi

# You'll assign x,y position using this sorted array
echo "Variable is: ${var}"
echo "${sortVar[@]}"
len=${#sortVar[@]}
echo "Consists of ${len} simulations"

for sim in $(ls pa*.gsd)
do

    if [ ${answer} == "n" ]; then
        # Get everything before the file extension
        pa=$(echo $sim | $sedtype 's/^.*pa\([0-9]*\)_.*/\1/')
        pb=$(echo $sim | $sedtype 's/^.*pb\([0-9]*\)_.*/\1/')
        xa=$(echo $sim | $sedtype 's/^.*xa\([0-9]*\)..*/\1/')
    elif [ ${answer} == "Y" ]; then
        pa=$(echo $sim | $sedtype 's/^.*pa\([0-9]*\)_.*/\1/')
        pb=$(echo $sim | $sedtype 's/^.*pb\([0-9]*\)_.*/\1/')
        xa=$(echo $sim | $sedtype 's/^.*xa\([0-9]*\)_.*/\1/')
        ep=$(echo $sim | $sedtype 's/^.*ep\([0-9]*\)..*/\1/')
    fi

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
    image=effSigma_pa${pa}_pb${pb}_xa${xa}_ep${ep}.png

    # Check if wall frame exists, if not, create it
    wallFrame="wall.png"
    if [ ! -f ${wallFrame} ]; then
        python ${script_path}/makeWallFrame.py ${wallFrame}
    fi
    # Add each png to corresponding wall frame
    python ${script_path}/placeStillFrame.py ${image} ${wallFrame} ${pos}

done
