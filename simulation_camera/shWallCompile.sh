#!/bin/sh

script_path="/nas/longleaf/home/kolbt/whingdingdilly/simulation_camera"
#script_path="/Users/kolbt/Desktop/compiled/whingdingdilly/simulation_camera"
sedtype='sed'
#sedtype='gsed'

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

    sbatch ${script_path}/sbatchWallCompile.sh ${pe} ${ep} ${phi} ${pos} ${script_path}
#    sh ${script_path}/sbatchWallCompile.sh ${pe} ${ep} ${phi} ${pos} ${script_path}

done
