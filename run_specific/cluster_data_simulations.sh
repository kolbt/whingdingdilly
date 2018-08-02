#!/bin/sh

current=$( date "+%m_%d_%y" )
this_path=$( pwd )

echo "Are you running on Longleaf (y/n)?"
read answer

if [ $answer == "y" ]; then
    hoomd_path='/nas/longleaf/home/kolbt/programs/cpu-hoomd/hoomd-blue/build'
    gsd_path='/nas/longleaf/home/kolbt/programs/gsd/build'
    script_path='/nas/longleaf/home/kolbt/whingdingdilly/run.sh'
    sedtype='sed'
    submit='sbatch'
else
    hoomd_path='/Users/kolbt/Desktop/compiled/hoomd-blue_11_8_17/hoomd-blue/build'
    gsd_path='/Users/kolbt/Desktop/compiled/gsd/build'
    script_path='/Users/kolbt/Desktop/compiled/whingdingdilly/run.sh'
    sedtype='gsed'
    submit='sh'
fi

echo "GPU (y/n)?"
read gpu

if [ $gpu == "y" ]; then
    hoomd_path='/nas/longleaf/home/kolbt/programs/hoomd_2.2.1/hoomd-blue/build'
    script_path='/nas/longleaf/home/kolbt/whingdingdilly/run_gpu.sh'
fi

echo "Is epsilon in the filename (y/n)?"
read isEps

mkdir ${current}_parent
cd ${current}_parent

for gsd in $(ls pa*.gsd)
do
    # Reads parameters from filename
    if [ $isEps == "n" ]; then
        pa=($(echo $gsd | $sedtype 's/^.*pa\([0-9]*\)_.*/\1/'))
        pb=($(echo $gsd | $sedtype 's/^.*pb\([0-9]*\)_.*/\1/'))
        xa=($(echo $gsd | $sedtype 's/^.*xa\([0-9]*\)..*/\1/'))
        ep=1    # default epsilon behavior
        pyFile="pa${pa}_pb${pb}_xa${xa}.py"
    else
        pa=($(echo $gsd | $sedtype 's/^.*pa\([0-9]*\)_.*/\1/'))
        pb=($(echo $gsd | $sedtype 's/^.*pb\([0-9]*\)_.*/\1/'))
        xa=($(echo $gsd | $sedtype 's/^.*xa\([0-9]*\)_.*/\1/'))
        ep=($(echo $gsd | $sedtype 's/^.*ep\([0-9]*\)..*/\1/'))
        pyFile="pa${pa}_pb${pb}_xa${xa}_eps{$ep}.py"
    fi

    # Needs to read seeds from old python file
    seed1=($(cat $gsd | $sedtype 's/^.*seed1 = \([0-9]*\)_.*/\1/'))


    # Things that need to be passed to the run.sh file:
    # 1.    Path to HOOMd
    # 2.    Path to GSD
    # 3.    Path to scripts
    # 4-7.  Params: pa, pb, xa, eps
    # 8-12. Seeds
    # 13.   Frame to load
    $submit $script_path $infile
done
