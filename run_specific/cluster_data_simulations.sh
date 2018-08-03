#!/bin/sh

current=$( date "+%m_%d_%y" )
this_path=$( pwd )

echo "Are you running on Longleaf (y/n)?"
read answer

if [ $answer == "y" ]; then
    hoomdPath='/nas/longleaf/home/kolbt/programs/cpu-hoomd/hoomd-blue/build'
    gsdPath='/nas/longleaf/home/kolbt/programs/gsd/build'
    script_path='/nas/longleaf/home/kolbt/whingdingdilly/run_specific/run_nuc.sh'
    sedtype='sed'
    submit='sbatch'
    inFile='/nas/longleaf/home/kolbt/whingdingdilly/run_specific/load_gsd.py'
else
    hoomdPath='/Users/kolbt/Desktop/compiled/hoomd-blue_11_8_17/hoomd-blue/build'
    gsdPath='/Users/kolbt/Desktop/compiled/gsd/build'
    script_path='/Users/kolbt/Desktop/compiled/whingdingdilly/run_specific/run_nuc.sh'
    sedtype='gsed'
    submit='sh'
    inFile='/Users/kolbt/Desktop/compiled/whingdingdilly/run_specific/load_gsd.py'
fi

echo "GPU (y/n)?"
read gpu

if [ $gpu == "y" ]; then
    hoomdPath='/nas/longleaf/home/kolbt/programs/hoomd_2.2.1/hoomd-blue/build'
    script_path='/nas/longleaf/home/kolbt/whingdingdilly/run_specific/run_gpu_nuc.sh'
fi

echo "Is epsilon in the filename (y/n)?"
read isEps

#mkdir ${current}_parent
#cd ${current}_parent

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
    seed1=($($sedtype -n -e 's/^seed1 = \([0-9]*\) .*$/\1/p' $pyFile))
    seed2=($($sedtype -n -e 's/^seed2 = \([0-9]*\) .*$/\1/p' $pyFile))
    seed3=($($sedtype -n -e 's/^seed3 = \([0-9]*\) .*$/\1/p' $pyFile))
    seed4=($($sedtype -n -e 's/^seed4 = \([0-9]*\) .*$/\1/p' $pyFile))
    seed5=($($sedtype -n -e 's/^seed5 = \([0-9]*\) .*$/\1/p' $pyFile))

    echo "What frame should be loaded?"
    read myFrame

    # Things that need to be passed to the run.sh file:
    # 1.    Path to HOOMd               GOT IT
    # 2.    Path to GSD                 DONE
    # 3.    Path to scripts             WORD UP
    # 4-7.  Params: pa, pb, xa, eps     PIECE OF CAKE
    # 8-12. Seeds                       NAILED IT
    # 13.   Frame to load               JUST ASK

    # Submit the nucleation simulation
    $submit $script_path $inFile $hoomdPath $gsdPath $pa $pb $xa $ep $seed1 $seed2 $seed3 $seed4 $seed5 $myFrame

done
