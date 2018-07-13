#!/bin/sh

current=$( date "+%m_%d_%y" )
this_path=$( pwd )

echo "Are you running on Longleaf (y/n)?"
read answer

if [ $answer == "y" ]; then
    hoomd_path='/nas/longleaf/home/kolbt/programs/cpu-hoomd/hoomd-blue/build'
    gsd_path='/nas/longleaf/home/kolbt/programs/gsd/build'
    script_path='/nas/longleaf/home/kolbt/whingdingdilly/run.sh'
    template='/nas/longleaf/home/kolbt/whingdingdilly/run_specific/vary_v_and_eps.py'
    sedtype='sed'
    submit='sbatch'
else
    hoomd_path='/Users/kolbt/Desktop/compiled/hoomd-blue_11_8_17/hoomd-blue/build'
    gsd_path='/Users/kolbt/Desktop/compiled/gsd/build'
    script_path='/Users/kolbt/Desktop/compiled/whingdingdilly/run.sh'
    template='/Users/kolbt/Desktop/compiled/whingdingdilly/run_specific/vary_v_and_eps.py'
    sedtype='gsed'
    submit='sh'
fi

echo "GPU (y/n)?"
read gpu

if [ $gpu == "y" ]; then
    hoomd_path='/nas/longleaf/home/kolbt/programs/hoomd_2.2.1/hoomd-blue/build'
    script_path='/nas/longleaf/home/kolbt/whingdingdilly/run_gpu.sh'
fi

echo "Epsilon in filename? (y/n)"
read iseps

# Run from the directory where the slurm files are
if [ $iseps == "n" ]; then
    for file in $(ls slurm-*)
    do

        count=$(grep -c "run complete" $file)
        if [ $count -eq 1 ]; then
            # Grab the python file that you need to resubmit
            pa=$(grep -a -m 1 -h 'pa' $file | $sedtype 's/^.*pa\([0-9]*\)_.*/\1/')
            pb=$(grep -a -m 1 -h '_pb' $file | $sedtype 's/^.*_pb\([0-9]*\)_.*/\1/')
            xa=$(grep -a -m 1 -h '_xa' $file | $sedtype 's/^.*_xa\([0-9]*\)..*/\1/')
            infile=pa${pa}_pb${pb}_xa${xa}.py
            $submit $script_path $infile
            rm $file
        fi

    done
fi

# If epsilon is in the filename
if [ $iseps == "y" ]; then
    for file in $(ls slurm-*)
    do

        count=$(grep -c "run complete" $file)
        if [ $count -eq 1 ]; then
            # Grab the python file that you need to resubmit
            pa=$(grep -a -m 1 -h 'pa' $file | $sedtype 's/^.*pa\([0-9]*\)_.*/\1/')
            pb=$(grep -a -m 1 -h '_pb' $file | $sedtype 's/^.*_pb\([0-9]*\)_.*/\1/')
            xa=$(grep -a -m 1 -h '_xa' $file | $sedtype 's/^.*_xa\([0-9]*\)_.*/\1/')
            ep=$(grep -a -m 1 -h '_eps' $file | $sedtype 's/^.*_eps\([0-9]*\)..*/\1/')
            infile=pa${pa}_pb${pb}_xa${xa}_eps${ep}.py
            $submit $script_path $infile
            rm $file
        fi

    done
fi


