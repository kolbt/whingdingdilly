#!/bin/bash

echo "Are you on longleaf?"
read answer

echo "Epsilon in filename? (y/n)"
read iseps

if [ $answer == 'y' ]; then
    echo "Clusters are dope"
    sedtype='sed'
elif [ $answer == 'n' ]; then
    echo "Classic desktop... laptop?"
    sedtype='gsed'
fi

mkdir "resubmits"

# Run from the directory where the slurm files are
if [ $iseps == "n" ]
then
    for file in $(ls slurm-*)
    do

        count=$(grep -c "run complete" $file)
        if [ $count -eq 1 ]; then
            # Grab the python file that you need to resubmit
            pa=$(grep -a -m 1 -h 'pa' $file | $sedtype 's/^.*pa\([0-9]*\)_.*/\1/')
            pb=$(grep -a -m 1 -h '_pb' $file | $sedtype 's/^.*_pb\([0-9]*\)_.*/\1/')
            xa=$(grep -a -m 1 -h '_xa' $file | $sedtype 's/^.*_xa\([0-9]*\)..*/\1/')
            pyfile=pa${pa}_pb${pb}_xa${xa}.py
            echo "Moving files... python: $pyfile   slurm: $file"
            mv $pyfile resubmits
            rm $file
        fi
    done

else
    for file in $(ls slurm-*)
    do

        count=$(grep -c "run complete" $file)
        if [ $count -eq 1 ]; then
            # Grab the python file that you need to resubmit
            pa=$(grep -a -m 1 -h 'pa' $file | $sedtype 's/^.*pa\([0-9]*\)_.*/\1/')
            pb=$(grep -a -m 1 -h '_pb' $file | $sedtype 's/^.*_pb\([0-9]*\)_.*/\1/')
            xa=$(grep -a -m 1 -h '_xa' $file | $sedtype 's/^.*_xa\([0-9]*\)_.*/\1/')
            ep=$(grep -a -m 1 -h '_eps' $file | $sedtype 's/^.*_eps\([0-9]*\)..*/\1/')
            pyfile=pa${pa}_pb${pb}_xa${xa}_eps${ep}.py
            echo "Moving python: $pyfile   Removing slurm: $file"
            mv $pyfile resubmits
            rm $file
        fi
    done

fi
