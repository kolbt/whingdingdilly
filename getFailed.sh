#!/bin/bash

echo "Are you on longleaf?"
read answer

if [ $answer == 'y' ]; then
    echo "Clusters are dope"
    sedtype='gsed'
elif [ $answer == 'n' ]; then
    echo "Classic desktop... laptop?"
    sedtype='sed'
fi

for file in $(ls slurm-*)
do

    count=$(grep -c "run complete" $file)
    if [ $count -eq 1 ]; then
        # Grab the python file that you need to resubmit
        pa=$(grep -a -m 1 -h 'pa' $file | $sedtype 's/^.*pa\([0-9]*\)_.*/\1/')
        pb=$(grep -a -m 1 -h '_pb' $file | $sedtype 's/^.*_pb\([0-9]*\)_.*/\1/')
        xa=$(grep -a -m 1 -h '_xa' $file | $sedtype 's/^.*_xa\([0-9]*\)..*/\1/')
        pyfile=pa${pa}_pb${pb}_xa${xa}.py
        echo "python: $pyfile   slurm: $file"
    fi

done
