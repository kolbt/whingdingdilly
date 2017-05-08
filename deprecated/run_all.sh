#!/bin/sh

script_path=$1
answer=$2

if [ $answer == "y" ]; then
    for foldername in $( ls )   # loop through all batch dirs
    do

        cd "$foldername"                    # enter each batch dir
        sbatch $script_path/run.sh          # submit your job (specifics in run.sh)
        cd ..                               # leave the batch dir

    done

else
    for foldername in $( ls )   # loop through all batch dirs
    do

        cd "$foldername"                    # enter each batch dir
        sh $script_path/run.sh              # run locally
        cd ..                               # leave the batch dir

    done

fi
