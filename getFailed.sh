#!/bin/bash

for file in $(ls slurm-*)
do

    count = $(( grep -c "run complete" $file))
    if [ $count -eq 1 ]; then
        echo $file
    fi

done
