#!/bin/sh

count=1
dircount=1

mkdir -p "subdir_$dircount" # sort into batches

for filename in $( ls *py ) # look for all infiles.py
do

    if [ $(( $count % 5 )) == 0 ]; then     # sort into batches of 5
        dircount=$(( $dircount + 1 ))       # increase directory number
        mkdir -p "subdir_$dircount"         # make the batch directory
    fi
        count=$(( $count + 1 ))             # increment the counter
        mv $filename subdir_$dircount       # move the file for batch submission

done
