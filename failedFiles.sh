#!/bin/sh

string1="Error"

for infile in $(ls slurm-*);
do
    # The slurm file is a failed file
	if tail -n 1 $infile | grep -Eq $string1;then
        # Get the python filename
        grep -o -m1 'pe[^"]*.py' $infile
    fi
done
