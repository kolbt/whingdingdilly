#!/bin/sh

# Grab computed filenames
string1="complete"
pyFiles=()

for infile in $(ls slurm-*);
do
    # The slurm file is a completed file
    if tail -n 1 $infile | grep -Eq $string1;then
        # Get the python filename, add to list
        pyFiles+=( "$(grep -o -m1 'pe[^"]*.py' $infile)" )
    fi
done

# Loop through the completed python files
for i in ${pyFiles[@]}
do
    # Grab parameters
    inFile=$i
    echo "$i"
    echo $inFile
    # Submit for analysis
    

done
