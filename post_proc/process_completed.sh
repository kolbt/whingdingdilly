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

echo "Are you running on Longleaf (y/n)?"
read answer

if [ $answer == "y" ]; then
    script_path='/nas/longleaf/home/kolbt/whingdingdilly/post_proc'
    sedtype='sed'
    submit='sbatch'
    module add python/3.6.6
else
    script_path='/Users/kolbt/Desktop/compiled/whingdingdilly/post_proc'
    sedtype='gsed'
    submit='sh'
fi

# Loop through the completed python files
for i in ${pyFiles[@]}
do
    # Grab parameters
    pe=$(echo $i | $sedtype 's/^.*pe\([0-9]*\)_.*/\1/')
    phi=$(echo $i | $sedtype 's/^.*phi\([0-9]*\)..*/\1/')
    # Submit for analysis
    for j in $(ls cluster*pe${}_phi${phi}*)
    do
        echo $j
    done
    
done
