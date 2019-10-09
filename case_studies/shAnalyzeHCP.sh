#!/bin/sh

echo "Are you running on Longleaf (y/n)?"
read answer

if [ $answer == "y" ]; then
    script_path='/nas/longleaf/home/kolbt/whingdingdilly/case_studies'
    sedtype='sed'
    submit='sbatch'
else
    script_path='/Users/kolbt/Desktop/compiled/whingdingdilly/case_studies'
    sedtype='gsed'
    submit='sh'
fi

searchString="active_"

for filename in $( ls *.gsd )
do
    
    # Check if simulation is ballistic or active
    if [[ $filename == ${searchString}* ]]; then
        ball='n'
    else
        ball='y'
    fi
    
    pe=$(echo $filename | $sedtype 's/^.*[^0-9]\([0-9]*\.[0-9]*\)_lattice.*$/\1/')
    lat=$(echo $filename | $sedtype 's/^.*[^0-9]\([0-9]*\.[0-9]*\).gsd.*$/\1/')
    $submit $script_path/sbAnalyzeHCP.sh ${ball} ${pe} ${lat} ${script_path}
    
done


