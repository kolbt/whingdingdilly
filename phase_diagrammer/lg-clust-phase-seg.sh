#!/bin/sh

echo "Longleaf?"
read longleaf

txtFiles=()

if [ $longleaf == 'y' ]; then
    script='/nas/longleaf/home/kolbt/whingdingdilly/phase_diagrammer'
    submit='sbatch'
else
    script='/Users/kolbt/Desktop/compiled/whingdingdilly/phase_diagrammer'
    submit='sh'
fi

# Get all text data
for txt in $(ls diam*.txt)
do

    # Add to the array that you'll pass to the py file
    txtFiles+=($txt)

done

$submit $script/txt-to-py.sh $script ${txtFiles[@]}
