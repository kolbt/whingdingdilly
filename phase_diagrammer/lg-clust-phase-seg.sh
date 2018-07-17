#!/bin/sh

txtFiles=()
script='/nas/longleaf/home/kolbt/whingdingdilly/phase_diagrammer/txt-to-py.sh'

# Get all text data
for txt in $(ls *.txt)
do

    # Add to the array that you'll pass to the py file
    txtFiles+=($txt)

done

sbatch $script ${txtFiles[@]}
