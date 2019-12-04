#!/bin/sh

script_path='/nas/longleaf/home/kolbt/whingdingdilly/case_studies/full-video-analysis.py'
sedtype='sed'

for filename in $( ls *.gsd )
do
    
    # Input file is always of the form pe${peInner}_${peOuter}
    pein=$(echo $filename | $sedtype "s/^.*pe\([0-9]*\)_.*/\1/")
    peout=$(echo $filename | $sedtype "s/^.*pe${pein}_\([0-9]*\)_.*/\1/")
    sh /nas/longleaf/home/kolbt/whingdingdilly/case_studies/sbfullVideoAnalysis.sh $filename $script_path $pein $peout
    
done
