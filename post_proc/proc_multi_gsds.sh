#!/bin/sh

echo "Are you running on Longleaf (y/n)?"
read answer

if [ $answer == "y" ]; then
    hoomd_path='/nas/longleaf/home/kolbt/programs/hoomd-blue/build'
    gsd_path='/nas/longleaf/home/kolbt/programs/gsd/build'
    script_path='/nas/longleaf/home/kolbt/whingdingdilly/post_proc'
    sedtype='sed'
    submit='sbatch'
else
    hoomd_path='/Users/kolbt/Desktop/compiled/hoomd-blue_11_8_17/hoomd-blue/build'
    gsd_path='/Users/kolbt/Desktop/compiled/gsd/build'
    script_path='/Users/kolbt/Desktop/compiled/whingdingdilly/post_proc'
    sedtype='gsed'
    submit='sh'
fi

echo "Are you running on your laptop (y/n)?"
read answer

if [ $answer == "y" ]; then
    hoomd_path='/Users/kolbt/Desktop/compiled/hoomd-blue/build'
fi

gsdFiles=()

for filename in $( ls pa*.gsd )
do

    gsdFiles+=($filename)

done

$submit $script_path/analyze_multi_gsds.sh $script_path $gsd_path ${gsdFiles[@]}

