#!/bin/sh

echo "Are you running on Longleaf (y/n)?"
read answer

if [ $answer == "y" ]; then
    script_path='/nas/longleaf/home/kolbt/whingdingdilly/post_proc'
    sedtype='sed'
    submit='sbatch'
else
    script_path='/Users/kolbt/Desktop/compiled/whingdingdilly/post_proc'
    sedtype='gsed'
    submit='sh'
fi

echo "Are you running on your laptop (y/n)?"
read answer

if [ $answer == "y" ]; then
    hoomd_path='/Users/kolbt/Desktop/compiled/hoomd-blue/build'
fi

txtFiles=()

for filename in $( ls *pa*.txt )
do

    txtFiles+=($filename)

done

$submit $script_path/analyzeTxt.sh $script_path ${txtFiles[@]}

