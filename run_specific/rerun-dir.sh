#!/bin/sh

echo "Are you running on Longleaf (y/n)?"
read answer

if [ $answer == "y" ]; then
    script_path='/nas/longleaf/home/kolbt/whingdingdilly/run.sh'
    submit='sbatch'
else
    script_path='/Users/kolbt/Desktop/compiled/whingdingdilly/run.sh'
    submit='sh'
fi

echo "GPU (y/n)?"
read gpu

if [ $gpu == "y" ]; then
    script_path='/nas/longleaf/home/kolbt/whingdingdilly/run_gpu.sh'
fi

# Run from the directory where the slurm files are
for file in $(ls *.py)
do

    $submit $script_path $file

done





