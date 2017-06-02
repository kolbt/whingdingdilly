#!/bin/sh
#SBATCH -p general                          # partition to run on
#SBATCH -n 1                                # number of cores
#SBATCH -t 0-10:00                          # time (D-HH:MM)


filename=$1
pa=$2
pb=$3
xa=$4

python /Users/kolbt/Desktop/compiled/whingdingdilly/post_proc/post_proc.py $filename $pa $pb $xa # run infiles with python
