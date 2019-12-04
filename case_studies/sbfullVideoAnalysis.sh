#!/bin/sh
#SBATCH -p general                          # partition to run on
#SBATCH -n 1                                # number of cores
#SBATCH -t 3-00:00                          # time (D-HH:MM)

# Command to increase memory allocated --mem=100g

filename=$1
script_path=$2
pein=$3
peout=$4

python $script_path $filename $pein $peout

ffmpeg -framerate 10 -i peIn${pein}_peOut${peout}_fm%04d.jpg\
 -vcodec libx264 -s 2596x1498 -pix_fmt yuv420p -threads 1\
 peIn${pein}_peOut${peout}.mp4
