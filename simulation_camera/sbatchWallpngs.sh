#!/bin/sh
#SBATCH -p general                          # partition to run on
#SBATCH -n 1                                # number of cores
#SBATCH -t 11-00:00                         # time (D-HH:MM)

# Command to increase memory allocated --mem=100g

sim=$1
script_path=$2

python ${script_path}/movieWSigma.py ${sim}

# Get everything before the file extension
pe=$(echo $sim | sed 's/^.*pe\([0-9]*\)_.*/\1/')
ep=$(echo $sim | sed 's/^.*ep\([0-9]*\)_.*/\1/')
phi=$(echo $sim | sed 's/^.*phi\([0-9]*\)..*/\1/')
# Make an individual ffmpeg movie
ffmpeg -framerate 10 -i pe${pe}_ep${ep}_phi${phi}_frame_%04d.png\
 -vcodec libx264 -s 2000x2000 -pix_fmt yuv420p\
 -threads 1 pe${pe}_ep${ep}_phi${phi}.mp4


