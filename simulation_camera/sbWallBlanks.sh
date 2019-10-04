#!/bin/sh
#SBATCH -p general                          # partition to run on
#SBATCH -n 1                                # number of cores
#SBATCH -t 3-00:00                         # time (D-HH:MM)
#SBATCH -o wallBlank.out

# Command to increase memory allocated --mem=100g

wallFrame=$1
script_path=$2

python ${script_path}/makeWallFrame.py ${wallFrame}
