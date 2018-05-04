#!/bin/sh
#SBATCH -p general                          # partition to run on
#SBATCH -n 1                                # number of cores
#SBATCH --mem=20000                         # increase memory request (a gig)
#SBATCH -t 1-00:00                          # time (D-HH:MM)

script_path=$1
txtFiles=("$@")

python $script_path/compareTexts.py ${txtFiles[@]}

