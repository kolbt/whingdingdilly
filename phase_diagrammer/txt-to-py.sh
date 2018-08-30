#!/bin/sh
#SBATCH -p general                          # partition to run on
#SBATCH -n 1                                # number of cores
#SBATCH --mem=20000                         # increase memory request (a gig)
#SBATCH -t 1-00:00                          # time (D-HH:MM)

script=$1
txtFiles=("$@")

#python $script/write-phase-txt.py ${txtFiles[@]}
#python $script/renormalize_experiment.py ${txtFiles[@]}
python $script/lennard-jones_diameter_overlay_mono.py ${txtFiles[@]}
