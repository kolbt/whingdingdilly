#!/bin/sh
#SBATCH -p general                          # partition to run on
#SBATCH -n 1                                # number of cores
#SBATCH --mem=20000                         # increase memory request (a gig)
#SBATCH -t 1-00:00                          # time (D-HH:MM)

txtFiles=("$@")
script='/nas/longleaf/home/kolbt/whingdingdilly/phase_diagrammer/write-phase-txt.py'

python $script ${txtFiles[@]}
