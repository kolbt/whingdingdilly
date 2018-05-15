#!/bin/sh
#SBATCH -p general                          # partition to run on
#SBATCH -n 1                                # number of cores
#SBATCH -t 11-00:00                         # time (D-HH:MM)
#SBATCH --mem=64G                           # memory requested


filename=$1

python $filename                            # run infiles with python
