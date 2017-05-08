#!/bin/sh
#SBATCH -p gpu                              # partition to run on
#SBATCH -n 1                                # number of cores
#SBATCH -t 0-04:00                          # time (D-HH:MM)
#SBATCH --mail-type=END,FAIL                # send me email when job finishes or fails
#SBATCH --mail-user=kolbt@live.unc.edu      # address to send to

filename=$1

python $filename                            # run infiles with python
