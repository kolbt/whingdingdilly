#!/bin/sh
#SBATCH -p general                          # partition to run on
#SBATCH -n 1                                # number of cores
#SBATCH -t 0-10:00                          # time (D-HH:MM)


filename=$1

python post_script.py $filename             # run infiles with python
