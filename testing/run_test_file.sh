#!/bin/sh
#SBATCH --qos gpu_access                    # quality of service for gpu job
#SBATCH -p gpu                              # partition to run on
#SBATCH -n 1                                # number of cores
#SBATCH -t 0-10:00                          # time (D-HH:MM)
#SBATCH --mail-type=END,FAIL                # send me email when job finishes or fails
#SBATCH --mail-user=kolbt@live.unc.edu      #address to send to

python "longleaf_test.py"
