#!/bin/sh
#SBATCH --qos gpu_access                    # quality of service
#SBATCH --gres=gpu:1                        # I want one gpus
#SBATCH --partition=gpu                     # partition to run on
#SBATCH --time=6-00:00                      # time (D-HH:MM)
#SBATCH --exclude=g0605

filename=$1

python $filename --mode=gpu # I want one gpu
