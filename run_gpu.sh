#!/bin/sh
#SBATCH --qos gpu_access                    # quality of service
#SBATCH --gres=gpu:1                        # I want one gpus
#SBATCH --partition=gpu                     # partition to run on
#SBATCH --nodes=1                           # run on one node
#SBATCH --time=6-00:00                      # time (D-HH:MM)
#SBATCH --mem=20g
#SBATCH --exclude=g0601

filename=$1

python $filename --mode=gpu # I want one gpu
