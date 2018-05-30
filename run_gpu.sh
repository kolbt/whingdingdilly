#!/bin/sh
#SBATCH --qos gpu_access                    # quality of service
#SBATCH --gres=gpu:1                        # I want one gpus
#SBATCH --partition=gpu                     # partition to run on
#SBATCH --nodes=1                           # run on one node
#SBATCH --nodelist=g0603                    # specific gpu to run on
#SBATCH --time=11-00:00                     # time (D-HH:MM)

filename=$1

python $filename --mode=gpu # I want one gpu
