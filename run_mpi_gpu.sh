#!/bin/sh
#SBATCH --qos gpu_access                    # quality of service
#SBATCH --gres=gpu:2                        # I want two gpus
#SBATCH --partition=gpu                     # partition to run on
#SBATCH --mem=32G                           # memory requested
#SBATCH --nodes=1                           # run on one node
#SBATCH --ntasks=2                          # task count
#SBATCH --ntasks-per-core=1                 # one task each gpu
#SBATCH --cpus-per-task=1                   # core count
#SBATCH --time=11-00:00                     # time (D-HH:MM)


filename=$1

mpirun -n 2 python $filename --mode=gpu     # I want two gpus
