#!/bin/sh
#SBATCH --qos gpu_access                    # quality of service
#SBATCH --gres=gpu:2                        # I want two gpus
#SBATCH --partition=gpu                     # partition to run on
#SBATCH --mem=8G                            # memory requested
#SBATCH --nodes=1                           # run on one node
#SBATCH --ntasks=1                          # task count
#SBATCH --cpus-per-task=4                   # core count
#SBATCH --time=11-00:00                     # time (D-HH:MM)
#SBATCH --output=slurm-%j.out               # get Slurm output


filename=$1

python $filename                            # run infiles with python
