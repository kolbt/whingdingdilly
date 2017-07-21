#!/bin/sh
#SBATCH --qos gpu_access                    # quality of service
#SBATCH --gres=gpu:2                        # I want two gpus
#SBATCH -p gpu                              # partition to run on
#SBATCH -N 1                                # keep to one node
#SBATCH -t 11-00:00                         # time (D-HH:MM)


filename=$1

python $filename                            # run infiles with python
