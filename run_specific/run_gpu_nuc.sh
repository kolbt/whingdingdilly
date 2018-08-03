#!/bin/sh
#SBATCH --qos gpu_access                    # quality of service
#SBATCH --gres=gpu:1                        # I want one gpus
#SBATCH --partition=gpu                     # partition to run on
#SBATCH --nodes=1                           # run on one node
#SBATCH --time=11-00:00                     # time (D-HH:MM)

inFile=$1
hoomdPath=$2
gsdPath=$3
pa=$4
pb=$5
xa=$6
ep=$7
seed1=$8
seed2=$9
seed3=${10}
seed4=${11}
seed5=${12}
myFrame=${13}

python $inFile $hoomdPath $gsdPath $pa $pb $xa $ep $seed1 $seed2 $seed3 $seed4 $seed5 $myFrame --mode=gpu
