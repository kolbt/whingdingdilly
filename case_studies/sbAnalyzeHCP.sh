#!/bin/sh
#SBATCH -p general                          # partition to run on
#SBATCH -n 1                                # number of cores
#SBATCH -t 3-00:00                         # time (D-HH:MM)

# Command to increase memory allocated --mem=100g

ball=$1
pe=$2
lat=$3
script_path=$4

python $script_path/analyzeHCP.py ${ball} ${pe} ${lat}
