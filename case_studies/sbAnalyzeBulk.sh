#!/bin/sh
#SBATCH -p general                          # partition to run on
#SBATCH -n 1                                # number of cores
#SBATCH -t 3-00:00                          # time (D-HH:MM)
#SBATCH --mem=10g

# Command to increase memory allocated --mem=100g

pe=$1
lat=$2
script_path=$3

python $script_path/analyzeBulkVelocity.py ${pe} ${lat}
