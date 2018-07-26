#!/bin/sh
#SBATCH -p general                          # partition to run on
#SBATCH -n 1                                # number of cores
#SBATCH --mem=20000                         # increase memory request (a gig)
#SBATCH -t 1-00:00                          # time (D-HH:MM)

script_path=$1
gsd_path=$2
gsdFiles=("$@")

#python $script_path/hist-mode-phiEffective.py ${gsd_path} ${gsdFiles[@]}
python $script_path/position_w_eff_radii.py ${gsd_path} ${gsdFiles[@]}

