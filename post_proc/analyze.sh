#!/bin/sh
#SBATCH -p general                          # partition to run on
#SBATCH -n 1                                # number of cores
#SBATCH -t 2-00:00                          # time (D-HH:MM)

pa=$1
pb=$2
xa=$3
hoomd_path=$4
gsd_path=$5
script_path=$6

python $script_path/heatmap.py $pa $pb $xa $hoomd_path $gsd_path
#python $script_path/post_proc.py $pa $pb $xa $hoomd_path $gsd_path
#python $script_path/pp_msd_perc_A.py $pa $pb $xa $hoomd_path $gsd_path


