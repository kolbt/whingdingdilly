#!/bin/sh
#SBATCH -p general                          # partition to run on
#SBATCH -n 1                                # number of cores
#SBATCH -t 0-10:00                          # time (D-HH:MM)

pa=$1
pb=$2
xa=$3
hoomd_path=$4
gsd_path=$5
script_path=$6

#python $script_path/post_proc.py $pa $pb $xa $hoomd_path $gsd_path
python $script_path/test_method.py $pa $pb $xa $hoomd_path $gsd_path

