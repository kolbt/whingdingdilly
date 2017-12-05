#!/bin/sh
#SBATCH -p general                          # partition to run on
#SBATCH -n 1                                # number of cores
#SBATCH --mem=20000                         # increase memory request (a gig)
#SBATCH -t 2-00:00                          # time (D-HH:MM)

pa=$1
pb=$2
xa=$3
hoomd_path=$4
gsd_path=$5
script_path=$6

#python $script_path/nearest_neigh_small_array.py $pa $pb $xa $hoomd_path $gsd_path
#python $script_path/nearest_neigh.py $pa $pb $xa $hoomd_path $gsd_path
#python $script_path/heatmap.py $pa $pb $xa $hoomd_path $gsd_path
#python $script_path/post_proc.py $pa $pb $xa $hoomd_path $gsd_path
#python $script_path/pp_msd_perc_A.py $pa $pb $xa $hoomd_path $gsd_path
#python $script_path/pp_msdten_perc_A.py $pa $pb $xa $hoomd_path $gsd_path
#python $script_path/MCS.py $pa $pb $xa $hoomd_path $gsd_path
#python $script_path/MCSten.py $pa $pb $xa $hoomd_path $gsd_path
#python $script_path/voronoi.py $pa $pb $xa $hoomd_path $gsd_path
python $script_path/active_pressure.py $pa $pb $xa $hoomd_path $gsd_path

# Orientation specific scripts
#myfile=$(pwd)
#mkdir "pa${pa}_pb${pb}_xa${xa}_images"
#cd "pa${pa}_pb${pb}_xa${xa}_images"
#python $script_path/orientations.py $pa $pb $xa $hoomd_path $gsd_path $myfile



