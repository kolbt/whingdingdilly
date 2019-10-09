#!/bin/sh
#SBATCH -p general                          # partition to run on
#SBATCH -n 1                                # number of cores
#SBATCH -t 11-00:00                         # time (D-HH:MM)

# Command to increase memory allocated --mem=100g

pa=$1
pb=$2
xa=$3
hoomd_path=$4
gsd_path=$5
script_path=$6
ep=$7
phi=$8

#python $script_path/nearest_neigh_small_array.py $pa $pb $xa $hoomd_path $gsd_path
#python $script_path/nearest_neigh.py $pa $pb $xa $hoomd_path $gsd_path
#python $script_path/heatmap.py $pa $pb $xa $hoomd_path $gsd_path
#python $script_path/post_proc.py $pa $pb $xa $hoomd_path $gsd_path
#python $script_path/pp_msd_perc_A.py $pa $pb $xa $hoomd_path $gsd_path
#python $script_path/pp_msdten_perc_A.py $pa $pb $xa $hoomd_path $gsd_path
#python $script_path/MCS.py $pa $pb $xa $hoomd_path $gsd_path
#python $script_path/MCSten.py $pa $pb $xa $hoomd_path $gsd_path
#python $script_path/voronoi.py $pa $pb $xa $hoomd_path $gsd_path
#python $script_path/meshed_output.py $pa $pb $xa $hoomd_path $gsd_path
#python $script_path/per_particle_output.py $pa $pb $xa $hoomd_path $gsd_path
#python $script_path/gtar_pressure.py $pa $pb $xa $hoomd_path $gsd_path
#python $script_path/phase_types.py $pa $pb $xa $hoomd_path $gsd_path
#python $script_path/dense_CoM.py $pa $pb $xa $hoomd_path $gsd_path
#python $script_path/number_densities.py $pa $pb $xa $hoomd_path $gsd_path
#python $script_path/force_diff_sources.py $pa $pb $xa $hoomd_path $gsd_path
#python $script_path/histogram_distance.py $pa $pb $xa $hoomd_path $gsd_path $ep
#python $script_path/histogram_output_txt.py $pa $pb $xa $hoomd_path $gsd_path $ep
#python $script_path/plotNumberDensities.py $pa $pb $xa
#python $script_path/pairCorrelationRelations.py $pa $pb $xa $hoomd_path $gsd_path
#python $script_path/mesh_nearest_neighbor.py $pa $pb $xa $hoomd_path $gsd_path
#python $script_path/computeRDF.py $pa $pb $xa $hoomd_path $gsd_path $ep
#python $script_path/heatmapType.py $pa $pb $xa $hoomd_path $gsd_path $ep
#python $script_path/extrinsic_txt.py $pa $pb $xa $hoomd_path $gsd_path $ep
#python $script_path/extrinsic_all_time.py $pa $pb $xa $hoomd_path $gsd_path $ep
#python $script_path/phi_extrinsic_all_time.py $pa $pb $xa $hoomd_path $gsd_path $ep $phi
#python $script_path/sim_to_txt_no_cutoff.py $pa $pb $xa $hoomd_path $gsd_path $ep
#python $script_path/convergence_analysis.py $pa $pb $xa $hoomd_path $gsd_path $ep $num
#python $script_path/extrinsic_all_time_back_compat.py $pa $pb $xa $hoomd_path $gsd_path $ep
#python $script_path/edge_detection.py $pa $pb $xa $hoomd_path $gsd_path
#python $script_path/edge_detection_v2.py $pa $pb $xa $hoomd_path $gsd_path
#python $script_path/check_cluster_alg.py $pa $pb $xa $hoomd_path $gsd_path $ep
#python $script_path/diffHeatmapType.py $pa $pb $xa $hoomd_path $gsd_path $ep
#python $script_path/orientation_snapshots.py $pa $pb $xa $hoomd_path $gsd_path $ep
#python $script_path/binnedNetActivity.py $pa $pb $xa $hoomd_path $gsd_path $ep
#python $script_path/analyze_alpha.py $pa $pb $xa $hoomd_path $gsd_path $ep $al
#python $script_path/alpha_diameter_histogram.py $pa $pb $xa $hoomd_path $gsd_path $ep $al
#python $script_path/voronoi_tessellation.py $pa $pb $xa $hoomd_path $gsd_path $ep
#python $script_path/mesh_nearest_neighbors.py $pa $pb $xa $hoomd_path $gsd_path $ep
#python $script_path/mesh_nearest_neighbors_periodic.py $pa $pb $xa $hoomd_path $gsd_path $ep
#python $script_path/delta_spatial.py $pa $pb $xa $hoomd_path $gsd_path $ep
#python $script_path/soft_nearest_neighbors.py $pa $pb $xa $hoomd_path $gsd_path $ep
python $script_path/compute_phase_area.py $pa $pb $xa $hoomd_path $gsd_path $ep $phi


## Movie for defects
#ffmpeg -start_number 0 -framerate 6 -i defects_pa${pa}_pb${pb}_xa${xa}_ep${ep}_frame%04d.png\
# -vcodec libx264 -s 1600x1200 -pix_fmt yuv420p -threads 1\
# defects_pa${pa}_pb${pb}_xa${xa}_ep${ep}.mp4

# Movie for RDF
#ffmpeg -framerate 10 -i RDF_pa${pa}_pb${pb}_xa${xa}_ep${ep}_fm%d.png\
# -vcodec libx264 -s 1000x1000 -pix_fmt yuv420p -threads 1\
# RDF_pa${pa}_pb${pb}_xa${xa}_ep${ep}.mp4
#rm RDF_pa${pa}_pb${pb}_xa${xa}_ep${ep}_fm*.png

# Movies for net activity and other force relations
#ffmpeg -start_number 450 -framerate 10 -i forces_pa${pa}_pb${pb}_xa${xa}_eps${ep}_fm%d.png\
# -vcodec libx264 -s 2000x2000 -pix_fmt yuv420p -threads 1\
# forces_pa${pa}_pb${pb}_xa${xa}_ep${ep}.mp4
#rm forces_pa${pa}_pb${pb}_xa${xa}_eps${ep}_fm*.png
#
#ffmpeg -start_number 450 -framerate 10 -i netPe_pa${pa}_pb${pb}_xa${xa}_eps${ep}_fm%d.png\
# -vcodec libx264 -s 2000x2000 -pix_fmt yuv420p -threads 1\
# netPe_pa${pa}_pb${pb}_xa${xa}_ep${ep}.mp4
#rm netPe_pa${pa}_pb${pb}_xa${xa}_eps${ep}_fm*.png

# Movie for heatmap by type
#ffmpeg -start_number 450 -framerate 10 -i heatType_pa${pa}_pb${pb}_xa${xa}_ep${ep}_fm%d.png\
# -vcodec libx264 -s 2000x2000 -pix_fmt yuv420p -threads 1\
# heatType_pa${pa}_pb${pb}_xa${xa}_ep${ep}.mp4
#rm heatType_pa${pa}_pb${pb}_xa${xa}_ep${ep}_fm*.png

# Command to make movie for checking cluster algorithm
#ffmpeg -framerate 10 -i pa${pa}_pb${pb}_xa${xa}_ep${ep}_frame%d.png\
# -vcodec libx264 -s 2000x2000 -pix_fmt yuv420p -threads 1\
# clust_alg_pa${pa}_pb${pb}_xa${xa}_ep${ep}.mp4
#
#rm pa${pa}_pb${pb}_xa${xa}_ep${ep}_frame*.png

# Center of mass movie
#ffmpeg -framerate 10 -i mvy_pa${pa}_pb${pb}_xa${xa}_%d.png\
# -vcodec libx264 -s 1000x1000 -pix_fmt yuv420p -threads 1\
# CoM_pa${pa}_pb${pb}_xa${xa}.mp4

# Orientation specific scripts
#myfile=$(pwd)
#mkdir "pa${pa}_pb${pb}_xa${xa}_images"
#cd "pa${pa}_pb${pb}_xa${xa}_images"

#python $script_path/orientations.py $pa $pb $xa $hoomd_path $gsd_path $myfile $ep
#python $script_path/orientationsCentered.py $pa $pb $xa $hoomd_path $gsd_path $myfile

#ffmpeg -start_number 450 -framerate 10 -i orientation_pa${pa}_pb${pb}_xa${xa}_ep${ep}_fm%d.png\
# -vcodec libx264 -s 2000x2000 -pix_fmt yuv420p -threads 1\
# orientation_pa${pa}_pb${pb}_xa${xa}_ep${ep}.mp4
#rm orientation_pa${pa}_pb${pb}_xa${xa}_ep${ep}_fm*.png

# Move the movie once it's been made
#mv orientation_pa${pa}_pb${pb}_xa${xa}.mp4 ../orientation*

#ffmpeg -framerate 10 -i tot_press_pa${pa}_pb${pb}_xa${xa}_mvout_%d.png\
# -vcodec libx264 -s 1000x1000 -pix_fmt yuv420p -threads 1\
# pressure_pa${pa}_pb${pb}_xa${xa}.mp4

# Movies for binned vector force
#ffmpeg -framerate 10 -i nBins100_pa${pa}_pb${pb}_xa${xa}_step_%d.png\
# -vcodec libx264 -s 1000x1000 -pix_fmt yuv420p -threads 1\
# nBins100_pa${pa}_pb${pb}_xa${xa}.mp4
#
#rm nBins100_pa${pa}_pb${pb}_xa${xa}*.png

