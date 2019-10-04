#!/bin/sh
#SBATCH -p general                          # partition to run on
#SBATCH -n 1                                # number of cores
#SBATCH -t 3-00:00                          # time (D-HH:MM)
#SBATCH -o movie.out

# Command to increase memory allocated --mem=100g

ffmpeg -framerate 10 -i wall_%04d.png\
 -vcodec libx264 -s 6000x4500 -pix_fmt yuv420p\
 -threads 1 wallOutLarge.mp4

ffmpeg -framerate 10 -i wall_%04d.png\
 -vcodec libx264 -s 2000x1500 -pix_fmt yuv420p\
 -threads 1 wallOut.mp4

ffmpeg -framerate 10 -i wall_%04d.png\
 -vcodec libx264 -s 1600x1200 -pix_fmt yuv420p\
 -threads 1 wallOutmedium.mp4

ffmpeg -framerate 10 -i wall_%04d.png\
 -vcodec libx264 -s 1000x750 -pix_fmt yuv420p\
 -threads 1 wallOutsmall.mp4
