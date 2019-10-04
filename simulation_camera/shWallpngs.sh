#!/bin/sh

script_path="/nas/longleaf/home/kolbt/whingdingdilly/simulation_camera"
sedtype='sed'

for sim in $(ls pe*.gsd)
do

    # Make png series for a simulation
    sbatch ${script_path}/sbWallpngs.sh ${sim} ${script_path}
     
done
