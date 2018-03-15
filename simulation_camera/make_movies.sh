#!/bin/sh

in_path=$1
script_path="/Users/kolbt/Desktop/compiled/whingdingdilly/simulation_camera"

for sim in $(ls ${in_path}/*.gsd)
do

    ovitos ${script_path}/make_movie.py ${sim}

done
