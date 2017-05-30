#!/bin/sh


for filename in $( ls *.gsd )
do

    ovitos /Users/kolbt/Desktop/compiled/whingdingdilly/simulation_camera/snap_final_tstep.py $filename

done
