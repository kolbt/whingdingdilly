#!/bin/sh

camPath='/Users/kolbt/Desktop/compiled/whingdingdilly/simulation_camera'
flip=0
size=1000

for gsd in $(ls *gsd);
do

    ovitos ${camPath}/png_final_tstep.py "${gsd}" $flip $size

done

