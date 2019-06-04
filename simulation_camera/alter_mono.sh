#!/bin/sh

# You can only run this file locally
hoomd_path='/Users/kolbt/Desktop/compiled/hoomd-blue_11_8_17/hoomd-blue/build'
gsd_path='/Users/kolbt/Desktop/compiled/gsd/build'
fixPath='/Users/kolbt/Desktop/compiled/whingdingdilly/simulation_camera/mono_to_mix.py'
parmin=$(( 0 ))
parmax=$(( 100 ))

# Loop through all gsd files in the directory
for gsd in $(ls pa*gsd);
do

    pa=$(echo $gsd | gsed 's/^.*pa\([0-9]*\)_.*/\1/')
    parfrac=$(( $parmin ))

    while [ $parfrac -le $parmax ];
    do
        # Write final frame with different composition
        python ${fixPath} ${gsd} ${pa} ${parfrac} ${hoomd_path} ${gsd_path}
        # Increment counter
        parfrac=$(( $parfrac + 10 ))

    done

done
