#!/bin/sh

fixPath='/Users/kolbt/Desktop/compiled/whingdingdilly/simulation_camera/mono_to_mix.py'
parmin=$(( 0 ))
parmax=$(( 100 ))

# Loop through all gsd files in the directory
for gsd in $(ls *gsd);
do

    parfrac=$(( $parmin ))

    while [ $parfrac -le $parmax ];
    do
        # Write final frame with different composition
        python ${fixPath} ${parfrac}
        # Increment counter
        parfrac=$(( $parfrac + 10 ))

    done

done

