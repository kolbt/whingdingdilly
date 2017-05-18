#!/bin/sh

pa=0
pb=150

while [ $pa -le 150 ]
do

    cd pe$pa-$pb

    for filename in $( ls pa${pa}_pb${pb}*.gsd )
    do

        ovitos /Users/kolbt/Desktop/snap_final_tstep.py $filename

    done

    pa=$(( $pa + 10 ))
    cd ..

done
