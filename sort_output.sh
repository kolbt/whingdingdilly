#!/bin/sh

pa=0
pb=50

while [ $pa -le 150 ]
do

    mkdir pe$pa-$pb

    for filename in $( ls *pa${pa}_pb${pb}* )
    do

        mv $filename "pe$pa-$pb"

    done

    pa=$(( $pa + 10 ))

done
