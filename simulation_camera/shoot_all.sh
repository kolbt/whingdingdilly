#!/bin/sh

# starting at 10 cause I already did 0
# run from Hagrid
#cd /Volumes/Hagrid/clust_1000

pb=0
while [ $pb -le 150 ]
do

    cd pe${pb}B
    cd *
    sh /Users/kolbt/Desktop/compiled/whingdingdilly/simulation_camera/snap_composite.sh $pb
    cd ../..
    pb=$(( $pb + 10 ))

done
