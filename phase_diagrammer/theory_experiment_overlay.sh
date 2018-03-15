#!/bin/sh

cd
cd /Volumes/Hagrid/clust_1000

PeB=0
while [ $PeB -le 150 ]
do

    cd pe${PeB}B
    cd *

    python /Users/kolbt/Desktop/compiled/whingdingdilly/phase_diagrammer/phase_diag_plane.py pe${PeB}B_sep.txt ${PeB}
    mv phase_diagram_pb${PeB}.png /Users/kolbt/Desktop/phase_diagrams

    PeB=$(( $PeB + 10 ))
    cd ../..

done

# Rename files for movie
cd /Users/kolbt/Desktop/phase_diagrams
mkdir phase_movie

PeB=0
count=0
while [ $PeB -le 150 ]
do

    cp phase_diagram_pb${PeB}.png phase_movie/phase_diagram_pb${count}.png
    PeB=$(( $PeB + 10 ))
    count=$(( $count + 1 ))

done

cd phase_movie

ffmpeg -framerate 1 -i phase_diagram_pb%d.png\
 -vcodec libx264 -s 1000x1000 -pix_fmt yuv420p -threads 1\
 phase_diagram.mp4
