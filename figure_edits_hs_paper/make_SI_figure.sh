#!/bin/bash

run='/Users/kolbt/Desktop/compiled/whingdingdilly/figure_edits_hs_paper/make_SI_figure.py'
pb=0
pbmax=150
width=640
height=480

# Make the file that I'll paste to
python '/Users/kolbt/Desktop/compiled/whingdingdilly/figure_edits_hs_paper/make_composite.py' $width $height

while [ $pb -le $pbmax ];
do

    # Paste both the MIPS and raw images
    python $run $pb $width $height

    # Increment the counter
    pb=$(( $pb + 10 ))

done
