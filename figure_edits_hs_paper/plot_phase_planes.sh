#!/bin/bash

run='/Users/kolbt/Desktop/compiled/whingdingdilly/figure_edits_hs_paper/plot_planes_and_MIPS.py'

for planes in $(ls pb*png);
do

    pb=$(echo $planes | gsed 's/^.*pb\([0-9]*\)..*/\1/')
    python $run $pb

done
