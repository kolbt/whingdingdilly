#!/bin/sh
camPath="/Users/kolbt/Desktop/compiled/whingdingdilly/simulation_camera"

flip=0
size=1000

for gsd in $(ls pa*gsd);
do

    pa=$(echo $gsd | gsed 's/^.*pa\([0-9]*\)_.*/\1/')
    pb=$(echo $gsd | gsed 's/^.*pb\([0-9]*\)_.*/\1/')
    xa=$(echo $gsd | gsed 's/^.*xa\([0-9]*\)_.*/\1/')
    ep=$(echo $gsd | gsed 's/^.*ep\([0-9]*\)..*/\1/')
    base="pa${pa}_pb${pb}_xa${xa}_ep${ep}"
    ovitos ${camPath}/png_final_tstep.py "${base}.gsd" $flip $size

done

exit 0
