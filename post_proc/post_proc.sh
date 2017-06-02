#!/bin/sh

for filename in $( ls *.gsd )
do

    # pull parameters from filename

    pa=$(echo $filename | gsed 's/^.*pa\([0-9]*\)_.*/\1/')
    pb=$(echo $filename | gsed 's/^.*pb\([0-9]*\)_.*/\1/')
    xa=$(echo $filename | gsed 's/^.*xa\([0-9]*\)..*/\1/')

    echo $pa
    echo $pb
    echo $xa

    sh /Users/kolbt/Desktop/compiled/whingdingdilly/post_proc/analyze.sh $filename $pa $pb $xa

done
