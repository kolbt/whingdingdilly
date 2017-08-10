#!/bin/sh

hoomd_path='/Users/kolbt/Desktop/compiled/hoomd-blue/build'
gsd_path='/Users/kolbt/Desktop/compiled/gsd/build'
script_path='/Users/kolbt/Desktop/compiled/whingdingdilly/optimization'
sedtype='gsed'

size_min=1500

for filename in $( ls pa*.gsd )
do

    # pull parameters from filename
    pa=$(echo $filename | $sedtype 's/^.*pa\([0-9]*\)_.*/\1/')
    pb=$(echo $filename | $sedtype 's/^.*pb\([0-9]*\)_.*/\1/')
    xa=$(echo $filename | $sedtype 's/^.*xa\([0-9]*\)..*/\1/')

    while [ $size_min -le 5000 ]
    do

        echo "Running minimum cluster size of ${size_min}... "
        # Run python script to generate .png's
        python ${script_path}/optimize_size_min.py\
        $pa $pb $xa $hoomd_path $gsd_path $size_min
        # Run video creation command
        ffmpeg -y -framerate 20 -i sm_${size_min}_opt_%03d.png -vcodec libx264\
        -pix_fmt yuv420p -vf "crop=in_w-180:in_h" sm_${size_min}_opt.mov
        # Increment size_min
        size_min=$(( $size_min + 500 ))

    done

    size_min=$(( 1500 ))

done
