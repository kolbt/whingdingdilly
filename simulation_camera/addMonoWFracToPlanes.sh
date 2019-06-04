#!/bin/sh
camPath="/Users/kolbt/Desktop/compiled/whingdingdilly/simulation_camera"

# Keep xA as a variable (since we can use it now)
# Find where EITHER
# 1. xA = 0 or 100
# 2. peA = peB

# Resolution of an individual simulation frame
size=600
# Start value for PeB plane
pea=0

# Loop through B particle activity
while [ ${pea} -le 150 ]; do

    # Frame to place the monodisperse simulation on
    composite="pb${pea}.png"
    # x-position for this frame
    x=$(bc <<< "scale=0;((($pea/10))*$size)")
    # Get the name of the gsd to use
    for getfile in $(ls pa${pea}_pb0_xa100*gsd); do
        file=$getfile
    done
    # Get rid of .gsd suffix
    base=${file%.gsd}
    # Create the frame (with green particles)
    ovitos ${camPath}/png_final_tstep.py $file 1 $size


    # This places the monodisperse simulation in the vertical column (PeA=PeB)
    xa=10
    while [ ${xa} -le 90 ]; do
        # y-position on composite
        y=$(bc <<< "scale=0;((($xa/10))*$size)")
        # Place image on composite
        python ${camPath}/placeOnComp.py $composite "${base}.png" $x $y
        # Increment xA
        xa=$(( $xa + 10 ))
    done


    # This places the bottom row of frames (xA=0)
    xvar=0
    while [ ${xvar} -le 150 ]; do
        # Get the x position
        x=$(bc <<< "scale=0;((($xvar/10))*$size)")
        # y-position on composite (bottom row)
        y=$(bc <<< "scale=0;(((100/10))*$size)")
        # Place image on composite
        python ${camPath}/placeOnComp.py $composite "${base}.png" $x $y
        # Increment x-position
        xvar=$(( $xvar + 10 ))
    done


    # Remove B-colored image
    rm ${base}.png
    # Create A-colored image
    ovitos ${camPath}/png_final_tstep.py $file 0 $size
    # Place image on each PeB frame in top row
    frames=0
    while [ $frames -le 150 ]; do
        # x-position depends on frame
        x=$(bc <<< "scale=0;((($pea/10))*$size)")
        # y-position is top of composite
        y=0
        # Place (pink) frame on composite
        python ${camPath}/placeOnComp.py "pb${frames}.png" "${base}.png" $x $y
        # Increment frame to place on
        frames=$(( $frames + 10 ))
    done

    # Remove the A-colored image
    rm ${base}.png
    # Increment peB plane counter
    pea=$(( $pea + 10 ))

done

exit 0
