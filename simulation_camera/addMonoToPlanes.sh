#!/bin/sh
camPath="/Users/kolbt/Desktop/compiled/whingdingdilly/simulation_camera"

# Resolution of an individual simulation frame
size=600
# Start value for PeB plane
pea=0

# Loop through B particle activity
while [ ${pea} -le 20 ]; do

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
    # Create the frame
    ovitos ${camPath}/png_final_tstep.py $file 1 $size

    xa=10
    # This places the monodisperse simulation in the vertical column (PeA=PeB)
    while [ ${xa} -le 90 ]; do
        # y-position on composite
        y=$(bc <<< "scale=0;((($xa/10))*$size)")
        # Place image on composite
        python ${camPath}/placeOnComp.py $composite "${base}.png" $x $y
        # Increment xA
        xa=$(( $xa + 10 ))
    done

    # This places the bottom row of frames (xA=0)
    while [ ${xa} -le 90 ]; do
        x=$(bc <<< "scale=0;((($pea/10))*$size)")
        # y-position on composite
        y=0
        # Place image on composite
        python ${camPath}/placeOnComp.py $composite "${base}.png" $x $y
        # Increment xA
        xa=$(( $xa + 10 ))
    done

    pea=$(( $pea + 10 ))

done

exit 0
