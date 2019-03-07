#!/bin/sh
camPath="/Users/kolbt/Desktop/compiled/whingdingdilly/simulation_camera"


# Set xA
xa=0
# Resolution of an individual simulation frame
size=600

# Loop through particle fraction (constant on each plane)
while [ ${xa} -le 100 ]; do

    # Set PeB
    peb=0

    # Loop through B particle activity
    while [ ${peb} -le 150 ]; do

        # Create the blank xA plane frame
        python ${camPath}/makeWall.py "xa${xa}.png" $size
        # Set PeA
        pea=0

        # Loop through A particle activity
        while [ ${pea} -le 150 ]; do

            flip=0
            base="pa${pea}_pb${peb}_xa${xa}"
            x=$(bc <<< "scale=0;((($pea/10))*$size)")
            y=$(bc <<< "scale=0;((($peb/10))*$size)")

            # Check if file exists (block executes if file does NOT exist)
            if [ ! -f "${base}.gsd" ]; then
                flip=1
                tmp=$x
                x=$y
                y=$tmp
                base="pa${peb}_pb${pea}_xa${tmp}"
            fi

            # Columns where PeA < PeB (on composite plot) are inverted
            if [ $peb -gt $pea ]; then
                tmp=$(( 100 - $xa ))
            fi

            # If either filename DOES exist, execute this block
            if [ -f "${base}.gsd" ]; then
                # Write the final timestep png
                ovitos ${camPath}/png_final_tstep.py "${base}.gsd" $flip $size
                # Place it on the composite image
                python ${camPath}/placeOnComp.py "xa${xa}.png" "${base}.png" $x $y
                # Remove the source png
                rm "${base}.png"
            fi

            # Increment PeA
            pea=$(( $pea + 10 ))

        done

        # Increment PeA
        peb=$(( $peb + 10 ))

    done

    # Increment xA
    xa=$(( $xa + 10 ))

done

exit 0
