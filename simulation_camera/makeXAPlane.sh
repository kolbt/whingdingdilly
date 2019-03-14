#!/bin/sh
camPath="/Users/kolbt/Desktop/compiled/whingdingdilly/simulation_camera"

# It's important to note that the upper LEFT corner is 0, 0
# in the python PILlow pixel coordinate system (x, y)

# Set xA
xa=0
# Resolution of an individual simulation frame
size=500
# This is to adjust how y coordinates are implemented
l_box=$(( $size * 15 ))

# Loop through particle fraction (constant on each plane)
while [ ${xa} -le 50 ]; do

    # Create the blank xA plane frame (16 x 16 image)
    python ${camPath}/makeWall.py "xa${xa}.png" $size
    # Set PeB
    peb=0

    # Loop through B particle activity
    while [ ${peb} -le 150 ]; do

        # Set PeA
        pea=0

        # Loop through A particle activity
        while [ ${pea} -le 150 ]; do

            # The coordinates are predefined (color may flip)
            x=$(bc <<< "scale=0;(($pea/10)*$size)")
            y=$(bc <<< "scale=0;($l_box-(($peb/10)*$size))")
            flip=0
            base="pa${pea}_pb${peb}_xa${xa}"

            # Check if file exists (block executes if file does NOT exist)
            if [ ! -f "${base}.gsd" ]; then
                # Flip the color scheme
                flip=1
                # Get the new xa value
                var=$(( 100 - $xa ))
                # Search for the inverse file
                base="pa${peb}_pb${pea}_xa${var}"
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
