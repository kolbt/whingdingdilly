#!/bin/sh

exec_path="/Users/kolbt/Desktop/compiled/whingdingdilly/simulation_camera"

pb=$1
pa=0
xa=0

if [ -z "$1" ]; then
    echo "You haven't passed a value for PeB, set one now."
    read pb
fi

x_loc=0
y_loc=1800

# write an empty image
python $exec_path/write_new.py peB${pb}.png

theory_counter=4

while [ $pa -le 150 ]
do

    while [ $xa -le 100 ]
    do
        #store name in variable
        filename="pa${pa}_pb${pb}_xa${xa}"
        # check for file existence
        if [ -f "${filename}.gsd" ]; then
            # look at the theory text file so you can color the background appropriately
            theory_binary=$(sed "${theory_counter}q;d" "/Volumes/Hagrid/theory_txts/peB${pb}_theory.mtx")
            # take snapshot
            ovitos $exec_path/snap_final_tstep.py ${filename}.gsd ${theory_binary}
            # crop and paste
            python $exec_path/crop_place.py final_tstep_${filename}.png peB${pb}.png $x_loc $y_loc
        fi
        xa=$(( $xa + 10 ))
        y_loc=$(( $y_loc - 180 ))
        theory_counter=$(( $theory_counter + 1 ))

    done

    pa=$(( $pa + 10 ))
    xa=$(( 0 ))
    x_loc=$(( $x_loc + 180 ))
    y_loc=$(( 1800 ))

done
