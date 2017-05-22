#!/bin/sh

exec_path="/Users/kolbt/Desktop/compiled/whingdingdilly/simulation_camera"

pa=0
pb=100
xa=0

x_loc=0
y_loc=1800

# write an empty image
python $exec_path/write_new.py peB${pb}.png

while [ $pa -le 150 ]
do

    while [ $xa -le 100 ]
    do
        #store name in variable
        filename="pa${pa}_pb${pb}_xa${xa}"
        # take snapshot
        ovitos $exec_path/snap_final_tstep.py ${filename}.gsd
        # crop and paste
        python $exec_path/crop_place.py final_tstep_${filename}.png peB${pb}.png $x_loc $y_loc

        xa=$(( $xa + 10 ))
        y_loc=$(( $y_loc - 180 ))

    done

    pa=$(( $pa + 10 ))
    xa=$(( 0 ))
    x_loc=$(( $x_loc + 180 ))
    y_loc=$(( 1800 ))

done
