#!/bin/sh

exec_path="/Users/kolbt/Desktop/compiled/whingdingdilly/simulation_camera"

pa=50
pb=150
xa=60

x_loc=0
y_loc=0

# write an empty image
python $exec_path/lg_write_new.py peB${pb}.png

while [ $pa -le 50 ]
do

    while [ $xa -le 100 ]
    do
        #store name in variable
        filename="pa${pa}_pb${pb}_xa${xa}"
        # take snapshot
        #ovitos $exec_path/snap_final_tstep.py ${filename}.gsd
        # crop and paste
        python $exec_path/lg_crop_place.py final_tstep_${filename}.png peB${pb}.png $x_loc $y_loc

        xa=$(( $xa + 2 ))
        x_loc=$(( $x_loc + 1800 ))

    done

    pa=$(( $pa + 10 ))
    xa=$(( 110 ))                   # run loop once
#    x_loc=$(( $x_loc + 180 ))
#    y_loc=$(( 1800 ))

done
