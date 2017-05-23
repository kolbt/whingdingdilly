#!/bin/sh

exec_path="/Users/kolbt/Desktop/compiled/whingdingdilly/simulation_camera"

# dir names "pe#-#"

pa=0
pb=100
xa=0

x_loc=0
y_loc=18000

# write an empty image
python $exec_path/lg_write_new.py lg_peB${pb}.png

while [ $pa -le 150 ]
do

    mv lg_peB${pb}.png pe${pa}-${pb}
    cd pe${pa}-${pb}

    while [ $xa -le 100 ]
    do
        #store name in variable
        filename="pa${pa}_pb${pb}_xa${xa}"
        # take snapshot
        #ovitos $exec_path/snap_final_tstep.py ${filename}.gsd
        # crop and paste
        python $exec_path/lg_crop_place.py final_tstep_${filename}.png lg_peB${pb}.png $x_loc $y_loc

        xa=$(( $xa + 10 ))
        y_loc=$(( $y_loc - 1800 ))

    done

    mv lg_peB${pb}.png ../
    cd ..

    pa=$(( $pa + 10 ))
    xa=$(( 0 ))
    x_loc=$(( $x_loc + 1800 ))
    y_loc=$(( 18000 ))

done
