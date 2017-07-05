#!/bin/sh

pb=0

if [ -e "phase_3D.txt" ]; then
    rm "phase_3D.txt"
fi

touch "phase_3D.txt"

n=0

while [ $pb -le 150 ]
do

    if [ -e pe${pb}B_sep.txt ]; then
        paste -d' ' "phase_3D.txt" pe${pb}B_sep.txt > tmpout && mv tmpout "phase_3D.txt"
    else
        paste -d' ' "phase_3D.txt" zeros.txt > tmpout && mv tmpout "phase_3D.txt"
    fi

    pb=$(( $pb + 10 ))

done

python /Users/kolbt/Desktop/compiled/whingdingdilly/phase_diagrammer/phase_3D.py phase_3D.txt
