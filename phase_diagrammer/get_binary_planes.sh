#!/bin/sh

pb=0
path=$(pwd)
diag_path="${path}/phase_3D.txt"
echo ${diag_path}

if [ -e "phase_3D.txt" ]; then
    rm "phase_3D.txt"
fi

touch "phase_3D.txt"

n=0

while [ $pb -le 150 ]
do

    cd pe${pb}B
    cd *
    echo $(pwd)

    if [ -e pe${pb}B_sep.txt ]; then
        paste -d' ' "${diag_path}" pe${pb}B_sep.txt > tmpout && mv tmpout "${diag_path}"
    else
        paste -d' ' "phase_3D.txt" zeros.txt > tmpout && mv tmpout "phase_3D.txt"
    fi

    pb=$(( $pb + 10 ))

    cd ../..

done

python /Users/kolbt/Desktop/compiled/whingdingdilly/phase_diagrammer/phase_3D.py phase_3D.txt
