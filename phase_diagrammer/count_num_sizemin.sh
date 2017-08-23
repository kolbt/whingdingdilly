#!/bin/sh

#size_min=5300           # maybe like 2000?
size_min=4500           # this reflects the optimal parameter (gf=0.3)
occurrences=0
occur_min=200           # maybe like 100 tsteps?
output=0

pa=0
pb=0
xa=0

while [ $pb -le 150 ]
do

    cd "pe${pb}B"
    cd *
    echo $(pwd)

    if [ -e pe${pb}B_sep.txt ]; then
        rm pe${pb}B_sep.txt             # Delete file if it exists
    fi
    touch pe${pb}B_sep.txt

    while [ $pa -le 150 ]
    do

        while [ $xa -le 100 ]
        do
            #####################################
            ### Read the largest clust files ####
            for i in $(cat largest_pa${pa}_pb${pb}_xa${xa}.txt | awk '{print $1}')
            do
                if [ $i -ge $size_min ]; then           # are we past the size threshold?
                    occurrences=$(( $occurrences + 1 ))
                fi
            done

            if [ $occurrences -ge $occur_min ]; then    # do we pass that threshold enough?
                output=$(( 1 ))                         # if so, we have phase sep (=1)
            fi

            echo $output >> pe${pb}B_sep.txt
            ###                               ###
            #####################################

            occurrences=$(( 0 ))    # reset count
            output=$(( 0 ))         # reset output value
            xa=$(( $xa + 10 ))      # increment in-counter
        done

        xa=$(( 0 ))                 # reset in-counter
        pa=$(( $pa + 10 ))          # increment out-counter
    done

    python /Users/kolbt/Desktop/compiled/whingdingdilly/phase_diagrammer/plot_binary_phase.py pe${pb}B_sep.txt ${pb}

    pa=$(( 0 ))
    xa=$(( 0 ))
    pb=$(( $pb + 10 ))
    cd ../..

done
