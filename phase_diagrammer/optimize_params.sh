#!/bin/sh

# Run from Hagrid
cd /Volumes/Hagrid
# This script gives optimal agreement for,
#       kappa = 2.220
#       gf    = 0.40

part_num=15000
kappa=2000
gf=30
pb=0

while [ $pb -le 150 ]
do

    # get to the files
    cd pe${pb}B
    cd *

    echo "PeB is: ${pb}"

    while [ $kappa -le 2500 ]
    do
        echo "Kappa is: ${kappa}"

        while [ $gf -le 60 ]
        do

            size_min=$(echo "${part_num}*${gf}/100" |bc)
            occurrences=0
            occur_min=200           # maybe like 100 tsteps?
            output=0

            pa=0
#            pb=100
            xa=0

            if [ -e pe${pb}B_sep.txt ]; then
                rm pe${pb}B_sep.txt             # Delete file if it exists
            fi
            touch pe${pb}B_sep.txt

            echo "Gas Fraction is: ${gf}"
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

            python /Users/kolbt/Desktop/compiled/whingdingdilly/phase_diagrammer/optimize_gasfrac.py\
            pe${pb}B_sep.txt ${pb} ${gf} ${kappa}

            pa=$(( 0 ))
            gf=$(( ${gf} + 5 ))
        done

        gf=$(( 30 ))

        kappa=$(( ${kappa} + 100 ))
    done

    kappa=$(( 2000 ))

    pb=$(( ${pb} + 10 ))
    cd ../..
done
