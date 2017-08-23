#!/bin/sh

# What this will do:
# 1. initialize counter for set gf and kap
# 2. Count number of mismatches over ALL pb planes
# 3. print the total number of mismatches

cd /Volumes/Hagrid/clust_1000
part_num=15000
occur_min=200
all_counter=0
kappa=1000
gf=30
pb=0
pa=0
xa=0

while [ $kappa -le 2000 ]
do

    while [ $gf -le 60 ]
    do

        size_min=$(echo "${part_num}*${gf}/100" |bc)
        occurrences=0
        output=0

        while [ $pb -le  150 ]
        do

            cd pe${pb}B
            cd *

            if [ -e pe${pb}B_sep.txt ]; then
                rm pe${pb}B_sep.txt                     # Delete file if it exists
            fi
            touch pe${pb}B_sep.txt

            while [ $pa -le 150 ]
            do

                while [ $xa -le 100 ]
                do

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

                    occurrences=$(( 0 ))    # reset count
                    output=$(( 0 ))         # reset output value
                    xa=$(( $xa + 10 ))

                done
                xa=$(( 0 ))
                pa=$(( $pa + 10 ))

            done

            python /Users/kolbt/Desktop/compiled/whingdingdilly/optimization/optimize_gasfrac.py\
            pe${pb}B_sep.txt ${pb} ${gf} ${kappa}

            tmp=$?
            all_counter=$(( $all_counter + $tmp ))

            pa=$(( 0 ))
            pb=$(( $pb + 10 ))
            cd ../..

        done
        pb=$(( 0 ))
        # output allcounter for kappa and gf here
        echo "For kappa = ${kappa} and gf = ${gf} the number of mismatches is ${all_counter}"
        gf=$(( $gf + 5 ))
        all_counter=$(( 0 ))

    done
    gf=$(( 30 ))
    kappa=$(( $kappa + 50 ))

done
