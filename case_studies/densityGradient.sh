#!/bin/sh

current=$( date "+%m_%d_%y" )
this_path=$( pwd )

echo "Are you running on Longleaf (y/n)?"
read answer

if [ $answer == "y" ]; then
    script_path='/nas/longleaf/home/kolbt/whingdingdilly/run.sh'
    script_path='/nas/longleaf/home/kolbt/whingdingdilly/run_gpu.sh'
    template='/nas/longleaf/home/kolbt/whingdingdilly/case_studies/densityGradient.py'
    sedtype='sed'
    submit='sbatch'
else
    script_path='/Users/kolbt/Desktop/compiled/whingdingdilly/run.sh'
    template='/Users/kolbt/Desktop/compiled/whingdingdilly/case_studies/densityGradient.py'
    sedtype='gsed'
    submit='sh'
fi

# Ballistic particle force
swimStart=$(( 50 ))
swimCount=$(( 50 ))
swimSpace=$(( 50 ))
swimMax=$(( 500 ))

# This script is for the specific types of jobs (one row or column of a plane)
loop="n"
while [ $loop == "n" ]
do

    echo "Swim force start is ${swimStart}, input new value or same value."
    read input
    if ! [[ -z "$input" ]]; then
        swimStart=$input
        swimCount=$input
    fi
    
    echo "Swim force spacer is ${swimSpace}, input new value or same value."
    read input
    if ! [[ -z "$input" ]]; then
        swimSpace=$input
    fi

    echo "Swim force max is ${swimMax}, input new value or same value."
    read input
    if ! [[ -z "$input" ]]; then
        swimMax=$input
    fi

    echo "Swim force start is ${swimStart}"
    echo "Swim force spacer is ${swimSpace}"
    echo "Swim force max is ${swimMax}"

    echo "Are these values okay (y/n)?"
    read loop

done

mkdir ${current}_parent
cd ${current}_parent

# This is a simple way to get different seeds
seeds=(2937 1949 1929 719 3784)
# Initialize count to 0
count=$(( 0 ))

for i in ${seeds[@]}
do
    # Run at next seed
    inSeed=${seeds[${count}]}
    # Reset swim force counter
    swimCount=$(( $swimStart ))

    while [ $swimCount -le $swimMax ]
    do

        infile=density_gradient_pe${swimCount}_seed${inSeed}.py                   # set unique infile name
        #'s/\${replace_in_text_File}/'"${variable_to_replace_with}"'/g'
        $sedtype -e 's/\${swimCount}/'"${swimCount}"'/g' $template > $infile    # write time in tau to infile
        $sedtype -i 's/\${inSeed}/'"${inSeed}"'/g' $infile

        $submit $script_path $infile

        swimCount=$(( $swimCount + $swimSpace ))

    done
    
    # Increment counter
    count=$(( $count + 1 ))
    
done
