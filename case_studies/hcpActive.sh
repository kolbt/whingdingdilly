#!/bin/sh

current=$( date "+%m_%d_%y" )
this_path=$( pwd )

echo "Are you running on Longleaf (y/n)?"
read answer

if [ $answer == "y" ]; then
    script_path='/nas/longleaf/home/kolbt/whingdingdilly/run.sh'
    script_path='/nas/longleaf/home/kolbt/whingdingdilly/run_gpu.sh'
    template='/nas/longleaf/home/kolbt/whingdingdilly/case_studies/hcpActive.py'
    sedtype='sed'
    submit='sbatch'
else
    script_path='/Users/kolbt/Desktop/compiled/whingdingdilly/run.sh'
    template='/Users/kolbt/Desktop/compiled/whingdingdilly/case_studies/hcpActive.py'
    sedtype='gsed'
    submit='sh'
fi

# Default values for simulations
nCol=$(( 50 ))
nRow=$(( 50 ))
# Ballistic particle force
swimStart=$(( 50 ))
swimCount=$(( 50 ))
swimSpace=$(( 10 ))
swimMax=$(( 500 ))
# Lattice spacing
latStart=$(( 70 ))
latCount=$(( 70 ))
latSpace=$(( 5 ))
latMax=$(( 100 ))

# This script is for the specific types of jobs (one row or column of a plane)
loop="n"
while [ $loop == "n" ]
do

    echo "HCP columns is ${nCol}, input new value or same value."
    read input
    if ! [[ -z "$input" ]]; then
        nCol=$input
    fi
    
    echo "HCP rows is ${nRow}, input new value or same value."
    read input
    if ! [[ -z "$input" ]]; then
        nRow=$input
    fi

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
    
    echo "Lattice spacing start is ${latStart}, input new value or same value."
    read input
    if ! [[ -z "$input" ]]; then
        latStart=$input
        latCount=$input
    fi
    
    echo "Lattice spacing spacer is ${latSpace}, input new value or same value."
    read input
    if ! [[ -z "$input" ]]; then
        latSpace=$input
    fi

    echo "Lattice spacing max is ${latMax}, input new value or same value."
    read input
    if ! [[ -z "$input" ]]; then
        latMax=$input
    fi

    echo "HCP columns is ${nCol}"
    echo "HCP rows is ${nRow}"
    echo "Swim force start is ${swimStart}"
    echo "Swim force spacer is ${swimSpace}"
    echo "Swim force max is ${swimMax}"
    echo "Lattice spacing start is ${latStart}"
    echo "Lattice spacing spacer is ${latSpace}"
    echo "Lattice spacing max is ${latMax}"

    echo "Are these values okay (y/n)?"
    read loop

done

mkdir ${current}_parent
cd ${current}_parent

# Loop through the lattice spacing
while [ $latCount -le $latMax ]
do

    # Reset swim force counter
    swimCount=$(( $swimStart ))

    while [ $swimCount -le $swimMax ]
    do

        infile=active_pe${swimCount}_lattice${latCount}.py          # set unique infile name
        #'s/\${replace_in_text_File}/'"${variable_to_replace_with}"'/g'
        $sedtype -e 's/\${nCol}/'"${nCol}"'/g' $template > $infile  # write path to infile (delimit with @)
        $sedtype -i 's/\${nRow}/'"${nRow}"'/g' $infile              # write particle number
        $sedtype -i 's/\${latCount}/'"${latCount}"'/g' $infile      # write particle number
        $sedtype -i 's/\${swimCount}/'"${swimCount}"'/g' $infile    # write time in tau to infile

        $submit $script_path $infile

        swimCount=$(( $swimCount + $swimSpace ))

    done
    
    # Increment latttic spacing counter
    latCount=$(( $latCount + $latSpace ))
    
done
