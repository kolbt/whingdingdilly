#!/bin/sh

current=$( date "+%m_%d_%y" )
this_path=$( pwd )

echo "Are you running on Longleaf (y/n)?"
read answer

if [ $answer == "y" ]; then
    hoomd_path='/nas/longleaf/home/kolbt/programs/cpu-hoomd/hoomd-blue/build'
    gsd_path='/nas/longleaf/home/kolbt/programs/gsd/build'
    script_path='/nas/longleaf/home/kolbt/whingdingdilly/run.sh'
    tempOne='/nas/longleaf/home/kolbt/whingdingdilly/run_specific/soft_particles.py'
    tempTwo='/nas/longleaf/home/kolbt/whingdingdilly/run_specific/soft_clusters.py'
    sedtype='sed'
    submit='sbatch'
else
    hoomd_path='/Users/kolbt/Desktop/compiled/hoomd-blue/build'
    gsd_path='/Users/kolbt/Desktop/compiled/gsd/build'
    script_path='/Users/kolbt/Desktop/compiled/whingdingdilly/run.sh'
    tempOne='/Users/kolbt/Desktop/compiled/whingdingdilly/run_specific/soft_particles.py'
    tempTwo='/Users/kolbt/Desktop/compiled/whingdingdilly/run_specific/soft_clusters.py'
    sedtype='gsed'
    submit='sh'
fi

echo "GPU (y/n)?"
read gpu

if [ $gpu == "y" ]; then
    hoomd_path='/nas/longleaf/home/kolbt/programs/hoomd_2.2.1/hoomd-blue/build'
    script_path='/nas/longleaf/home/kolbt/whingdingdilly/run_gpu.sh'
fi

# Default values for simulations
part_num=$(( 100000 ))
runfor=$(( 100 ))
dump_freq=$(( 6 ))
pe_start=$(( 100 ))
pe_count=$(( 100 ))
pe_spacer=$(( 100 ))
pe_max=$(( 500 ))
phi_start=$(( 45 ))
phi_count=$(( 45 ))
phi_spacer=$(( 10 ))
phi_max=$(( 65 ))
eps=0.001

# This script is for the specific types of jobs (one row or column of a plane)
loop="n"
while [ $loop == "n" ]
do

    echo "part_num is ${part_num}, input new value or same value."
    read input
    if ! [[ -z "$input" ]]; then
        part_num=$input
    fi

    echo "runfor is ${runfor} tau, input new value or same value."
    read input
    if ! [[ -z "$input" ]]; then
        runfor=$input
    fi

    echo "dump_freq is ${dump_freq}, input new value or same value."
    read input
    if ! [[ -z "$input" ]]; then
        dump_freq=$input
    fi

    echo "pe_start is ${pe_count}, input new value or same value."
    read input
    if ! [[ -z "$input" ]]; then
        pe_start=$input
        pe_count=$input
    fi

    echo "pe_max is ${pe_max}, input new value or same value."
    read input
    if ! [[ -z "$input" ]]; then
        pe_max=$input
    fi

    echo "pe_spacer is ${pe_spacer}, input new value or same value."
    read input
    if ! [[ -z "$input" ]]; then
        pe_spacer=$input
    fi

    echo "phi_start is ${phi_count}, input new value or same value."
    read input
    if ! [[ -z "$input" ]]; then
        phi_start=$input
        phi_count=$input
    fi

    echo "phi_max is ${phi_max}, input new value or same value."
    read input
    if ! [[ -z "$input" ]]; then
        phi_max=$input
    fi

    echo "phi_spacer is ${phi_spacer}, input new value or same value."
    read input
    if ! [[ -z "$input" ]]; then
        phi_spacer=$input
    fi
    
    echo "epsilon is ${eps}, input new value or same value."
    read input
    if ! [[ -z "$input" ]]; then
        eps=$input
    fi

    # this shows how long the simulation will run
    tsteps=$(bc <<< "scale=2;$runfor/0.000001")

    echo "part_num is ${part_num}"
    echo "runfor is ${runfor}"
    echo "dump_freq is ${dump_freq}"
    echo "pe_start is ${pe_count}"
    echo "pe_max is ${pe_max}"
    echo "pe_spacer is ${pe_spacer}"
    echo "phi_start is ${phi_count}"
    echo "phi_max is ${phi_max}"
    echo "phi_spacer is ${phi_spacer}"
    echo "epsilon is ${eps}"
    echo "Simulation will run for ${tsteps} timesteps"

    echo "Are these values okay (y/n)?"
    read loop

done

echo "Time to set some seeds!"
echo "Positional seed"
read seed1
echo "Equilibration seed"
read seed2
echo "Orientational seed"
read seed3
echo "A activity seed"
read seed4


mkdir ${current}_parent
cd ${current}_parent

# Set system density to smallest value
phi_count=$(( $phi_start ))
# Loop through system density
while [ $phi_count -le $phi_max ]
do
    # Reset activity
    pe_count=$(( $pe_start ))
    
    # Loop through activity
    while [ $pe_count -le $pe_max ]
    do

        # Submit randomly initialized simulations
        random=pe${pe_count}_phi${phi_count}_eps${eps}.py                       # set unique infile name
        #'s/\${replace_in_text_File}/'"${variable_to_replace_with}"'/g'
        $sedtype -e 's/\${part_num}/'"${part_num}"'/g' $tempOne > $random       # write particle number
        $sedtype -i 's/\${phi}/'"${phi_count}"'/g' $random                      # write particle number
        $sedtype -i 's/\${runfor}/'"${runfor}"'/g' $random                      # write time in tau to infile
        $sedtype -i 's/\${dump_freq}/'"${dump_freq}"'/g' $random                # write dump frequency to infile
        $sedtype -i 's/\${pe}/'"${pe_count}"'/g' $random                        # write activity to infile
        $sedtype -i 's/\${eps}/'"${eps}"'/g' $random                            # write epsilon to infile
        $sedtype -i 's/\${seed1}/'"${seed1}"'/g' $random                        # set your seeds
        $sedtype -i 's/\${seed2}/'"${seed2}"'/g' $random
        $sedtype -i 's/\${seed3}/'"${seed3}"'/g' $random
        $sedtype -i 's/\${seed4}/'"${seed4}"'/g' $random

        $submit $script_path $random
        
        # Now submit cluster simulation with same parameters
        clust=clust_pe${pe_count}_phi${phi_count}.py                            # set unique infile name
        #'s/\${replace_in_text_File}/'"${variable_to_replace_with}"'/g'
        $sedtype -e 's/\${part_num}/'"${part_num}"'/g' $tempTwo > $clust        # write particle number
        $sedtype -i 's/\${phi}/'"${phi_count}"'/g' $clust                       # write particle number
        $sedtype -i 's/\${runfor}/'"${runfor}"'/g' $clust                       # write time in tau to infile
        $sedtype -i 's/\${dump_freq}/'"${dump_freq}"'/g' $clust                 # write dump frequency to infile
        $sedtype -i 's/\${pe}/'"${pe_count}"'/g' $clust                         # write activity to infile
        $sedtype -i 's/\${eps}/'"${eps}"'/g' $clust                             # write epsilon to infile
        $sedtype -i 's/\${seed1}/'"${seed1}"'/g' $clust                         # set your seeds
        $sedtype -i 's/\${seed2}/'"${seed2}"'/g' $clust
        $sedtype -i 's/\${seed3}/'"${seed3}"'/g' $clust

        $submit $script_path $clust

        # Increment activity
        pe_count=$(( $pe_count + $pe_spacer ))
        
    done
    
    # Increment system density
    phi_count=$(( $phi_count + $phi_spacer ))

done
