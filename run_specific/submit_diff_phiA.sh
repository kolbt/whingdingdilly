#!/bin/sh

current=$( date "+%m_%d_%y" )
this_path=$( pwd )

echo "Are you running on Longleaf (y/n)?"
read answer

if [ $answer == "y" ]; then
    hoomd_path='/nas/longleaf/home/kolbt/programs/cpu-hoomd/hoomd-blue/build'
    gsd_path='/nas/longleaf/home/kolbt/programs/gsd/build'
    script_path='/nas/longleaf/home/kolbt/whingdingdilly/run.sh'
    template='/nas/longleaf/home/kolbt/whingdingdilly/run_specific/logAndNormalHS.py'
    sedtype='sed'
    submit='sbatch'
else
    hoomd_path='/Users/kolbt/Desktop/compiled/hoomd-blue_11_8_17/hoomd-blue/build'
    gsd_path='/Users/kolbt/Desktop/compiled/gsd/build'
    script_path='/Users/kolbt/Desktop/compiled/whingdingdilly/run.sh'
    template='/Users/kolbt/Desktop/compiled/whingdingdilly/run_specific/logAndNormalHS.py'
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
part_num=$(( 50000 ))
runfor=$(( 100 ))
dump_freq=$(( 20000 ))
phi_spacer=$(( 6 ))
x_a_spacer=$(( 10 ))
pe_a_spacer=$(( 10 ))
pe_b=$(( 0 ))

phi_count=$(( 6 ))
phi_max=$(( 54 ))
x_count=$(( 100 ))  # start at lowest desired fraction of a
x_max=$(( 100 ))    # end at highest desired fraction of a
pe_start=$(( 500 )) # minimum pe value
pe_max=$(( 500 ))   # maximum pe value


# This script is for the specific types of jobs (one row or column of a plane)
loop="n"
while [ $loop == "n" ]
do

    echo "part_num is ${part_num}, input new value or same value."
    read input
    if ! [[ -z "$input" ]]; then
        part_num=$input
    fi

    echo "density, as a percent, is ${phi}."
    read input
    if ! [[ -z "$input" ]]; then
        phi=$input
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

    echo "phi_start is ${phi_count}, input new value or same value."
    read input
    if ! [[ -z "$input" ]]; then
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

    echo "x_start is ${x_count}, input new value or same value."
    read input
    if ! [[ -z "$input" ]]; then
        x_count=$input
    fi

    echo "x_max is ${x_max}, input new value or same value."
    read input
    if ! [[ -z "$input" ]]; then
        x_max=$input
    fi

    echo "x_a_spacer is ${x_a_spacer}, input new value or same value."
    read input
    if ! [[ -z "$input" ]]; then
        x_a_spacer=$input
    fi

    echo "pe_start is ${pe_start}, input new value or same value."
    read input
    if ! [[ -z "$input" ]]; then
        pe_start=$input
    fi

    echo "pe_max is ${pe_max}, input new value or same value."
    read input
    if ! [[ -z "$input" ]]; then
        pe_max=$input
    fi
    echo "pe_a_spacer is ${pe_a_spacer}, input new value or same value."
    read input
    if ! [[ -z "$input" ]]; then
        pe_a_spacer=$input
    fi

    echo "PeB is ${pe_b}, input new value or same value."
    read input
    if ! [[ -z "$input" ]]; then
        pe_b=$input
    fi

    # this shows how long the simulation will run
    tsteps=$(bc <<< "scale=2;$runfor/0.00001")

    echo "part_num is ${part_num}"
    echo "runfor is ${runfor}"
    echo "dump_freq is ${dump_freq}"
    echo "phi_start is ${phi_count}"
    echo "phi_max is ${phi_max}"
    echo "phi_spacer is ${phi_spacer}"
    echo "x_start is ${x_count}"
    echo "x_max is ${x_max}"
    echo "xa_spacer is ${x_a_spacer}"
    echo "pe_start is ${pe_start}"
    echo "pe_max is ${pe_max}"
    echo "pe_a_spacer is ${pe_a_spacer}"
    echo "PeB is ${pe_b}"
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
echo "B activity seed"
read seed5


mkdir ${current}_parent
cd ${current}_parent


# this segment of code writes the infiles
while [ $x_count -le $x_max ]       # loop through particle fraction
do

    pe_count=$(( $pe_start ))    # start value for each set of fixed x_a simulations

    while [ $pe_count -le $pe_max ] # loop through activity at constant particle fraction
    do

        phi_count=$(( $phi_start )) # start value for each set of fixed pe_a and x_a simulations

        while [ $phi_count -le $phi_max ] # loop through activity at constant particle fraction
        do

#            infile=pa${pe_count}_pb${pe_b}_xa${x_count}.py                          # set unique infile name
#            $sedtype -e 's@\${hoomd_path}@'"${hoomd_path}"'@g' $template > $infile  # write path to infile (delimit with @)
#            $sedtype -i 's/\${part_num}/'"${part_num}"'/g' $infile                  # write particle number
#            $sedtype -i 's/\${phi}/'"${phi}"'/g' $infile                            # write particle number
#            $sedtype -i 's/\${runfor}/'"${runfor}"'/g' $infile                      # write time in tau to infile
#            $sedtype -i 's/\${dump_freq}/'"${dump_freq}"'/g' $infile                # write dump frequency to infile
#            $sedtype -i 's/\${part_frac_a}/'"${x_count}"'/g' $infile                # write particle fraction to infile
#            $sedtype -i 's/\${pe_a}/'"${pe_count}"'/g' $infile                      # write activity of A to infile
#            $sedtype -i 's/\${pe_b}/'"${pe_b}"'/g' $infile                          # write activity of B to infile
#            $sedtype -i 's@\${gsd_path}@'"${gsd_path}"'@g' $infile                  # set gsd path variable
#            $sedtype -i 's/\${seed1}/'"${seed1}"'/g' $infile                        # set your seeds
#            $sedtype -i 's/\${seed2}/'"${seed2}"'/g' $infile
#            $sedtype -i 's/\${seed3}/'"${seed3}"'/g' $infile
#            $sedtype -i 's/\${seed4}/'"${seed4}"'/g' $infile
#            $sedtype -i 's/\${seed5}/'"${seed5}"'/g' $infile
#
#            $submit $script_path $infile

            echo $phi_count

            phi_count=$(( $phi_count + $phi_spacer ))

        done

        pe_count=$(( $pe_count + $pe_a_spacer ))

    done

    x_count=$(( $x_count + $x_a_spacer ))

done

