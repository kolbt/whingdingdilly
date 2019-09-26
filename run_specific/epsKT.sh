#!/bin/sh

current=$( date "+%m_%d_%y" )
this_path=$( pwd )

echo "Are you running on Longleaf (y/n)?"
read answer

if [ $answer == "y" ]; then
    hoomd_path='/nas/longleaf/home/kolbt/programs/cpu-hoomd/hoomd-blue/build'
    gsd_path='/nas/longleaf/home/kolbt/programs/gsd/build'
    script_path='/nas/longleaf/home/kolbt/whingdingdilly/run.sh'
    template='/nas/longleaf/home/kolbt/whingdingdilly/run_specific/epsKT.py'
    sedtype='sed'
    submit='sbatch'
else
    hoomd_path='/Users/kolbt/Desktop/compiled/hoomd-blue/build'
    gsd_path='/Users/kolbt/Desktop/compiled/gsd/build'
    script_path='/Users/kolbt/Desktop/compiled/whingdingdilly/run.sh'
    template='/Users/kolbt/Desktop/compiled/whingdingdilly/run_specific/epsKT.py'
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
dump_freq=$(( 10 ))
pe_count=$(( 0 ))
pe_spacer=$(( 10 ))
pe_max=$(( 500 ))
phi=$(( 60 ))

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

    echo "phi is ${phi}, input new value or same value."
    read input
    if ! [[ -z "$input" ]]; then
        phi=$input
    fi

    # this shows how long the simulation will run
    tsteps=$(bc <<< "scale=2;$runfor/0.00001")

    echo "part_num is ${part_num}"
    echo "runfor is ${runfor}"
    echo "dump_freq is ${dump_freq}"
    echo "pe_start is ${pe_count}"
    echo "pe_max is ${pe_max}"
    echo "pe_spacer is ${pe_spacer}"
    echo "phi is ${phi}"
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


# this segment of code writes the infiles
while [ $pe_count -le $pe_max ]       # loop through particle fraction
do

    infile=pe${pe_count}_phi${phi}.py                                 # set unique infile name
    #'s/\${replace_in_text_File}/'"${variable_to_replace_with}"'/g'
    $sedtype -e 's@\${hoomd_path}@'"${hoomd_path}"'@g' $template > $infile  # write path to infile (delimit with @)
    $sedtype -i 's/\${part_num}/'"${part_num}"'/g' $infile                  # write particle number
    $sedtype -i 's/\${phi}/'"${phi}"'/g' $infile                            # write particle number
    $sedtype -i 's/\${runfor}/'"${runfor}"'/g' $infile                      # write time in tau to infile
    $sedtype -i 's/\${dump_freq}/'"${dump_freq}"'/g' $infile                # write dump frequency to infile
    $sedtype -i 's/\${pe}/'"${pe_count}"'/g' $infile                        # write activity to infile
    $sedtype -i 's@\${gsd_path}@'"${gsd_path}"'@g' $infile                  # set gsd path variable
    $sedtype -i 's/\${seed1}/'"${seed1}"'/g' $infile                        # set your seeds
    $sedtype -i 's/\${seed2}/'"${seed2}"'/g' $infile
    $sedtype -i 's/\${seed3}/'"${seed3}"'/g' $infile
    $sedtype -i 's/\${seed4}/'"${seed4}"'/g' $infile

    $submit $script_path $infile

    pe_count=$(( $pe_count + $pe_spacer ))

done
