#!/bin/sh

current=$( date "+%m_%d_%y" )
this_path=$( pwd )

echo "Are you running on Longleaf (y/n)?"
read answer

if [ $answer == "y" ]; then
    hoomd_path='/nas/longleaf/home/kolbt/programs/cpu-hoomd/hoomd-blue/build'
    gsd_path='/nas/longleaf/home/kolbt/programs/gsd/build'
    script_path='/nas/longleaf/home/kolbt/whingdingdilly/run.sh'
    template='/nas/longleaf/home/kolbt/whingdingdilly/run_specific/vary_v_and_eps.py'
    sedtype='sed'
    submit='sbatch'
else
    hoomd_path='/Users/kolbt/Desktop/compiled/hoomd-blue_11_8_17/hoomd-blue/build'
    gsd_path='/Users/kolbt/Desktop/compiled/gsd/build'
    script_path='/Users/kolbt/Desktop/compiled/whingdingdilly/run.sh'
    template='/Users/kolbt/Desktop/compiled/whingdingdilly/run_specific/vary_v_and_eps.py'
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
part_num=$(( 30000 ))
phi=$(( 60 ))
runfor=$(( 100 ))
dump_freq=$(( 20000 ))
xA_spacer=$(( 10 ))
peA_spacer=$(( 10 ))
peB_spacer=$(( 10 ))

xA_start=$(( 10 ))  # exclude monodisperse B
xA_max=$(( 90 ))    # and monodisperse A
peA_start=$(( 0 ))  # minimum pe value
peA_max=$(( 150 ))  # maximum pe value
peB_start=$(( 0 ))  # minimum peB value
peB_max=$(( 150 ))  # maximum peB value

# Ask for input values
#loop="n"
#while [ $loop == "n" ]
#do
#
#    echo "part_num is ${part_num}, input new value or same value."
#    read part_num
#    echo "density, as a percent, is ${phi}."
#    read phi
#    echo "runfor is ${runfor} tau, input new value or same value."
#    read runfor
#    echo "dump_freq is ${dump_freq}, input new value or same value."
#    read dump_freq
#    echo "xA_start is ${xA_start}, input new value or same value."
#    read xA_start
#    echo "xA_max is ${xA_max}, input new value or same value."
#    read xA_max
#    echo "xA_spacer is ${xA_spacer}, input new value or same value."
#    read xA_spacer
#    echo "peA_start is ${peA_start}, input new value or same value."
#    read peA_start
#    echo "peA_max is ${peA_max}, input new value or same value."
#    read peA_max
#    echo "peA_spacer is ${peA_spacer}, input new value or same value."
#    read peA_spacer
#    echo "peB_start is ${peB_start}, input new value or same value."
#    read peB_start
#    echo "peB_max is ${peB_max}, input new value or same value."
#    read peB_max
#    echo "peB_spacer is ${peB_spacer}, input new value or same value."
#    read peB_spacer
#
#    # this shows how long the simulation will run
#    tsteps=$(bc <<< "scale=2;$runfor/0.000001")
#
#    echo "part_num is ${part_num}"
#    echo "runfor is ${runfor}"
#    echo "dump_freq is ${dump_freq}"
#    echo "xA_start is ${xA_start}"
#    echo "xA_max is ${xA_max}"
#    echo "xA_spacer is ${xA_spacer}"
#    echo "peA_start is ${peA_start}"
#    echo "peA_max is ${peA_max}"
#    echo "peA_spacer is ${peA_spacer}"
#    echo "peB_start is ${peB_start}"
#    echo "peB_max is ${peB_max}"
#    echo "peB_spacer is ${peB_spacer}"
#    echo "Simulation will run for ${tsteps} timesteps"
#
#    echo "Are these values okay (y/n)?"
#    read loop
#
#done

# Ask for seeds
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

# Submission routine (no duplicates)
peB_current=$(( $peB_start ))
# While loop for faster activity species
while [ $peB_current -le $peB_max ]
do
    # Reset activity of species A
    peA_current=$(( $peA_start ))
    # While loop for slower activity species (ends at current faster particle activity)
    while [ $peA_current -lt $peB_current ]
    do
        # Reset particle fraction
        xA_current=$(( $xA_start ))
        # While loop for particle fraction (not monodisperse)
        while [ $xA_current -le $xA_max ]
        do

            echo "PeB = ${peB_current}, PeA = ${peA_current}, xA = ${xA_current}"
            # Submit simulations with given parameters
            infile=pa${peA_current}_pb${peB_current}_xa${xA_current}.py                          # set unique infile name
            $sedtype -e 's@\${hoomd_path}@'"${hoomd_path}"'@g' $template > $infile  # write path to infile (delimit with @)
            $sedtype -i 's/\${part_num}/'"${part_num}"'/g' $infile                  # write particle number
            $sedtype -i 's/\${phi}/'"${phi}"'/g' $infile                            # write particle number
            $sedtype -i 's/\${runfor}/'"${runfor}"'/g' $infile                      # write time in tau to infile
            $sedtype -i 's/\${dump_freq}/'"${dump_freq}"'/g' $infile                # write dump frequency to infile
            $sedtype -i 's/\${part_frac_a}/'"${xA_current}"'/g' $infile             # write particle fraction to infile
            $sedtype -i 's/\${pe_a}/'"${peA_current}"'/g' $infile                     # write activity of A to infile
            $sedtype -i 's/\${pe_b}/'"${peB_current}"'/g' $infile                     # write activity of B to infile
            $sedtype -i 's@\${gsd_path}@'"${gsd_path}"'@g' $infile                  # set gsd path variable
            $sedtype -i 's/\${seed1}/'"${seed1}"'/g' $infile                        # set your seeds
            $sedtype -i 's/\${seed2}/'"${seed2}"'/g' $infile
            $sedtype -i 's/\${seed3}/'"${seed3}"'/g' $infile
            $sedtype -i 's/\${seed4}/'"${seed4}"'/g' $infile
            $sedtype -i 's/\${seed5}/'"${seed5}"'/g' $infile
            $submit $script_path $infile

            # Increment particle fraction
            xA_current=$(( $xA_current + $xA_spacer ))

        done
        # Increment activity of species A
        peA_current=$(( peA_current + $peA_spacer ))

    done
    # Increment activity of species B
    peB_current=$(( $peB_current + $peB_spacer ))

done
