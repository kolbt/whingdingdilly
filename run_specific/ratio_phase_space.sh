#!/bin/sh

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
part_num=$(( 15000 ))
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

# Ask for seeds
#echo "Time to set some seeds!"
#echo "Positional seed"
#read seed1
#echo "Equilibration seed"
#read seed2
#echo "Orientational seed"
#read seed3
#echo "A activity seed"
#read seed4
#echo "B activity seed"
#read seed5

# Submission routine (no duplicates)
peB_current=$(( $peB_start ))
count=$(( 0 ))
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
            # Submit simulations with given parameters
            echo "PeB = ${peB_current}, PeA = ${peA_current}, xA = ${xA_current}"
            # Increment particle fraction
            xA_current=$(( $xA_current + $xA_spacer ))
            count=$(( $count + 1 ))

        done
        # Increment activity of species A
        peA_current=$(( peA_current + $peA_spacer ))

    done
    # Increment activity of species B
    peB_current=$(( $peB_current + $peB_spacer ))

done

echo "Ran ${count} times."
