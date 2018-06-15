#!/bin/sh

# Just gonna hard code the values I'm using

current=$( date "+%m_%d_%y" )
this_path=$( pwd )

echo "Are you running on Longleaf (y/n)?"
read answer

if [ $answer == "y" ]; then
    #    hoomd_path='/nas/longleaf/home/kolbt/programs/hoomd-blue/build'
    hoomd_path='/nas/longleaf/home/kolbt/programs/hoomd_2.2.1/hoomd-blue/build'
    gsd_path='/nas/longleaf/home/kolbt/programs/gsd/build'
    script_path='/nas/longleaf/home/kolbt/whingdingdilly/run.sh'
    #    template='/nas/longleaf/home/kolbt/whingdingdilly/run_specific/template_spec.py'
    template='/nas/longleaf/home/kolbt/whingdingdilly/run_specific/changeSoftness.py'
    sedtype='sed'
    submit='sbatch'
else
    #    hoomd_path='/Users/kolbt/Desktop/compiled/hoomd-blue/build'
    hoomd_path='/Users/kolbt/Desktop/compiled/hoomd-blue_11_8_17/hoomd-blue/build'
    gsd_path='/Users/kolbt/Desktop/compiled/gsd/build'
    script_path='/Users/kolbt/Desktop/compiled/whingdingdilly/run.sh'
    #    template='/Users/kolbt/Desktop/compiled/whingdingdilly/run_specific/template_spec.py'
    template='/Users/kolbt/Desktop/compiled/whingdingdilly/run_specific/changeSoftness.py'
    sedtype='gsed'
    submit='sh'
fi

echo "GPU (y/n)?"
read gpu

if [ $gpu == "y" ]; then
    hoomd_path='/nas/longleaf/home/kolbt/programs/hoomd_2.2.1/hoomd-blue/build'
    script_path='/nas/longleaf/home/kolbt/whingdingdilly/run_gpu.sh'
fi

# General param settings
part_num=$(( 20000 ))
phi=$(( 60 ))
runfor=$(( 200 ))
dump_freq=$(( 20000 ))

# Random number seeds
seed1=$(( 476 ))
seed2=$(( 2984 ))
seed3=$(( 227 ))
seed4=$(( 92 ))
seed5=$(( 9784 ))

# Epsilon loop
eps_start=$(( 1 ))
eps_count=$(( $eps_start ))
eps_space=$(( 1 ))
eps_stop=$(( 1 ))

# Fast swimmer activity loop
peF_start=$(( 100 ))
peF_count=$(( $peF_start ))
peF_space=$(( 200 ))
peF_stop=$(( 500 ))

# Slow swimmer activity (as multiple of fast) loop
peS_start=$(( 0 ))
peS_count=$(( $peS_start ))
peS_space=$(( 25 ))
peS_stop=$(( 75 ))

# Slow particle fraction loop
xS_start=$(( 25 ))
xS_count=$(( $xS_start ))
xS_space=$(( 25 ))
xS_stop=$(( 100 ))

# Loop through epsilon
while [ $eps_count -le $eps_stop ]
do

    # Reset the counter
    peF_count=$(( $peF_start ))

    # Loop through PeFast = [100, 300, 500]
    while [ $peF_count -le $peF_stop ]
    do

        # Reset the counter
        peS_count=$(( $peS_start ))

        # Loop through PeSlow = [0, 0.25, 0.5, 0.75] * PeFast
        while [ $peS_count -le $peS_stop ]
        do

            # Compute the slow activity wih basic calculator
            peS=$(bc <<< "scale=0;($peF_count*$peS_count)/100")
            # Reset the counter
            xS_count=$(( $xS_start ))

            # Loop through (slow) particle fraction = [0.25, 0.50, 0.75, 1.0]
            while [ $xS_count -le $xS_stop ]
            do

                # Go ahead and submit your simulations
                infile=pa${peS}_pb${peF_count}_xa${xS_count}_eps${eps_count}.py
                $sedtype -e 's@\${hoomd_path}@'"${hoomd_path}"'@g' $template > $infile
                $sedtype -i 's/\${part_num}/'"${part_num}"'/g' $infile
                $sedtype -i 's/\${phi}/'"${phi}"'/g' $infile
                $sedtype -i 's/\${runfor}/'"${runfor}"'/g' $infile
                $sedtype -i 's/\${dump_freq}/'"${dump_freq}"'/g' $infile
                $sedtype -i 's/\${part_frac_a}/'"${xS_count}"'/g' $infile
                $sedtype -i 's/\${pe_a}/'"${peS}"'/g' $infile
                $sedtype -i 's/\${pe_b}/'"${peF_count}"'/g' $infile
                $sedtype -i 's@\${gsd_path}@'"${gsd_path}"'@g' $infile
                $sedtype -i 's/\${eps}/'"${eps_count}"'/g' $infile
                $sedtype -i 's/\${seed1}/'"${seed1}"'/g' $infile
                $sedtype -i 's/\${seed2}/'"${seed2}"'/g' $infile
                $sedtype -i 's/\${seed3}/'"${seed3}"'/g' $infile
                $sedtype -i 's/\${seed4}/'"${seed4}"'/g' $infile
                $sedtype -i 's/\${seed5}/'"${seed5}"'/g' $infile

                $submit $script_path $infile

#                # Quick way to check that everything is copacetic
#                echo "Epsilon: $eps_count"
#                echo "Fast activity: $peF_count"
#                echo "Slow activity: $peS"
#                echo "Particle Fraction: $xS_count"

                xS_count=$(( $xS_count + $xS_space ))
            done
            peS_count=$(( $peS_count + $peS_space ))
        done
        peF_count=$(( $peF_count + $peF_space ))
    done
    eps_count=$(( $eps_count + $eps_space ))
done










