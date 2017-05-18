#!/bin/sh

current=$( date "+%m_%d_%y" )
this_path=$( pwd )

echo "Are you running on Longleaf (y/n)?"
read answer

if [ $answer == "y" ]; then
    hoomd_path='/nas/longleaf/home/kolbt/programs/hoomd-blue/build'
    gsd_path='/nas/longleaf/home/kolbt/programs/gsd/build'
    script_path='/nas/longleaf/home/kolbt/whingdingdilly'
    template='/nas/longleaf/home/kolbt/whingdingdilly/template.py'
    sedtype='sed'
    submit='sbatch'
else
    hoomd_path='/Users/kolbt/Desktop/compiled/hoomd-blue/build'
    gsd_path='/Users/kolbt/Desktop/compiled/gsd/build'
    script_path='/Users/kolbt/Desktop/whingdingdilly'
    template='/Users/kolbt/Desktop/whingdingdilly/template.py'
    sedtype='gsed'
    submit='sh'
fi

# This is my comment
#This file is explicitly designed to:
#    1.) Write infiles which vary
#        a.) Particle fraction of species A
#        b.) Activity of species A
#    2.) Submit these as batch jobs to SLURM manager
#    3.) Collect the output
#        a.) infiles
#        b.) .gsd files
#        c.) analysis data (.pngs?, .txt?)

mkdir ${current}_parent
cd ${current}_parent

#echo "How many timesteps are in these simulations?"
#read tsteps
#echo "How often is a dumpfile written?"
#read dump_freq
#echo "What spacing should be used for the particle fraction? (as percent, suggested 10)"
#read x_a_spacer
#echo "What spacing should be used for the activity?"
#read pe_a_spacer
#echo "What is the activity of species B?"
#read pe_b

# these are good run settings
tsteps=$(( 10000000 ))
dump_freq=$(( 20000 ))
x_a_spacer=$(( 10 ))
pe_a_spacer=$(( 10 ))
pe_b=$(( 150 ))

# these are good debug settings
#tsteps=$(( 50000 ))
#dump_freq=$(( 10 ))
#x_a_spacer=$(( 50 ))
#pe_a_spacer=$(( 50 ))
#pe_b=$(( 150 ))

x_count=$(( 0 ))    # monodisperse species b
x_max=$(( 100 ))    # monodisperse species a
pe_max=$(( 150 ))   # maximum pe value

# this segment of code writes the infiles
while [ $x_count -le $x_max ]       # loop through particle fraction
do

    pe_count=$(( 0 ))               # start value for each set of fixed x_a simulations

    while [ $pe_count -le $pe_max ] # loop through activity at constant particle fraction
    do

        infile=pa${pe_count}_pb${pe_b}_xa${x_count}.py                          # set unique infile name
        $sedtype -e 's@\${hoomd_path}@'"${hoomd_path}"'@g' $template > $infile  # write path to infile (delimit with @)
        $sedtype -i 's/\${tsteps}/'"${tsteps}"'/g' $infile                      # write tsteps to infile
        $sedtype -i 's/\${dump_freq}/'"${dump_freq}"'/g' $infile                # write dump frequency to infile
        $sedtype -i 's/\${part_frac_a}/'"${x_count}"'/g' $infile                # write particle fraction to infile
        $sedtype -i 's/\${pe_a}/'"${pe_count}"'/g' $infile                      # write activity of A to infile
        $sedtype -i 's/\${pe_b}/'"${pe_b}"'/g' $infile                          # write activity of B to infile
        $sedtype -i 's@\${gsd_path}@'"${gsd_path}"'@g' $infile                  # set gsd path variable

        $submit $script_path/run.sh $infile

        pe_count=$(( $pe_count + $pe_a_spacer ))

    done

    x_count=$(( $x_count + $x_a_spacer ))

done

#sh $script_path/sort_files.sh                       # sort the infiles into batches
#sh $script_path/run_all.sh $script_path $answer     # run each batch

