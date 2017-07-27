#!/bin/sh

echo "Transfer TO desktop?"
read which_ssh

echo "What is the path to the parent directory?"
read parent_path

if [ $which_sftp -eq "y" ]; then
    sftp_cmd="kolbt@longleaf.unc.edu"
fi

sftp ${sftp_cmd}
cd ${parent_path}

for tarname in $( ls *.gz )
do

    get ${tarname}

done

#declare -a tar_reps=()
#count=0
#total=0
#for tarname in $( ls *.gz )
#do
#
#tar_reps[$count]="${tarname}"
#count=$(( $count + 1 ))
#total=$(( $total + 1 ))
#
#done
#
#count=$(( 0 ))
#while [ $count -lt $total ]
#do
#
#echo "${tar_reps[$count]}"
#count=$(( $count + 1 ))
#
#done
#
#find | grep post_proc.sh
