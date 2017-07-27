#!/bin/sh

echo "Transfer TO desktop?"
read which_sftp

echo "What is the path to the parent directory?"
read parent_path

if [ $which_sftp == "y" ]; then
    sftp_cmd="kolbt@longleaf.unc.edu"
fi

# This is how you run commands into sftp
# alternatively use -b batch_file.my flag
sftp -R 1 ${sftp_cmd} << EOF
cd ${parent_path}
get *.gz
EOF

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
