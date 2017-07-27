#!/bin/sh

echo "Transfer TO desktop?"
read which_sftp

if [ $which_sftp == "y" ]; then
sftp_cmd="kolbt@longleaf.unc.edu"
#else
#    sftp_cmd="not_a_private_repo..."
fi

echo "Do you know the filename?"
read just_find

if [ $just_find == "y" ]; then
    echo "What is it?"
    read fname
    size=$(printf "%s" "$fname" | wc -m)
else
    echo "Okay, what is the path to the parent directory (on remote)?"
    read parent_path
fi

ssh ${sftp_cmd} << EOF
tmp_path=$(find | grep ${fname})
echo "${tmp_path}"
EOF

tmp_path=$(ssh ${sftp_cmd} 'cd ms; tmp_path=$(find | $(grep ${fname})); echo $tmp_path')
#parent_path=${tmp_path::-${size}}
echo "${tmp_path}"
#echo "${parent_path}"

# This is how you run commands into sftp
# alternatively use -b batch_file.my flag
sftp -R 1 ${sftp_cmd} << EOF
cd ${parent_path}
get *.gz
EOF

#jabref_dir=$( ssh my_server 'jabref_exe=$(which jabref); jabref_dir=$(dirname $jabref_exe); java -jar $jabref_dir/../jabref.jar > /dev/null; echo $jabref_dir' )

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
