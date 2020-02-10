#!/bin/sh

# Rename files with epsilon
eps=0.1
for file in $(ls *txt)
do

    mv ${file} "cluster_${file}"

done
