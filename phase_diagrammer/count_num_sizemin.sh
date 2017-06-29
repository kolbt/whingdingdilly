#!/bin/sh

size_min=??? # maybe like 2000?
occurrences=0
occur_min=??? # maybe like 100 tsteps?
output=0

for i in "cat myfile.txt | awk '{print $1}'"
do
    if [ $i -ge $size_min ]; then
        occurrences=$(( $occurrences + 1 ))
    fi
done

if [ $occurrences -ge $occur_min ]; then
    output=$(( 1 ))

echo $output >> outfile.txt
