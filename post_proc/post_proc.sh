#!/bin/sh

for filename in $( ls *.gsd )
do

    sh analyze.sh $filename

done
