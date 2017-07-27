#!/bin/sh

# tar all directories in a given directory
for directories in $( ls -d */ )
do

    tarname=${directories%?}
    tar -cf ${tarname}.tar.gz ${directories}

done
