#!/bin/sh

for tarname in $( ls *.gz )
do

    tar -xf ${tarname}
    rm ${tarname}

done
