#!/bin/sh

for dirname in $( ls -d */ )
do

    cd $dirname
    cd *
    sh /nas/longleaf/home/kolbt/whingdingdilly/post_proc/post_proc.sh
    cd ../..

done
