#!/bin/sh

exec_path="/Users/kolbt/Desktop/compiled/whingdingdilly/developing"
python $exec_path/write_new.py peB150.png                       # create blank image
x=0
y=180

for filename in $( ls final*.png )
do

    python $exec_path/crop_place.py $filename peB150.png $x $y  # populate image with snaps

    x=$(( $x + 180 ))                                           # get next coords for paste
    if [ $x -eq 2880 ]; then                                    # check if at image edge
        x=$(( 0 ))                                              # move back to left edge
        y=$(( $y + 180 ))                                       # move to next row
    fi

done

# you have to hand the placement file a position from this loop
# I'll have to pass the image files in a more intelligent way
# (so that everything is pasted in it's rightful spot)
# also, no reason that I couldn't move up columns
