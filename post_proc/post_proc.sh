#!/bin/sh

echo "Are you running on Longleaf (y/n)?"
read answer

if [ $answer == "y" ]; then
#    hoomd_path='/nas/longleaf/home/kolbt/programs/hoomd_2.2.1/hoomd-blue/build'
    hoomd_path='/nas/longleaf/home/kolbt/programs/hoomd-blue/build'
    gsd_path='/nas/longleaf/home/kolbt/programs/gsd/build'
    script_path='/nas/longleaf/home/kolbt/whingdingdilly/post_proc'
    sedtype='sed'
    submit='sbatch'
else
    hoomd_path='/Users/kolbt/Desktop/compiled/hoomd-blue_11_8_17/hoomd-blue/build'
    gsd_path='/Users/kolbt/Desktop/compiled/gsd/build'
    script_path='/Users/kolbt/Desktop/compiled/whingdingdilly/post_proc'
    sedtype='gsed'
    submit='sh'
fi

echo "Are you running on your laptop (y/n)?"
read answer

if [ $answer == "y" ]; then
    hoomd_path='/Users/kolbt/Desktop/compiled/hoomd-blue/build'
fi

echo "Epsilon in filename?"
read isEps

if [ $isEps == 'n' ]; then
    ep=$(( 1 ))
    for filename in $( ls pa*.gsd )
    #for filename in $( ls *pa*_0.png )
    #for filename in $( ls all_pa*.txt )
    do

        # pull parameters from filename

        pa=$(echo $filename | $sedtype 's/^.*pa\([0-9]*\)_.*/\1/')
        pb=$(echo $filename | $sedtype 's/^.*pb\([0-9]*\)_.*/\1/')
        xa=$(echo $filename | $sedtype 's/^.*xa\([0-9]*\)..*/\1/')

        $submit $script_path/analyze.sh $pa $pb $xa $hoomd_path $gsd_path $script_path $ep

    done

else
    for filename in $( ls pa*.gsd )
    do

        # pull parameters from filename

        pa=$(echo $filename | $sedtype 's/^.*pa\([0-9]*\)_.*/\1/')
        pb=$(echo $filename | $sedtype 's/^.*pb\([0-9]*\)_.*/\1/')
        xa=$(echo $filename | $sedtype 's/^.*xa\([0-9]*\)_.*/\1/')
        ep=$(echo $filename | $sedtype 's/^.*ep\([0-9]*\)..*/\1/')
#        phi=$(echo $filename | $sedtype 's/^.*phi\([0-9]*\)..*/\1/')
#        al=$(echo $filename | $sedtype 's/^.*alpha\([0-9]*\)..*/\1/')
#        num=$(echo $filename | $sedtype 's/^.*N\([0-9]*\)..*/\1/')

#        $submit $script_path/analyze.sh $pa $pb $xa $hoomd_path $gsd_path $script_path $ep $num
        $submit $script_path/analyze.sh $pa $pb $xa $hoomd_path $gsd_path $script_path $ep

    done
fi

#pa=100
#pb=500
#xa=50
#$submit $script_path/analyze.sh $pa $pb $xa $hoomd_path $gsd_path $script_path

