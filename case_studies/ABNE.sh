#!/bin/sh

current=$( date "+%m_%d_%y" )
this_path=$( pwd )

echo "Are you running on Longleaf (y/n)?"
read answer

if [ $answer == "y" ]; then
    script_path='/nas/longleaf/home/kolbt/whingdingdilly/run.sh'
    script_path='/nas/longleaf/home/kolbt/whingdingdilly/run_gpu.sh'
    slowOut='/nas/longleaf/home/kolbt/whingdingdilly/case_studies/binaryBubbleSlowOut.py'
    fastOut='/nas/longleaf/home/kolbt/whingdingdilly/case_studies/binaryBubbleFastOut.py'
    sedtype='sed'
    submit='sbatch'
else
    script_path='/Users/kolbt/Desktop/compiled/whingdingdilly/run.sh'
    template='/Users/kolbt/Desktop/compiled/whingdingdilly/case_studies/densityGradient.py'
    sedtype='gsed'
    submit='sh'
fi

# Lists for simulations to submit
pes=(0 50 50 150)
pef=(150 150 500 500)
phiG=(10 8 4 3)

mkdir ${current}_parent
cd ${current}_parent

count=$(( 0 ))
for i in ${pes[@]}
do
    
    # Get variables
    curPeS=${pes[${count}]}
    curPeF=${pef[${count}]}
    curPhi=${phiG[${count}]}
    
    # Write variables to outer layer slow file
    slowFile=slow_out_pes${curPeS}_pef${curPeF}.py
    $sedtype -e 's/\${slowAct}/'"${curPeS}"'/g' $slowOut > $slowFile
    $sedtype -i 's/\${fastAct}/'"${curPeF}"'/g' $slowFile
    $sedtype -i 's/\${phiG}/'"${curPhi}"'/g' $slowFile
    # Run slow file
    $submit $script_path $slowFile
    
    # Write to outer layer fast file
    fastFile=fast_out_pes${curPeS}_pef${curPeF}.py
    $sedtype -e 's/\${slowAct}/'"${curPeS}"'/g' $fastOut > $fastFile
    $sedtype -i 's/\${fastAct}/'"${curPeF}"'/g' $fastFile
    $sedtype -i 's/\${phiG}/'"${curPhi}"'/g' $fastFile
    # Run fast file
    $submit $script_path $fastFile

    # Increment counter
    count=$(( $count + 1 ))

done
