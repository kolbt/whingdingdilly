#!bin/sh

path='/Users/kolbt/Desktop/compiled/whingdingdilly/gnuplot/'
xa=0
pb=0

#while [ $xa -le 100 ]
#do
#
#    cp ${path}MSD_template_gnuplot.txt MSD_gnuplot_xa${xa}.txt
#    cp ${path}MCS_template_gnuplot.txt MCS_gnuplot_xa${xa}.txt
#
#    for file in $(ls LOG*MCS*xa${xa}*.txt)
#    do
#
#        echo "\'$file\' using 1:2" >> MCS_gnuplot_xa${xa}.txt
#
#    done
#
#    for file in $(ls *MSD*ten*xa${xa}*.txt)
#    do
#
#        echo "\'$file\' using 1:2" >> MSD_gnuplot_xa${xa}.txt
#
#    done
#
#    xa=$(( $xa + 10 ))
#
#done

while [ $xa -le 100 ]
do

    cp ${path}MSD_temp.txt MSD_gnuplot_xa${xa}.txt
    gsed -i 's/\${xa}/'"${xa}"'/g' MSD_gnuplot_xa${xa}.txt
    gsed -i 's/\${pb}/'"${pb}"'/g' MSD_gnuplot_xa${xa}.txt
    xa=$(( $xa + 10 ))

done
