#!bin/sh

path1='/Users/kolbt/Desktop/compiled/whingdingdilly/gnuplot/'
path2=$(pwd)
xa=50
pb=500
pa=0

while [ $pa -le 500 ]
do

    rep_name="all_pa${pa}_pb${pb}_xa${xa}"

    cp ${path1}gnuplot_all.txt stats_pa${pa}.txt
    gsed -i 's/\${iname}/'"${rep_name}"'/g' stats_pa${pa}.txt

    gnuplot ${path2}/stats_pa${pa}.txt

    pa=$(( $pa + 50 ))

done

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

#while [ $xa -le 100 ]
#do
#
#    cp ${path}MSD_temp.txt MSD_gnuplot_xa${xa}.txt
#    gsed -i 's/\${xa}/'"${xa}"'/g' MSD_gnuplot_xa${xa}.txt
#    gsed -i 's/\${pb}/'"${pb}"'/g' MSD_gnuplot_xa${xa}.txt
#
#    cp ${path}MSD_temp_dense.txt MSD_dense_gnuplot_xa${xa}.txt
#    gsed -i 's/\${xa}/'"${xa}"'/g' MSD_dense_gnuplot_xa${xa}.txt
#    gsed -i 's/\${pb}/'"${pb}"'/g' MSD_dense_gnuplot_xa${xa}.txt
#
#    cp ${path}MSD_temp_dilute.txt MSD_dilute_gnuplot_xa${xa}.txt
#    gsed -i 's/\${xa}/'"${xa}"'/g' MSD_dilute_gnuplot_xa${xa}.txt
#    gsed -i 's/\${pb}/'"${pb}"'/g' MSD_dilute_gnuplot_xa${xa}.txt
#
#    cp ${path}MSD_temp_dense_typed.txt MSD_dense_typed_gnuplot_xa${xa}.txt
#    gsed -i 's/\${xa}/'"${xa}"'/g' MSD_dense_typed_gnuplot_xa${xa}.txt
#    gsed -i 's/\${pb}/'"${pb}"'/g' MSD_dense_typed_gnuplot_xa${xa}.txt
#
#    cp ${path}MSD_temp_dilute_typed.txt MSD_dilute_typed_gnuplot_xa${xa}.txt
#    gsed -i 's/\${xa}/'"${xa}"'/g' MSD_dilute_typed_gnuplot_xa${xa}.txt
#    gsed -i 's/\${pb}/'"${pb}"'/g' MSD_dilute_typed_gnuplot_xa${xa}.txt
#
#    cp ${path}MCS_temp.txt MCS_gnuplot_xa${xa}.txt
#    gsed -i 's/\${xa}/'"${xa}"'/g' MCS_gnuplot_xa${xa}.txt
#    gsed -i 's/\${pb}/'"${pb}"'/g' MCS_gnuplot_xa${xa}.txt
#
#    xa=$(( $xa + 10 ))
#
#done

