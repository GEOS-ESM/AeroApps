#!/bin/bash  

# usage: append_vars.sh yyyymmdd 
#

runtime=$1
y4=`echo $runtime | cut -c1-4`
y2=`echo $runtime | cut -c3-4`
m2=`echo $runtime | cut -c5-6`
d2=`echo $runtime | cut -c7-8`
path_store="/home/xinxinye/data/FIREX-AQ/compare/CAMS/cropped"
path_data=$path_store"/"$y4$m2$d2

cd $path_data
for fh in `seq -w 0 3 120`
do
  i=0
  for fname in `ls *sfc_${fh}*`
  do
    let i+=1
    if [ $i -eq '1' ];then
       ss=${fname/sfc/merged}
       fout=${ss/_aod1240/}
       echo $fout
       cp $fname $fout
    else
       ncks -h -A $fname $fout
    fi
  done    

  for fname in `ls *pl_${fh}*`
  do 
       ncks -h -A $fname $fout
  done

  for fname in `ls *ml_${fh}*`
  do 
       ncks -h -A $fname $fout
  done
done 


