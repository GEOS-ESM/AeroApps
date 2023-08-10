#!/bin/bash
# Post-processing of model forecasts 
# Usage: cams.sh [ccyymmdd] [hhmmss] [tau] [input dir] [output dir]
# ccyymmdd:  initial date
# hhmmss:    initial time 
# tau:       forecast hour 
# input dir: data storage dir 
# output dir
# ---------------------------------------------------------------------------

export PATH=/gpfsm/dnb04/projects/p109/FIREX-AQ/sample_forecasts2/scripts:$PATH

idate=$1
itime=$2 
tau=`seq -w $3 120 120`

y1=`echo $idate | cut -c1-4`
m1=`echo $idate | cut -c5-6`
d1=`echo $idate | cut -c7-8`
h1=`echo $itime | cut -c1-2`
m1=`echo $itime | cut -c3-4`

fh=`expr $3 \* 100 | cut -c1-3`
dir_in=$4
dir_out=$5

if [ ! -d $dir_out ]; then
    mkdir $dir_out
fi

# append the files to a single file 
  i=0
  cd $dir_in 

  for fname in `ls *sfc_${fh}*`
  do
    let i+=1
    if [ $i -eq '1' ];then
       ss=${fname/sfc/merged}
       fmerged=${ss/_aod1240/}
       ncks -O -d latitude,15.0,70.0 -d longitude,190.0,310.0 $fname $fmerged 
    else
       ncks -h -A -O -d latitude,15.0,70.0 -d longitude,190.0,310.0 $fname $fmerged
    fi
  done    

  for fname in `ls *pl_${fh}*`
  do 
       ncks -h -A -O -d latitude,15.0,70.0 -d longitude,190.0,310.0 $fname $fmerged
  done

  for fname in `ls *ml_${fh}*`
  do 
       ncks -h -A -O -d latitude,15.0,70.0 -d longitude,190.0,310.0 $fname $fmerged
  done

# interpolate to llp coordinate
 fout=${fmerged/merged/prep}
 cams.py -v -o $dir_out/$fout  $dir_in/$fmerged 



