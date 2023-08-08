#!/bin/sh  
# split the full outputs for all forecast hours to hourly files
# append 3-D pressure to the chemical files
# Usage :  yyyymmdd path_data
# -------------------------------------------------------------


runtime=$1
y4=`echo $runtime | cut -c1-4`
y2=`echo $runtime | cut -c3-4`
m2=`echo $runtime | cut -c5-6`
d2=`echo $runtime | cut -c7-8`
h2=12 

path_store=$2

cd $path_store

for dh in `seq 0 1 72`
do
    dtime=`date -d "$y4$m2$d2 $h2 $dh hour" +"%Y%m%d_%Hz"`
    ncks -O -d TSTEP,$dh,$dh ./aqm.$runtime.t12z.conc-select.ncf ./aqm.$runtime.t12z.conc-select.$dtime.nc
    ncks -A -d TSTEP,$dh,$dh -d LAY,0,19 -v PRES,QV,QC,ZF,TA ./aqm.$runtime.t12z.metcro3d-select.ncf ./aqm.$runtime.t12z.conc-select.$dtime.nc 
    #ncks -A -d TSTEP,$dh,$dh -d LAY,0,19 -v UWIND,VWIND ./aqm.$runtime.t12z.metdot3d-select.ncf ./aqm.$runtime.t12z.conc-select.$dtime.nc 
    ncks -A -d TSTEP,$dh,$dh -d LAY,0 -v PBL2 ./aqm.$runtime.t12z.metcro2d-select.ncf ./aqm.$runtime.t12z.conc-select.$dtime.nc 
done 

