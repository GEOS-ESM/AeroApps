#!/bin/sh

idate=$1

ML="p3089 pm2p5_3d tcom tcpm2p5"
PL="aermr01 aermr02 aermr03 aermr04 aermr05 aermr06 aermr07 aermr08 aermr09 aermr10 aermr11 co no no2 so2"
SFC="aod550 pm2p5 bcaod550 duaod550 omaod550 pm1 pm10 ssaod550 suaod550 tc_no tc_pan tcco tcno2"

i=0
while [ $i -le 120 ]; do

  tau=`expr $i + 1000 | cut -c2-4`

  for var in $PL; do
    file=z_cams_c_ecmf_${idate}000000_prod_fc_pl_${tau}_$var.nc 
    if [ ! -f "$file" ]; then echo $file; fi
  done

  for var in $SFC; do
    file=z_cams_c_ecmf_${idate}000000_prod_fc_sfc_${tau}_$var.nc
    if [ ! -f "$file" ]; then echo $file; fi
  done

# for var in $ML; do
#   file=z_cams_c_ecmf_${idate}000000_prod_fc_ml_${tau}_$var.nc
#   if [ ! -f "$file" ]; then echo $file; fi
# done

  i=`expr $i + 3`

done

exit 0