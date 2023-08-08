#!/bin/sh

OPTIONS='--lights_off'

#### GEOS-5 ####
#wxmap.py --config `pwd` --stream GEOSFP --plot aod --oname '$stream.$plot.$level.%Y%m%dT%H%M.png' --region usa --fcst_dt 20190321 --time_dt 20190322 $OPTIONS

#wxmap.py --config `pwd` --stream GEOSAN --plot aod --oname '$stream.$plot.$level.%Y%m%dT%H%M.png' --region usa --time_dt 20180807T12 $OPTIONS

#### CAMS #####
for var in  'tcpm25' 'tcoc' 'pm25' 'aod' 'pm25sfc' 'tcco' 'tcno' 'su' 'oc' 'bc' 'du' 'ss' 'co' 'no2' 'no1' 'so2'
do 
#wxmap.py --config `pwd` --stream CAMS --plot $var  --oname '$stream.$plot.$level.%Y%m%dT%H%M.png' --region usa2 --fcst_dt 20180902T00 --time_dt 20180902T06 $OPTIONS
wxmap.py --config `pwd` --stream CAMS --plot $var  --oname '$stream.$plot.$level.%Y%m%dT%H%M.png' --region usa2 --fcst_dt 20190501T00 --time_dt 20190501T06 $OPTIONS
done

for var in 'no2' 'su' 'oc' 'bc' 'du' 'ss' 'co' 'no2' 'no1' 'so2'
do 
##wxmap.py --config `pwd` --stream CAMS --plot ${var}_lml  --oname '$stream.$plot.$level.%Y%m%dT%H%M.png' --region usa2 --fcst_dt 20180902T00 --time_dt 20180902T06 $OPTIONS
wxmap.py --config `pwd` --stream CAMS --plot ${var}_lml  --oname '$stream.$plot.$level.%Y%m%dT%H%M.png' --region usa2 --fcst_dt 20190501T00 --time_dt 20190501T06 $OPTIONS
#done 

#### RAQMS ####
for var in 'tcoc' 'aod' 'su' 'oc' 'bc' 'du' 'ss' 'co' 'no1' 'no2' 'o3' 'so2' 'tcco' 
do 
wxmap.py --config `pwd` --stream RAQMS --plot $var --oname '$stream.$plot.$level.%Y%m%dT%H%M.png' --region usa2 --fcst_dt 20180807T12 --time_dt 20180807T12 $OPTIONS
done 

for var in 'su' 'oc' 'bc' 'du' 'ss' 'co' 'no1' 'no2' 'o3' 'so2'
do 
wxmap.py --config `pwd` --stream RAQMS --plot ${var}_lml --oname '$stream.$plot.$level.%Y%m%dT%H%M.png' --region usa2 --fcst_dt 20180807T12 --time_dt 20180807T12 $OPTIONS
done 

#### CAM Chem ####
for var in 'co' 'co_smoke' 'am' 'su' 'oc' 'bc' 'du' 'ss' 'no1' 'no2' 'o3' 
do 
wxmap.py --config `pwd` --stream CAMchem --plot $var  --oname '$stream.$plot.$level.%Y%m%dT%H%M.png' --region usa2 --fcst_dt 20180808T00 --time_dt 20180808T00 $OPTIONS
done 

for var in 'co_smoke' 'am' 'su' 'oc' 'bc' 'du' 'ss' 'co' 'co_smoke' 'no1' 'no2' 'o3' 
do 
wxmap.py --config `pwd` --stream CAMchem --plot ${var}_lml  --oname '$stream.$plot.$level.%Y%m%dT%H%M.png' --region usa2 --fcst_dt 20180808T00 --time_dt 20180808T00 $OPTIONS
done 

#### WRF Chem ####
#wxmap.py --config `pwd` --stream WRFchem --plot aod --oname '$stream.$plot.$level.%Y%m%dT%H%M.png' --region usa --fcst_dt 20180807T12 --time_dt 20180807T12 $OPTIONS
#wxmap.py --config `pwd` --stream WRFchem --plot pm25pl --oname '$stream.$plot.$level.%Y%m%dT%H%M.png' --region usa --fcst_dt 20180807T12 --time_dt 20180807T12 $OPTIONS
#wxmap.py --config `pwd` --stream WRFchem --plot pm25sfc --oname '$stream.$plot.$level.%Y%m%dT%H%M.png' --region usa --fcst_dt 20180807T12 --time_dt 20180807T12 $OPTIONS

#### HRRR #### 
#wxmap.py --config `pwd` --stream HRRR --plot aod --oname '$stream.$plot.$level.%Y%m%dT%H%M.png' --region usa --fcst_dt 20180807T12 --time_dt 20180808T07 $OPTIONS
#wxmap.py --config `pwd` --stream HRRRX --plot aod --oname '$stream.$plot.$level.%Y%m%dT%H%M.png' --region usa --fcst_dt 20190311T19 --time_dt 20190312T01 $OPTIONS
#wxmap.py --config `pwd` --stream HRRRX --plot pm25_smoke --oname '$stream.$plot.$level.%Y%m%dT%H%M.png' --region usa --fcst_dt 20190311T19 --time_dt 20190312T01 $OPTIONS
#wxmap.py --config `pwd` --stream RAPchem --plot aod --oname '$stream.$plot.$level.%Y%m%dT%H%M.png' --region usa --fcst_dt 20180807T00 --time_dt 20180807T06 $OPTIONS
#wxmap.py --config `pwd` --stream RAPchem --plot pm25pl --oname '$stream.$plot.$level.%Y%m%dT%H%M.png' --region usa --fcst_dt 20180807T00 --time_dt 20180807T06 $OPTIONS
#wxmap.py --config `pwd` --stream RAPchem --plot am --oname '$stream.$plot.$level.%Y%m%dT%H%M.png' --region usa --fcst_dt 20180807T00 --time_dt 20180807T06 $OPTIONS
#wxmap.py --config `pwd` --stream RAPchem --plot su --oname '$stream.$plot.$level.%Y%m%dT%H%M.png' --region usa --fcst_dt 20180807T00 --time_dt 20180807T06 $OPTIONS
#wxmap.py --config `pwd` --stream RAPchem --plot ni --oname '$stream.$plot.$level.%Y%m%dT%H%M.png' --region usa --fcst_dt 20180807T00 --time_dt 20180807T06 $OPTIONS
#wxmap.py --config `pwd` --stream RAPchem --plot oc --oname '$stream.$plot.$level.%Y%m%dT%H%M.png' --region usa --fcst_dt 20180807T00 --time_dt 20180807T06 $OPTIONS
#wxmap.py --config `pwd` --stream RAPchem --plot bc --oname '$stream.$plot.$level.%Y%m%dT%H%M.png' --region usa --fcst_dt 20180807T00 --time_dt 20180807T06 $OPTIONS
#wxmap.py --config `pwd` --stream RAPchem --plot du --oname '$stream.$plot.$level.%Y%m%dT%H%M.png' --region usa --fcst_dt 20180807T00 --time_dt 20180807T06 $OPTIONS
#wxmap.py --config `pwd` --stream RAPchem --plot ss --oname '$stream.$plot.$level.%Y%m%dT%H%M.png' --region usa --fcst_dt 20180807T00 --time_dt 20180807T06 $OPTIONS
#wxmap.py --config `pwd` --stream RAPchem --plot co --oname '$stream.$plot.$level.%Y%m%dT%H%M.png' --region usa --fcst_dt 20180807T00 --time_dt 20180807T06 $OPTIONS
#wxmap.py --config `pwd` --stream RAPchem --plot no2 --oname '$stream.$plot.$level.%Y%m%dT%H%M.png' --region usa --fcst_dt 20180807T00 --time_dt 20180807T06 $OPTIONS
#wxmap.py --config `pwd` --stream RAPchem --plot no1 --oname '$stream.$plot.$level.%Y%m%dT%H%M.png' --region usa --fcst_dt 20180807T00 --time_dt 20180807T06 $OPTIONS
#wxmap.py --config `pwd` --stream RAPchem --plot so2 --oname '$stream.$plot.$level.%Y%m%dT%H%M.png' --region usa --fcst_dt 20180807T00 --time_dt 20180807T06 $OPTIONS
#wxmap.py --config `pwd` --stream RAPchem --plot o3 --oname '$stream.$plot.$level.%Y%m%dT%H%M.png' --region usa --fcst_dt 20180807T00 --time_dt 20180807T06 $OPTIONS
#wxmap.py --config `pwd` --stream RAPchem --plot smoke --oname '$stream.$plot.$level.%Y%m%dT%H%M.png' --region usa --fcst_dt 20180807T00 --time_dt 20180807T06 $OPTIONS


#wxmap.py --config `pwd` --stream ARQI --plot aod  --oname '$stream.$plot.$level.%Y%m%dT%H%M.png' --region usa2 --fcst_dt 20180807T12 --time_dt 20180807T18 $OPTIONS

#wxmap.py --config `pwd` --stream NAAPS --plot so2  --oname '$stream.$plot.$level.%Y%m%dT%H%M.png' --region usa2 --fcst_dt 20190312T00 --time_dt 20190312T06 $OPTIONS
#wxmap.py --config `pwd` --stream NAAPS --plot su  --oname '$stream.$plot.$level.%Y%m%dT%H%M.png' --region usa2 --fcst_dt 20190312T00 --time_dt 20190312T06 $OPTIONS
#wxmap.py --config `pwd` --stream NAAPS --plot du  --oname '$stream.$plot.$level.%Y%m%dT%H%M.png' --region usa2 --fcst_dt 20190312T00 --time_dt 20190312T06 $OPTIONS
#wxmap.py --config `pwd` --stream NAAPS --plot ss  --oname '$stream.$plot.$level.%Y%m%dT%H%M.png' --region usa2 --fcst_dt 20190312T00 --time_dt 20190312T06 $OPTIONS

