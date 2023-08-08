#!/bin/sh

if [ $# -ge 1 ]; then
  idate=$1
else
  idate=`date "+%Y%m%d"`
fi

firex_plot.py $idate 0 models.rc models/GEOSFP.yml
firex_plot.py $idate 0 models.rc models/GEOSCF.yml
firex_plot.py $idate 0 models.rc models/CAMS.yml
firex_plot.py $idate 0 models.rc models/CAM-Chem.yml
firex_plot.py $idate 0 models.rc models/RAQMS.yml
firex_plot.py $idate 0 models.rc models/FireWork.yml
firex_plot.py $idate 0 models.rc models/ALL.yml
firex_plot.py $idate 0 models.rc models/HRRR-Smoke.yml
firex_plot.py $idate 0 models.rc models/WRFChem.yml
firex_plot.py $idate 0 models.rc models/UCLAWRFchem.yml
firex_plot.py $idate 0 models.rc models/UIOWAWRFchem.yml
firex_plot.py $idate 0 models.rc models/NCARWRFchem.yml
firex_plot.py $idate 0 models.rc models/ARQI.yml

exit 0
