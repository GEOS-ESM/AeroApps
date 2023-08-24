#!/bin/bash
# example commands to run the data preprocessing scripts
#

 export PATH=/gpfsm/dnb04/projects/p109/FIREX-AQ/sample_forecasts2/scripts:$PATH
 dirn=/gpfsm/dnb04/projects/p109/FIREX-AQ/sample_forecasts2/PrsLatLon 

# CAMS
# ----
 cams.py -v CAMS/20190501/z_cams_c_ecmf_20190501000000_prod_fc_merged_006.nc 

# RAQMS 
# -----
 raqms2p.py   -v -o $dirn/raqms.llp.20180807_12z.nc4    RAQMS/uwhyb_08_07_2018_12Z.chem.assim.firex.nc
#arqi2llp.py  -v -o $dirn/arqi.llp.2018093012_036.nc4     ARQI/2018093012_036_uconversion_particles.netcdf4.compressed
#cam2p.py     -v -o $dirn/cam.llp.20180808.nc4          CAM-Chem/forecast.cam.na.2018-08-08-00000.nc
#cam_waccm2p.py     -v -o $dirn/cam.llp.20190508.nc4          CAM-Chem/20190503/forecast.cam.na.2019-05-08-00000.nc
#wrf2llp.py   -v -o $dirn/wrfchem.llp.20190331_12z.nc4  WRFChem-GOCART/wrfout_d01_2019-03-31_12.conus.subsetted.nc

#alt_g2ctl 1907019000600 >hrrrx.ctl
#alt_gmp hrrrx.ctl 
#lats4d.sh  -v -o $dirn/hrrrx.llp.20190312_01z.nc4       -i ./HRRRX/hrrrx.ctl -format netcdf4      # 3D


#ncl 'file_in="RAP-Chem/wrfout_d01_2018-08-07_06:00:00"' 'file_out="RAP-Chem/wrfout_d01_2018-08-07_06:00:00.subsetted.nc"' scripts/rapchem_subset.ncl 
#rapchem2llp.py   -v -o $dirn/rapchem.llp.20180807_06z.nc4  RAP-Chem/wrfout_d01_2018-08-07_06:00:00.subsetted.nc

#ncl 'file_in="NAAPS/2019031200_006_000_nap043o_000.nc"' 'file_out="PrsLatLon/naaps.llp.2019031200_006.nc"' scripts/naaps_subset.ncl  

#ncl 'file_in="NCAR-WRFchem/20190422/wrfout_d01_2019-04-22_06:00:00"' 'file_out="NCAR-WRFchem/20190422/wrfout_d01_2019-04-22_06:00:00.subsetted.nc"' scripts/ncar_wrfchem_subset.ncl 
#ncarwrfchem2llp.py -v -o $dirn/ncarwrfchem.llp.20190422_06z.nc4 NCAR-WRFchem/20190422/wrfout_d01_2019-04-22_06:00:00.subsetted.nc


#naqfc_split_hourly.sh 20190507 NAQFC-CMAQ/20190507
#naqfc_cmaq2llp.py  -v -o $dirn/cmaq.llp.20190507_12z.nc4      NAQFC-CMAQ/20190507/aqm.20190507.t12z.conc-select.20190507_12z.nc



##lats4d.sh    -v -o $dirn/cmaq.sfc.2080707_13z       -i NAQFC-CMAQ/cmaq-aot.ctl -format netcdf4  # 2D 
##lats4d.sh    -v -o $dirn/hrrr.sfc.2080707_07z       -i HRRR-Smoke/hrrr.ctl -format netcdf4      # 2D 


