dset /nfs3m/portal/GMAO/gds_work/Contents/Analysis/Eval/x0034/ensdiag/%e/x0034.inst3_3d_wxme_Np.%y4%m2%d2_%h200z.nc4
options template
title 2d,1-Hourly,Instantaneous,Single-Level,Assimilation,Single-Level Diagnostics
undef 1e+15
dtype netcdf
xdef 288 linear -180 1.25
ydef 181 linear -90 1.0
zdef 7 levels 850 700 600 500 400 300 200
tdef 500 linear 06Z01apr2018 6hr
edef 32 names mem001 mem002 mem003 mem004 mem005 mem006 mem007 mem008 mem009 mem010
              mem011 mem012 mem013 mem014 mem015 mem016 mem017 mem018 mem019 mem020
              mem021 mem022 mem023 mem024 mem025 mem026 mem027 mem028 mem029 mem030
              mem031 mem032
vars 20
BCEXTTAU=>bcexttau  0  t,y,x  Black Carbon Extinction AOT [550 nm]
CLDHGH=>cldhgh  0  t,y,x  cloud_area_fraction_for_high_clouds
CLDLOW=>cldlow  0  t,y,x  cloud_area_fraction_for_low_clouds
CLDMID=>cldmid  0  t,y,x  cloud_area_fraction_for_middle_clouds
CLOUD=>cloud  7  t,z,y,x  cloud_fraction_for_radiation
DU=>du  7  t,z,y,x  Dust Mass Mixing Ratio
DUEXTTAU=>duexttau  0  t,y,x  Dust Extinction AOT [550 nm]
EPV=>epv  7  t,z,y,x  ertels_potential_vorticity
OCEXTTAU=>ocexttau  0  t,y,x  Organic Carbon Extinction AOT [550 nm] __ENSEMBLE__
OMEGA=>omega  7  t,z,y,x  vertical_pressure_velocity
PRECTOT=>prectot  0  t,y,x  total_precipitation
RH=>rh  7  t,z,y,x  relative_humidity_after_moist
SLP=>slp  0  t,y,x  sea_level_pressure
SUEXTTAU=>suexttau  0  t,y,x  SO4 Extinction AOT [550 nm] __ENSEMBLE__
T=>t  7  t,z,y,x  air_temperature
TAUHGH=>tauhgh  0  t,y,x  in_cloud_optical_thickness_of_high_clouds(EXPORT)
TAULOW=>taulow  0  t,y,x  in_cloud_optical_thickness_of_low_clouds
TAUMID=>taumid  0  t,y,x  in_cloud_optical_thickness_of_middle_clouds
U=>u  7  t,z,y,x  eastward_wind
V=>v  7  t,z,y,x  northward_wind
endvars
