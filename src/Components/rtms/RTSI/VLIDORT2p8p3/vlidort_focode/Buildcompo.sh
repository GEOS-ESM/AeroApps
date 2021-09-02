rm *.o *.mod *~
exit

gfortran -c -Wall ../vlidort_def/vlidort_pars.f90

gfortran -c -Wall FO_geometry_Generic.f90
gfortran -c -Wall FO_Taylor.f90
gfortran -c -Wall FO_Planckfunc.f90
gfortran -c -Wall FO_VectorSS_Spherfuncs.f90
gfortran -c -Wall FO_WPgeometry_Routines.f90
gfortran -c -Wall FO_SSWPgeometry_master.f90
gfortran -c -Wall FO_DTWPgeometry_master.f90
#exit

gfortran -c -Wall FO_Thermal_RTCalcs_I.f90
gfortran -c -Wall FO_VectorSS_RTCalcs_I.f90
gfortran -c -Wall VFO_Master.f90
#exit

gfortran -c -Wall FO_Thermal_RTCalcs_ILCS.f90
gfortran -c -Wall FO_Thermal_RTCalcs_ILPS.f90
gfortran -c -Wall FO_VectorSS_RTCalcs_ILCS.f90
gfortran -c -Wall FO_VectorSS_RTCalcs_ILPS.f90
gfortran -c -Wall VFO_LinMasters.f90
exit
