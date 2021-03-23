rm *.o *.mod
exit

gfortran -c -Wall ../vlidort_def/vlidort_pars.f90
gfortran -c -Wall ../vlidort_def/vlidort_inputs_def.f90
gfortran -c -Wall ../vlidort_def/vlidort_outputs_def.f90
gfortran -c -Wall ../vlidort_def/vlidort_setups_def.f90
#exit

gfortran -c -Wall FO_geometry_Generic.f90   
gfortran -c -Wall FO_WPgeometry_Routines.f90
gfortran -c -Wall FO_DTWPgeometry_master.f90

gfortran -c -Wall FO_SSWPgeometry_master.f90
gfortran -c -Wall FO_VectorSS_Spherfuncs.f90

gfortran -c -Wall VFO_Geometry_Master.f90

gfortran -c -Wall FO_ThermalDTRT.f90
gfortran -c -Wall FO_ThermalDTRT_ILCS.f90
gfortran -c -Wall FO_ThermalDTRT_ILPS.f90

gfortran -c -Wall FO_VectorSSRT.f90
gfortran -c -Wall FO_VectorSSRT_ILCS.f90
gfortran -c -Wall FO_VectorSSRT_ILPS.f90

gfortran -c -Wall VFO_RTCalc_Master.f90
gfortran -c -Wall VFO_RTCalc_LinMasters.f90

exit

