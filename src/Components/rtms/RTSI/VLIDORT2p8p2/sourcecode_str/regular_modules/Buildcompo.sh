#rm  *.o  *.mod
#exit

#gfortran -c -Wall vlidort_lp_converge.f90
#gfortran -c -Wall vlidort_ls_converge.f90
#gfortran -c -Wall vlidort_rtcalc_lcsmaster.f90
#gfortran -c -Wall vlidort_rtcalc_lpsmaster.f90

#gfortran -c -Wall vlidort_setup_master.f90
#gfortran -c -Wall vlidort_miscsetups.f90
#gfortran -c -Wall vlidort_vfo_interface.f90
gfortran -c -Wall vlidort_rtcalc_master.f90
exit

gfortran -c -Wall ../../vlidort_def/vlidort_pars.f90
gfortran -c -Wall ../../vlidort_def/vlidort_inputs_def.f90
gfortran -c -Wall ../../vlidort_def/vlidort_outputs_def.f90
gfortran -c -Wall ../../vlidort_def/vlidort_sup_sleave_def.f90
gfortran -c -Wall ../../vlidort_def/vlidort_sup_brdf_def.f90
gfortran -c -Wall ../../vlidort_def/vlidort_sup_ss_def.f90
gfortran -c -Wall ../../vlidort_def/vlidort_sup_def.f90
gfortran -c -Wall ../../vlidort_def/vlidort_work_def.f90
gfortran -c -Wall ../../vlidort_def/vlidort_setups_def.f90
gfortran -c -Wall ../../vlidort_def/vlidort_io_defs.f90

gfortran -c -Wall ../../fo_main_1p5_NEW/FO_geometry_Generic.f90   
gfortran -c -Wall ../../fo_main_1p5_NEW/FO_WPgeometry_Routines.f90
gfortran -c -Wall ../../fo_main_1p5_NEW/FO_DTWPgeometry_master.f90
gfortran -c -Wall ../../fo_main_1p5_NEW/FO_SSWPgeometry_master.f90
gfortran -c -Wall ../../fo_main_1p5_NEW/FO_VectorSS_Spherfuncs.f90
gfortran -c -Wall ../../fo_main_1p5_NEW/FO_VectorSSRT.f90
gfortran -c -Wall ../../fo_main_1p5_NEW/FO_ThermalDTRT.f90
gfortran -c -Wall ../../fo_main_1p5_NEW/VFO_Geometry_Master.f90
gfortran -c -Wall ../../fo_main_1p5_NEW/VFO_RTCalc_Master.f90

gfortran -c -Wall lapack_tools.f90
gfortran -c -Wall vlidort_aux.f90
gfortran -c -Wall vlidort_inputs.f90
gfortran -c -Wall vlidort_Taylor.f90
gfortran -c -Wall vlidort_solutions.f90
gfortran -c -Wall vlidort_miscsetups.f90
gfortran -c -Wall vlidort_multipliers.f90
gfortran -c -Wall vlidort_thermalsup.f90
gfortran -c -Wall vlidort_writemodules.f90
gfortran -c -Wall vlidort_bvproblem.f90
gfortran -c -Wall vlidort_intensity.f90
gfortran -c -Wall vlidort_converge.f90
gfortran -c -Wall vlidort_mediaprops.f90
gfortran -c -Wall vlidort_pack.f90
gfortran -c -Wall vlidort_unpack.f90
gfortran -c -Wall vlidort_vfo_interface.f90
gfortran -c -Wall vlidort_setup_master.f90
gfortran -c -Wall vlidort_rtcalc_master.f90

exit

