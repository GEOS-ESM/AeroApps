rm  *.o  *.mod
exit
gfortran -c -Wall vlidort_vfo_lcs_interface.f90
gfortran -c -Wall vlidort_rtcalc_lcsmaster.f90
gfortran -c -Wall vlidort_vfo_lps_interface.f90
gfortran -c -Wall vlidort_rtcalc_lpsmaster.f90

exit

#  DEFS

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

gfortran -c -Wall ../../vlidort_def/vlidort_lin_inputs_def.f90
gfortran -c -Wall ../../vlidort_def/vlidort_lin_outputs_def.f90
gfortran -c -Wall ../../vlidort_def/vlidort_lin_sup_sleave_def.f90
gfortran -c -Wall ../../vlidort_def/vlidort_lin_sup_brdf_def.f90
gfortran -c -Wall ../../vlidort_def/vlidort_lin_sup_ss_def.f90
gfortran -c -Wall ../../vlidort_def/vlidort_lin_sup_def.f90
gfortran -c -Wall ../../vlidort_def/vlidort_lin_work_def.f90
gfortran -c -Wall ../../vlidort_def/vlidort_lin_io_defs.f90

#  FO CODE

gfortran -c -Wall ../../fo_main_1p5_NEW/FO_geometry_Generic.f90   
gfortran -c -Wall ../../fo_main_1p5_NEW/FO_WPgeometry_Routines.f90
gfortran -c -Wall ../../fo_main_1p5_NEW/FO_DTWPgeometry_master.f90
gfortran -c -Wall ../../fo_main_1p5_NEW/FO_SSWPgeometry_master.f90
gfortran -c -Wall ../../fo_main_1p5_NEW/FO_VectorSS_Spherfuncs.f90
gfortran -c -Wall ../../fo_main_1p5_NEW/VFO_Geometry_Master.f90

gfortran -c -Wall ../../fo_main_1p5_NEW/FO_VectorSSRT.f90
gfortran -c -Wall ../../fo_main_1p5_NEW/FO_ThermalDTRT.f90
gfortran -c -Wall ../../fo_main_1p5_NEW/VFO_RTCalc_Master.f90

gfortran -c -Wall ../../fo_main_1p5_NEW/FO_VectorSSRT_ILCS.f90
gfortran -c -Wall ../../fo_main_1p5_NEW/FO_ThermalDTRT_ILCS.f90
gfortran -c -Wall ../../fo_main_1p5_NEW/FO_VectorSSRT_ILPS.f90
gfortran -c -Wall ../../fo_main_1p5_NEW/FO_ThermalDTRT_ILPS.f90
gfortran -c -Wall ../../fo_main_1p5_NEW/VFO_RTCalc_LinMasters.f90

#  REGULAR

gfortran -c -Wall ../regular_modules/lapack_tools.f90
gfortran -c -Wall ../regular_modules/vlidort_aux.f90
gfortran -c -Wall ../regular_modules/vlidort_inputs.f90
gfortran -c -Wall ../regular_modules/vlidort_Taylor.f90
gfortran -c -Wall ../regular_modules/vlidort_solutions.f90
gfortran -c -Wall ../regular_modules/vlidort_miscsetups.f90
gfortran -c -Wall ../regular_modules/vlidort_multipliers.f90
gfortran -c -Wall ../regular_modules/vlidort_thermalsup.f90
gfortran -c -Wall ../regular_modules/vlidort_writemodules.f90
gfortran -c -Wall ../regular_modules/vlidort_bvproblem.f90
gfortran -c -Wall ../regular_modules/vlidort_intensity.f90
gfortran -c -Wall ../regular_modules/vlidort_converge.f90
gfortran -c -Wall ../regular_modules/vlidort_mediaprops.f90
gfortran -c -Wall ../regular_modules/vlidort_pack.f90
gfortran -c -Wall ../regular_modules/vlidort_unpack.f90
gfortran -c -Wall ../regular_modules/vlidort_vfo_interface.f90
gfortran -c -Wall ../regular_modules/vlidort_setup_master.f90
gfortran -c -Wall ../regular_modules/vlidort_rtcalc_master.f90

#  Linearized

gfortran -c -Wall vlidort_lpc_solutions.f90
gfortran -c -Wall vlidort_lpc_bvproblem.f90
gfortran -c -Wall vlidort_la_miscsetups.f90
gfortran -c -Wall vlidort_l_inputs.f90
gfortran -c -Wall vlidort_l_pack.f90
gfortran -c -Wall vlidort_l_unpack.f90
gfortran -c -Wall vlidort_l_writemodules.f90
gfortran -c -Wall vlidort_l_thermalsup.f90
gfortran -c -Wall vlidort_lbbf_jacobians_vector.f90

gfortran -c -Wall vlidort_lc_solutions.f90
gfortran -c -Wall vlidort_lc_mediaprops.f90
gfortran -c -Wall vlidort_lc_miscsetups.f90
gfortran -c -Wall vlidort_lc_bvproblem.f90
gfortran -c -Wall vlidort_lc_pack.f90
gfortran -c -Wall vlidort_lc_unpack.f90
gfortran -c -Wall vlidort_lc_wfatmos.f90
gfortran -c -Wall vlidort_lcs_converge.f90

gfortran -c -Wall vlidort_lp_solutions.f90
gfortran -c -Wall vlidort_lp_mediaprops.f90
gfortran -c -Wall vlidort_lp_miscsetups.f90
gfortran -c -Wall vlidort_lp_bvproblem.f90
gfortran -c -Wall vlidort_lp_pack.f90
gfortran -c -Wall vlidort_lp_unpack.f90
gfortran -c -Wall vlidort_lp_wfatmos.f90
gfortran -c -Wall vlidort_lps_converge.f90

gfortran -c -Wall vlidort_ls_wfsleave.f90
gfortran -c -Wall vlidort_ls_wfsurface.f90

#  TOP LEVEL

gfortran -c -Wall vlidort_vfo_lcs_interface.f90
gfortran -c -Wall vlidort_vfo_lps_interface.f90

exit

gfortran -c -Wall vlidort_rtcalc_lcsmaster.f90
gfortran -c -Wall vlidort_rtcalc_lpsmaster.f90

exit
