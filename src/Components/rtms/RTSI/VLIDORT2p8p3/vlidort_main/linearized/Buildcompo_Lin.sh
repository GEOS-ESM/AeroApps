rm  *.o  *.mod
exit

#  MASTERS
#gfortran -c -Wall vlidort_lcs_masters_V2p8p3.f90
#gfortran -c -Wall vlidort_lps_masters_V2p8p3.f90
#exit

#  DEFS

gfortran -c -Wall ../../vlidort_def/vlidort_pars.f90
gfortran -c -Wall ../../vlidort_def/vlidort_inputs_def.f90
gfortran -c -Wall ../../vlidort_def/vlidort_outputs_def.f90
gfortran -c -Wall ../../vlidort_def/vlidort_sup_sleave_def.f90
gfortran -c -Wall ../../vlidort_def/vlidort_sup_brdf_def.f90
gfortran -c -Wall ../../vlidort_def/vlidort_sup_ss_def.f90
gfortran -c -Wall ../../vlidort_def/vlidort_sup_def.f90
gfortran -c -Wall ../../vlidort_def/vlidort_work_def.f90
gfortran -c -Wall ../../vlidort_def/vlidort_io_defs.f90

gfortran -c -Wall ../../vlidort_def/vlidort_lin_inputs_def.f90
gfortran -c -Wall ../../vlidort_def/vlidort_lin_outputs_def.f90
gfortran -c -Wall ../../vlidort_def/vlidort_lin_sup_sleave_def.f90
gfortran -c -Wall ../../vlidort_def/vlidort_lin_sup_brdf_def.f90
gfortran -c -Wall ../../vlidort_def/vlidort_lin_sup_ss_def.f90
gfortran -c -Wall ../../vlidort_def/vlidort_lin_sup_def.f90
gfortran -c -Wall ../../vlidort_def/vlidort_lin_work_def.f90
gfortran -c -Wall ../../vlidort_def/vlidort_lin_io_defs.f90

#exit

#  FO CODE

gfortran -c -Wall ../../vlidort_focode/FO_geometry_Generic.f90   
gfortran -c -Wall ../../vlidort_focode/FO_WPgeometry_Routines.f90
gfortran -c -Wall ../../vlidort_focode/FO_DTWPgeometry_master.f90
gfortran -c -Wall ../../vlidort_focode/FO_SSWPgeometry_master.f90
gfortran -c -Wall ../../vlidort_focode/FO_VectorSS_Spherfuncs.f90

gfortran -c -Wall ../../vlidort_focode/FO_VectorSS_RTCalcs_I.f90
gfortran -c -Wall ../../vlidort_focode/FO_Thermal_RTCalcs_I.f90
gfortran -c -Wall ../../vlidort_focode/VFO_Master.f90

gfortran -c -Wall ../../vlidort_focode/FO_VectorSS_RTCalcs_ILCS.f90
gfortran -c -Wall ../../vlidort_focode/FO_Thermal_RTCalcs_ILCS.f90
gfortran -c -Wall ../../vlidort_focode/FO_VectorSS_RTCalcs_ILPS.f90
gfortran -c -Wall ../../vlidort_focode/FO_Thermal_RTCalcs_ILPS.f90
gfortran -c -Wall ../../vlidort_focode/VFO_LinMasters.f90

#exit

#  REGULAR

gfortran -c -Wall ../regular/lapack_tools.f90
gfortran -c -Wall ../regular/vlidort_aux.f90
gfortran -c -Wall ../regular/vlidort_inputs.f90
gfortran -c -Wall ../regular/vlidort_Taylor.f90
gfortran -c -Wall ../regular/vlidort_geometry.f90
gfortran -c -Wall ../regular/vlidort_solutions.f90
gfortran -c -Wall ../regular/vlidort_miscsetups.f90
gfortran -c -Wall ../regular/vlidort_thermalsup.f90
gfortran -c -Wall ../regular/vlidort_writemodules.f90
gfortran -c -Wall ../regular/vlidort_bvproblem.f90
gfortran -c -Wall ../regular/vlidort_PostProcessing.f90
gfortran -c -Wall ../regular/vlidort_intensity.f90
gfortran -c -Wall ../regular/vlidort_converge.f90
gfortran -c -Wall ../regular/vlidort_mediaprops.f90
gfortran -c -Wall ../regular/vlidort_pack.f90
gfortran -c -Wall ../regular/vlidort_unpack.f90
gfortran -c -Wall ../regular/vlidort_vfo_interface.f90
gfortran -c -Wall ../regular/vlidort_masters_V2p8p3.f90

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
gfortran -c -Wall vlidort_lc_PostProcessing.f90
gfortran -c -Wall vlidort_lc_wfatmos.f90
gfortran -c -Wall vlidort_lcs_converge.f90

gfortran -c -Wall vlidort_lp_solutions.f90
gfortran -c -Wall vlidort_lp_mediaprops.f90
gfortran -c -Wall vlidort_lp_miscsetups.f90
gfortran -c -Wall vlidort_lp_bvproblem.f90
gfortran -c -Wall vlidort_lp_pack.f90
gfortran -c -Wall vlidort_lp_unpack.f90
gfortran -c -Wall vlidort_lp_PostProcessing.f90
gfortran -c -Wall vlidort_lp_wfatmos.f90
gfortran -c -Wall vlidort_lps_converge.f90

gfortran -c -Wall vlidort_ls_wfsleave.f90
gfortran -c -Wall vlidort_ls_wfsurface.f90

#  TOP LEVEL

gfortran -c -Wall vlidort_vfo_lcs_interface.f90
gfortran -c -Wall vlidort_vfo_lps_interface.f90

#exit

gfortran -c -Wall vlidort_lcs_masters_V2p8p3.f90
gfortran -c -Wall vlidort_lps_masters_V2p8p3.f90

exit
