rm  *.o  *.mod
exit

#gfortran -c -Wall vlidort_vfo_interface.f90
#gfortran -c -Wall vlidort_masters_V2p8p3.f90
#exit

gfortran -c -Wall ../../vlidort_def/vlidort_pars.f90
gfortran -c -Wall ../../vlidort_def/vlidort_inputs_def.f90
gfortran -c -Wall ../../vlidort_def/vlidort_outputs_def.f90
gfortran -c -Wall ../../vlidort_def/vlidort_sup_sleave_def.f90
gfortran -c -Wall ../../vlidort_def/vlidort_sup_brdf_def.f90
gfortran -c -Wall ../../vlidort_def/vlidort_sup_ss_def.f90
gfortran -c -Wall ../../vlidort_def/vlidort_sup_def.f90
gfortran -c -Wall ../../vlidort_def/vlidort_work_def.f90
gfortran -c -Wall ../../vlidort_def/vlidort_io_defs.f90

gfortran -c -Wall ../../vlidort_focode/FO_geometry_Generic.f90   
gfortran -c -Wall ../../vlidort_focode/FO_WPgeometry_Routines.f90
gfortran -c -Wall ../../vlidort_focode/FO_DTWPgeometry_master.f90
gfortran -c -Wall ../../vlidort_focode/FO_SSWPgeometry_master.f90
gfortran -c -Wall ../../vlidort_focode/FO_VectorSS_Spherfuncs.f90
gfortran -c -Wall ../../vlidort_focode/FO_VectorSS_RTCalcs_I.f90
gfortran -c -Wall ../../vlidort_focode/FO_Thermal_RTCalcs_I.f90
gfortran -c -Wall ../../vlidort_focode/VFO_Master.f90

gfortran -c -Wall lapack_tools.f90
gfortran -c -Wall vlidort_aux.f90
gfortran -c -Wall vlidort_inputs.f90
gfortran -c -Wall vlidort_Taylor.f90
gfortran -c -Wall vlidort_geometry.f90
gfortran -c -Wall vlidort_getplanck.f90
gfortran -c -Wall vlidort_solutions.f90
gfortran -c -Wall vlidort_miscsetups.f90
gfortran -c -Wall vlidort_bvproblem.f90
gfortran -c -Wall vlidort_PostProcessing.f90
gfortran -c -Wall vlidort_intensity.f90
gfortran -c -Wall vlidort_converge.f90
gfortran -c -Wall vlidort_thermalsup.f90
gfortran -c -Wall vlidort_writemodules.f90
gfortran -c -Wall vlidort_mediaprops.f90
gfortran -c -Wall vlidort_pack.f90
gfortran -c -Wall vlidort_unpack.f90
gfortran -c -Wall vlidort_vfo_interface.f90

gfortran -c -Wall vlidort_masters_V2p8p3.f90

exit

