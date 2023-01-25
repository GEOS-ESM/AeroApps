rm *.o *.mod
exit

gfortran -c -Wall vlidort_pars.f90
gfortran -c -Wall vlidort_inputs_def.f90
gfortran -c -Wall vlidort_outputs_def.f90
gfortran -c -Wall vlidort_greens_def.f90
gfortran -c -Wall vlidort_sup_brdf_def.f90
gfortran -c -Wall vlidort_sup_sleave_def.f90
gfortran -c -Wall vlidort_sup_ss_def.f90
gfortran -c -Wall vlidort_work_def.f90
gfortran -c -Wall vlidort_sup_def.f90
gfortran -c -Wall vlidort_io_defs.f90

gfortran -c -Wall vlidort_lin_greens_def.f90
gfortran -c -Wall vlidort_lin_inputs_def.f90
gfortran -c -Wall vlidort_lin_outputs_def.f90
gfortran -c -Wall vlidort_lin_sup_brdf_def.f90
gfortran -c -Wall vlidort_lin_sup_sleave_def.f90
gfortran -c -Wall vlidort_lin_sup_ss_def.f90
gfortran -c -Wall vlidort_lin_work_def.f90
gfortran -c -Wall vlidort_lin_sup_def.f90
gfortran -c -Wall vlidort_lin_io_defs.f90

exit
