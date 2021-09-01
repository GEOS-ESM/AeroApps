#rm *.o *.mod
#exit

gfortran -c -Wall ../../vlidort_def/vlidort_pars.f90
gfortran -c -Wall vsleave_findpar.f90
gfortran -c -Wall vsleave_sup_inputs_def.f90
gfortran -c -Wall vsleave_sup_outputs_def.f90
gfortran -c -Wall vsleave_sup_aux.f90
gfortran -c -Wall vsleave_sup_routines.f90
gfortran -c -Wall vsleave_sup_masters.f90
gfortran -c -Wall vsleave_sup_mod.f90
#exit

gfortran -c -Wall vsleave_lin_sup_inputs_def.f90
gfortran -c -Wall vsleave_lin_sup_outputs_def.f90
gfortran -c -Wall vsleave_lin_sup_routines.f90

gfortran -c -Wall vsleave_lin_sup_masters.f90
gfortran -c -Wall vsleave_lin_sup_mod.f90
exit
