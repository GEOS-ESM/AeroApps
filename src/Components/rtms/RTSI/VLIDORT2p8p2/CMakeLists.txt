esma_set_this ()

set (src_dirs DEF_directory fo_main sourcecode_str/regular_modules sourcecode_str/linearized_modules LAPACK_directory)
set (srcs VLIDORT_Mod.F90 VLIDORT_SurfaceMod.F90 VLIDORT_ScatMod.F90)
foreach (dir ${src_dirs})
    file (GLOB srcs_ CONFIGURE_DEPENDS ${dir}/*.[Ffc] ${dir}/*.[fF]90)
  list (APPEND srcs ${srcs_})
endforeach()

set(VACCESS_PATH vsup/vaccessories)
set(VBRDF_PATH vsup/vbrdf)
list (APPEND srcs
   ${VACCESS_PATH}/vlidort_vbrdf_sup_accessories.f90
   DEF_directory/vlidort_pars.f90
   ${VBRDF_PATH}/vbrdf_findpar.f90
   ${VBRDF_PATH}/vbrdf_sup_inputs_def.f90
   ${VBRDF_PATH}/vbrdf_sup_outputs_def.f90
   ${VBRDF_PATH}/vbrdf_sup_aux.f90
   ${VBRDF_PATH}/vbrdf_sup_kernels.f90
   ${VBRDF_PATH}/vbrdf_sup_routines.f90
   ${VBRDF_PATH}/vbrdf_sup_masters.f90
   ${VBRDF_PATH}/vbrdf_sup_mod.f90
   ${VBRDF_PATH}/vbrdf_lin_sup_inputs_def.f90
   ${VBRDF_PATH}/vbrdf_lin_sup_outputs_def.f90
   ${VBRDF_PATH}/vbrdf_lin_sup_kernels.f90
   ${VBRDF_PATH}/vbrdf_lin_sup_routines.f90
   ${VBRDF_PATH}/vbrdf_lin_sup_masters.f90
   ${VBRDF_PATH}/vbrdf_lin_sup_mod.f90
   ${VACCESS_PATH}/vlidort_vbrdf_lin_sup_accessories.f90
)

set(VSLEAVE_PATH vsup/vsleave)
list (APPEND srcs
   ${VSLEAVE_PATH}/vsleave_findpar.f90
   ${VSLEAVE_PATH}/vsleave_sup_inputs_def.f90
   ${VSLEAVE_PATH}/vsleave_sup_outputs_def.f90
   ${VSLEAVE_PATH}/vsleave_sup_aux.f90
   ${VSLEAVE_PATH}/vsleave_sup_routines.f90
   ${VSLEAVE_PATH}/vsleave_sup_masters.f90
   ${VSLEAVE_PATH}/vsleave_sup_mod.f90
   ${VACCESS_PATH}/vlidort_vsleave_sup_accessories.f90
   ${VSLEAVE_PATH}/vsleave_lin_sup_inputs_def.f90
   ${VSLEAVE_PATH}/vsleave_lin_sup_outputs_def.f90
   ${VSLEAVE_PATH}/vsleave_lin_sup_routines.f90
   ${VSLEAVE_PATH}/vsleave_lin_sup_masters.f90
   ${VSLEAVE_PATH}/vsleave_lin_sup_mod.f90
   ${VACCESS_PATH}/vlidort_vsleave_lin_sup_accessories.f90
)  

set(VFZMAT_PATH vsup/vfzmat)
list (APPEND srcs
    ${VFZMAT_PATH}/vfzmat_Numerical.f90
    ${VFZMAT_PATH}/vfzmat_DevelopCoeffs.f90
    ${VFZMAT_PATH}/vfzmat_ExpandCoeffs.f90
    ${VFZMAT_PATH}/vfzmat_Rotation.f90
    ${VFZMAT_PATH}/vfzmat_PhasMat.f90
    ${VFZMAT_PATH}/vfzmat_Rayleigh.f90
    ${VFZMAT_PATH}/vfzmat_Master.f90
)

esma_add_library (${this}
  SRCS ${srcs}
  )
