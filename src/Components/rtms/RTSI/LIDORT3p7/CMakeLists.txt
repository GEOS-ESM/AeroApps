esma_set_this ()

set (src_dirs DEF_directory sourcecode_str LAPACK_directory)
set (srcs LIDORT_Mod.F90 LIDORT_SurfaceMod.F90 LIDORT_ScatMod.F90)
foreach (dir ${src_dirs})
    file (GLOB srcs_ CONFIGURE_DEPENDS ${dir}/*.[Ffc] ${dir}/*.[fF]90)
  list (APPEND srcs ${srcs_})
endforeach()


set(BRDF_PATH sup/brdf)
list (APPEND srcs
   DEF_directory/lidort_pars.f90
   ${BRDF_PATH}/brdf_findpar.f90
   ${BRDF_PATH}/brdf_sup_inputs_def.f90
   ${BRDF_PATH}/brdf_sup_outputs_def.f90
   ${BRDF_PATH}/brdf_sup_aux.f90
   ${BRDF_PATH}/brdf_sup_kernels.f90
   ${BRDF_PATH}/brdf_sup_routines.f90
   ${BRDF_PATH}/brdf_sup_masters.f90
   ${BRDF_PATH}/brdf_sup_mod.f90
)

set(SLEAVE_PATH sup/sleave)
list (APPEND srcs
   ${SLEAVE_PATH}/sleave_findpar.f90
   ${SLEAVE_PATH}/sleave_sup_inputs_def.f90
   ${SLEAVE_PATH}/sleave_sup_outputs_def.f90
   ${SLEAVE_PATH}/sleave_sup_aux.f90
   ${SLEAVE_PATH}/sleave_sup_routines.f90
   ${SLEAVE_PATH}/sleave_sup_masters.f90
   ${SLEAVE_PATH}/sleave_sup_mod.f90
)  

esma_add_library (${this}
  SRCS ${srcs}
  )
