esma_set_this ()

set (src_dirs DEF_directory sourcecode_str/regular_modules sourcecode_str/linearized_modules LAPACK_directory)
set (srcs VLIDORT_Mod.F90 VLIDORT_SurfaceMod.F90 VLIDORT_ScatMod.F90)
foreach (dir ${src_dirs})
    file (GLOB srcs_ CONFIGURE_DEPENDS ${dir}/*.[Ffc] ${dir}/*.[fF]90)
  list (APPEND srcs ${srcs_})
endforeach()

set (FO_PATH fo_main)
list (APPEND srcs
  ${FO_PATH}/FO_geometry_Generic.f90  
  ${FO_PATH}/FO_WPgeometry_Routines.f90
  ${FO_PATH}/FO_DTWPgeometry_master.f90
  ${FO_PATH}/FO_SSWPgeometry_master.f90
  ${FO_PATH}/FO_Planckfunc.f90     
  ${FO_PATH}/FO_Taylor.f90
  ${FO_PATH}/FO_Thermal_RTCalcs_I.f90
  ${FO_PATH}/FO_Thermal_RTCalcs_ILCS.f90
  ${FO_PATH}/FO_Thermal_RTCalcs_ILPS.f90
  ${FO_PATH}/FO_VectorSS_Spherfuncs.f90
  ${FO_PATH}/FO_VectorSS_RTCalcs_I.f90
  ${FO_PATH}/VFO_Master.f90 
  ${FO_PATH}/FO_VectorSS_RTCalcs_ILCS.f90
  ${FO_PATH}/FO_VectorSS_RTCalcs_ILPS.f90
  ${FO_PATH}/VFO_LinMasters.f90
  )

set(VBRDF_PATH vsup/vbrdf)
list (APPEND srcs
   DEF_directory/vlidort_pars.f90
   ${VBRDF_PATH}/vbrdf_findpar.f90
   ${VBRDF_PATH}/vbrdf_sup_inputs_def.f90
   ${VBRDF_PATH}/vbrdf_sup_outputs_def.f90
   ${VBRDF_PATH}/vbrdf_sup_aux.f90
   ${VBRDF_PATH}/vbrdf_sup_kernels.f90
   ${VBRDF_PATH}/vbrdf_sup_routines.f90
   ${VBRDF_PATH}/vbrdf_sup_masters.f90
   ${VBRDF_PATH}/vbrdf_sup_mod.f90
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
)  

esma_add_library (${this}
  SRCS ${srcs}
  )
