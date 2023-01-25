esma_set_this ()

set (src_dirs vlidort_main/regular vlidort_main/linearized vlidort_test)
foreach (dir ${src_dirs})
    file (GLOB srcs_ CONFIGURE_DEPENDS ${dir}/*.[Ffc] ${dir}/*.[fF]90)
  list (APPEND srcs ${srcs_})
endforeach()

set (DEF_PATH vlidort_def)
list (APPEND srcs
  ${DEF_PATH}/rayleigh.f90
  ${DEF_PATH}/rd_fp_single.f90
  ${DEF_PATH}/rd_fp.f90
  ${DEF_PATH}/interp.f90
  ${DEF_PATH}/gravity.f90
  ${DEF_PATH}/vlidort_sup_brdf_def.f90     
  ${DEF_PATH}/vlidort_sup_ss_def.f90       
  ${DEF_PATH}/vlidort_sup_sleave_def.f90   
  ${DEF_PATH}/vlidort_sup_def.f90          
  ${DEF_PATH}/vlidort_outputs_def.f90      
  ${DEF_PATH}/vlidort_inputs_def.f90       
  ${DEF_PATH}/vlidort_io_defs.f90          
  ${DEF_PATH}/vlidort_work_def.f90         
  ${DEF_PATH}/create_atm_prfl.f90          
  ${DEF_PATH}/create_atm_prfl_single.f90
  ${DEF_PATH}/vlidort_pars.f90
  )
   
set (FO_PATH vlidort_focode)
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
   vlidort_def/vlidort_pars.f90
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

list (APPEND srcs vlidort_test_gas_aerosol_retrieval_single.f90)

ecbuild_add_executable (
  TARGET vlidort_test_gas_aerosol_retrieval_single.x
  SOURCES ${srcs}
  LIBS GMAO_mpeu GMAO_gfio_r4 MAPL_cfio_r4 Chem_Base MAPL hdf5hl_fortran hdf5_fortran ${NETCDF_LIBRARIES} MPI::MPI_Fortran
  INCLUDES ${INC_HDF5}
  NOINSTALL
  )

install (
    TARGETS vlidort_test_gas_aerosol_retrieval_single.x
    DESTINATION bin/missions/${this}
    )

#esma_add_library (${this}
#  SRCS ${srcs}
#  )
