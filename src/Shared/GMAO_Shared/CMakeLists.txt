esma_check_if_debug()

set (subdirectories   
  GMAO_mpeu    GMAO_etc   GEOS_Shared
  GMAO_hermes  GMAO_ods
  GMAO_psas    GMAO_pyobs   GMAO_aeropyobs
  Chem_Base    GMAO_pilgrim MAPL
  )
foreach (subdir ${subdirectories})
  esma_add_subdirectory (${subdir})
endforeach()
    

# Special case - GMAO_gfio is built twice with two different precisions.
if (EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/GMAO_gfio)
  add_subdirectory (GMAO_gfio GMAO_gfio_r4)
  add_subdirectory (GMAO_gfio GMAO_gfio_r8)
  add_dependencies (GMAO_gfio_r4 GMAO_gfio_r8)
endif ()

