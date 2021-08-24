esma_add_subdirectories (Expansion_Coefs)
esma_set_this ()
find_package(F2PY2 REQUIRED)

ecbuild_add_executable (
  TARGET leo_vlidort_cloud.x
  SOURCES leo_vlidort_cloud.F90 
  LIBS GMAO_mpeu GMAO_gfio_r4 MAPL_cfio_r4 Chem_Base MAPL VLIDORT2OS geo_vlidort ${NETCDF_LIBRARIES} MPI::MPI_Fortran
  )

# MAT: Note uses add_f2py2 because I can't figure out the import test for netCDF f2py code
add_f2py2_module (VLIDORT_POLAR_
  SOURCES VLIDORT_POLAR_py.F90
  DESTINATION lib/Python/${this}
  LIBRARIES geo_vlidort
  INCLUDEDIRS ${CMAKE_CURRENT_BINARY_DIR} ${CMAKE_BINARY_DIR}/lib ${include_${this}} ${esma_include}/geo_vlidort
  )
add_dependencies(VLIDORT_POLAR_ geo_vlidort)

# MAT: Note uses add_f2py2 because I can't figure out the import test for netCDF f2py code
add_f2py2_module (READ_NC_
  SOURCES READ_NC_py.F90
  DESTINATION lib/Python/${this}
  LIBRARIES geo_vlidort 
  INCLUDEDIRS ${CMAKE_CURRENT_BINARY_DIR} ${CMAKE_BINARY_DIR}/lib ${include_${this}} ${esma_include}/geo_vlidort 
  )
add_dependencies(READ_NC_ geo_vlidort)

esma_add_f2py2_module (LidarAngles_
    SOURCES LidarAngles_py.F90
  DESTINATION lib/Python/${this}
  LIBRARIES 
  INCLUDEDIRS ${CMAKE_CURRENT_BINARY_DIR} ${CMAKE_BINARY_DIR}/lib ${include_${this}}
  )


# these are executable scripts
#set (PYSCRIPTS
#    geo_vlidort_AI_lc2.py
#    geo_vlidort_blksurf_lc2.py
#    geo_vlidort_cloud_lc2.py
#    geo_vlidort_lc2.py
#    geo_vlidort_ocean_lc2.py
#    geo_warp.py
#    run_geo_vlidort_AI_lc2.py
#    run_geo_vlidort_blksurf_lc2.py
#    run_geo_vlidort_cloud_lc2.py
#    run_geo_vlidort_lc2.py
#    run_geo_vlidort_ocean_lc2.py
#    utAngles.py
#    )

#install (
#    PROGRAMS  ${PYSCRIPTS}
#    DESTINATION bin/${this}
#    )


