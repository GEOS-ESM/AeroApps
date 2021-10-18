esma_add_subdirectories (Expansion_Coefs)
esma_set_this ()
find_package(F2PY2 REQUIRED)

ecbuild_add_executable (
  TARGET leo_vlidort_cloud.x
  SOURCES leo_vlidort_cloud.F90 
  LIBS GMAO_mpeu GMAO_gfio_r4 MAPL_cfio_r4 Chem_Base MAPL VLIDORT2p8p3 geo_vlidort ${NETCDF_LIBRARIES} MPI::MPI_Fortran
  )

# MAT: Note uses add_f2py2 because I can't figure out the import test for netCDF f2py code
add_f2py2_module (VLIDORT_POLAR_
  SOURCES VLIDORT_POLAR_py.F90
  DESTINATION lib/Python/${this}
  LIBRARIES geo_vlidort VLIDORT2p8p3
  INCLUDEDIRS ${CMAKE_CURRENT_BINARY_DIR} ${CMAKE_BINARY_DIR}/lib ${include_${this}} ${esma_include}/geo_vlidort ${esma_include}/VLIDORT2p8p3
  )
add_dependencies(VLIDORT_POLAR_ geo_vlidort VLIDORT2p8p3)

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


file (GLOB python_files CONFIGURE_DEPENDS py_leo_vlidort/*.py)
install (FILES ${python_files} DESTINATION lib/Python/py_leo_vlidort)


# these are executable scripts
set (PYSCRIPTS
    batch_lidar_sampler.py
    benchmark_polar_vlidort_grasp.py
    benchmark_polar_vlidort_hyperspectral.py
    benchmark_polar_vlidort_iprt.py
    benchmark_polar_vlidort.py
    brdf_mcd43c_sampler.py
    landcover_mcd12c_sampler.py
    leo_vlidort_lc.py
    lidar_sampler.py
    mp_lidar_vlidort.py
    pace_test.py
    plot_diff_iprt.py
    plot_diff.py
    plot_sleave.py
    polar_aero_dist.py
    polar_vlidort.py
    run_ext_sampler.py
    run_lidar_sampler.py
    run_lidar_vlidort.py
    run_polar_vlidort.py
    trace_gas_test.py
    )

set (SCRIPTS
    run_lidar_sampler.j
    )

install (
    PROGRAMS  ${PYSCRIPTS} ${SCRIPTS}
    DESTINATION bin/${this}
    )

# these are non-executaables
set (PYFILES
    mcd12c.py
    )

set (RCFILES
    aer_Nv_opendap.rc
    aer_Nv.rc
    Aod_EOS.rc
    asm_Nx.rc
    Chem_MieRegistry.rc
    chm_Nv_opendap.rc
    chm_Nv.rc
    met_Nv_opendap.rc
    met_Nv.rc
    myd13c2.rc
    omi_ler.rc
    )

set (OTHERFILES
    polar_vlidort.pcf
    )
install (
    FILES ${PYFILES} ${RCFILES} ${OTHERFILES}
    DESTINATION bin/${this}
    )

install (
   DIRECTORY calipso
   DESTINATION bin/${this}
   USE_SOURCE_PERMISSIONS
   MESSAGE_NEVER
   )
