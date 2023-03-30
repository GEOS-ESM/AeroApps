esma_set_this ()
find_package(F2PY3 REQUIRED)

# MAT: Note uses add_f2py3 because I can't figure out the import test for netCDF f2py code
add_f2py3_module (VLIDORT_STN_
    SOURCES VLIDORT_STN_py.F90
  DESTINATION lib/Python/${this}
  LIBRARIES geo_vlidort
  INCLUDEDIRS ${CMAKE_CURRENT_BINARY_DIR} ${CMAKE_BINARY_DIR}/lib ${include_${this}} ${esma_include}/geo_vlidort
  )
add_dependencies(VLIDORT_STN_ geo_vlidort)

esma_add_f2py3_module (stnAngles_
    SOURCES stnAngles_py.F90
  DESTINATION lib/Python/${this}
  LIBRARIES
  INCLUDEDIRS ${CMAKE_CURRENT_BINARY_DIR} ${CMAKE_BINARY_DIR}/lib ${include_${this}}
  )

file (GLOB python_files pystn/*.py)
install (FILES ${python_files} DESTINATION lib/Python/pystn)

# these are executable scripts
set (PYSCRIPTS
    aeronet_vlidort.py
    directsun_vlidort.py
    ext_stn_vlidort_sampler.py
    g5nr_stn_sampler.py
    mcd12c_stn_sampler.py
    mcd43c_stn_sampler.py
    run_g5nr_stn_sampler_invariant.py
    run_g5nr_stn_sampler.py
    run_stn_vlidort.py
    stn_vlidort_aux.py 
    )

install (
    PROGRAMS  ${PYSCRIPTS}
    DESTINATION bin/${this}
    )

set (RCFILES
    Aod_EOS.rc
    Chem_MieRegistry.rc
    )

install (
    FILES ${RCFILES} 
    DESTINATION bin/${this}
    )
