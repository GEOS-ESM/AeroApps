esma_set_this ()

find_package(F2PY2 REQUIRED)

# MAT: Note uses add_f2py2 because I can't figure out the import test for netCDF f2py code
add_f2py2_module (Mie_
  SOURCES Mie_py.F90
  DESTINATION lib/Python/${this}
  LIBRARIES Chem_Base GMAO_mpeu MAPL GMAO_gfio_r4 ${NETCDF_LIBRARIES}
  INCLUDEDIRS ${CMAKE_CURRENT_BINARY_DIR} ${CMAKE_BINARY_DIR}/lib ${include_${this}} ${esma_include}/Chem_Base ${BASEDIR}/lib ${INC_NETCDF}
  USE_MPI
  )
add_dependencies(Mie_ Chem_Base GMAO_mpeu MAPL GMAO_gfio_r4)

# MAT: Note uses add_f2py2 because I can't figure out the import test for netCDF f2py code
add_f2py2_module (scat_MieObs_
  SOURCES scat_MieObs_py.F90
  DESTINATION lib/Python/${this}
  LIBRARIES Chem_Base GMAO_mpeu MAPL ${NETCDF_LIBRARIES}
  INCLUDEDIRS ${CMAKE_CURRENT_BINARY_DIR} ${CMAKE_BINARY_DIR}/lib ${include_${this}} ${esma_include}/Chem_Base ${BASEDIR}/lib ${INC_NETCDF}
  )
add_dependencies(scat_MieObs_ Chem_Base GMAO_mpeu MAPL)

esma_add_f2py2_module (OMI_
  SOURCES OMI_py.F90
  DESTINATION lib/Python/${this}
  LIBRARIES Chem_Base VLIDORT
  INCLUDEDIRS ${CMAKE_CURRENT_BINARY_DIR} ${CMAKE_BINARY_DIR}/lib ${include_${this}} ${esma_include}/Chem_Base ${esma_include}/VLIDORT
  )
add_dependencies(OMI_ Chem_Base VLIDORT)

esma_add_f2py2_module (VLIDORT_BRDF_
  SOURCES VLIDORT_BRDF_py.F90
  DESTINATION lib/Python/${this}
  LIBRARIES VLIDORT90
  INCLUDEDIRS ${CMAKE_CURRENT_BINARY_DIR} ${CMAKE_BINARY_DIR}/lib ${include_${this}} ${esma_include}/VLIDORT90)
add_dependencies(VLIDORT_BRDF_ VLIDORT90)

esma_add_f2py2_module (VLIDORT_OMI_
  SOURCES VLIDORT_OMI_py.F90
  DESTINATION lib/Python/${this}
  LIBRARIES VLIDORT90
  INCLUDEDIRS ${CMAKE_CURRENT_BINARY_DIR} ${CMAKE_BINARY_DIR}/lib ${include_${this}}  ${esma_include}/VLIDORT90)
add_dependencies(VLIDORT_OMI_ VLIDORT90)

esma_add_f2py2_module (VLIDORT_cloud_
  SOURCES VLIDORT_cloud_py.F90
  DESTINATION lib/Python/${this}
  LIBRARIES Chem_Base VLIDORT90
  INCLUDEDIRS ${CMAKE_CURRENT_BINARY_DIR} ${CMAKE_BINARY_DIR}/lib ${include_${this}} ${esma_include}/Chem_Base ${esma_include}/VLIDORT90
  )
add_dependencies(VLIDORT_cloud_ Chem_Base VLIDORT90)

# Using add_f2py2_module here as I can't get the built-in import test to work
esma_add_f2py2_module (qsat_
  SOURCES qsat_py.F90
  DESTINATION lib/Python/${this}
  LIBRARIES MAPL.base
  INCLUDEDIRS ${CMAKE_CURRENT_BINARY_DIR} ${CMAKE_BINARY_DIR}/lib ${include_${this}} ${esma_include}/MAPL.base
  )
add_dependencies(qsat_ MAPL.base)

esma_add_f2py2_module (scat_binObs_
  SOURCES scat_binObs_py.F
  DESTINATION lib/Python/${this}
  INCLUDEDIRS ${CMAKE_CURRENT_BINARY_DIR} ${CMAKE_BINARY_DIR}/lib ${include_${this}}
  )

# these are executable scripts
set (PYSCRIPTS 
    omi_l3a.py
    check_ods.py
    sample_omaeruv.py
    )

set (SCRIPTS
    killme.csh
    )

install (
    PROGRAMS  ${PYSCRIPTS} ${SCRIPTS}
    DESTINATION bin/${this}
    )

# these are non-executaables
set (PYFILES
    calipso_simul.py
    CATS_simul.py
    mieobs.py
    mie.py
    mxd04_simul.py
    mxd04_simul_vector.py
    omi.py
    omi_l2a_simul.py
    omi_simul.py
    pace_simul.py
    test_gfio.py
    test_pmom.py
    vlidort.py
    )

set (RCFILES
    Aod_EOS.rc
    Aod_MODIS.rc
    Aod_PACE.rc
    Chem_MieRegistry.rc
    )

set (OTHERFILES
    omi.j
    )

install (
    FILES ${PYFILES} ${RCFILES} ${OTHERFILES}
    DESTINATION bin/${this}
    )
