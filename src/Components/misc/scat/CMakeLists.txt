esma_set_this ()

add_f2py_module (Mie_
  SOURCES Mie_py.F90
  DESTINATION lib/Python/${this}
  INCLUDEDIRS ${CMAKE_CURRENT_BINARY_DIR} ${CMAKE_BINARY_DIR}/lib ${include_${this}}
  )

add_f2py_module (MieObs2_
  SOURCES MieObs_py.F90
  DESTINATION lib/Python/${this}
  INCLUDEDIRS ${CMAKE_CURRENT_BINARY_DIR} ${CMAKE_BINARY_DIR}/lib ${include_${this}}
  )

add_f2py_module (OMI_
  SOURCES OMI_py.F90
  DESTINATION lib/Python/${this}
  INCLUDEDIRS ${CMAKE_CURRENT_BINARY_DIR} ${CMAKE_BINARY_DIR}/lib ${include_${this}}
  )

add_f2py_module (VLIDORT2_BRDF_
  SOURCES VLIDORT_BRDF_py.F90
  DESTINATION lib/Python/${this}
  INCLUDEDIRS ${CMAKE_CURRENT_BINARY_DIR} ${CMAKE_BINARY_DIR}/lib ${include_${this}}
  )

add_f2py_module (VLIDORT_OMI_
  SOURCES VLIDORT_OMI_py.F90
  DESTINATION lib/Python/${this}
  INCLUDEDIRS ${CMAKE_CURRENT_BINARY_DIR} ${CMAKE_BINARY_DIR}/lib ${include_${this}}
  )

add_f2py_module (VLIDORT_cloud_
  SOURCES VLIDORT_cloud_py.F90
  DESTINATION lib/Python/${this}
  INCLUDEDIRS ${CMAKE_CURRENT_BINARY_DIR} ${CMAKE_BINARY_DIR}/lib ${include_${this}}
  )

add_f2py_module (qsat_
  SOURCES qsat_py.F90
  DESTINATION lib/Python/${this}
  INCLUDEDIRS ${CMAKE_CURRENT_BINARY_DIR} ${CMAKE_BINARY_DIR}/lib ${include_${this}}
  )

add_f2py_module (binObs2_
  SOURCES binObs_py.F90
  DESTINATION lib/Python/${this}
  INCLUDEDIRS ${CMAKE_CURRENT_BINARY_DIR} ${CMAKE_BINARY_DIR}/lib ${include_${this}}
  )

add_f2py_module (test_mieobs
  SOURCES test_miobs.F90
  DESTINATION lib/Python/${this}
  INCLUDEDIRS ${CMAKE_CURRENT_BINARY_DIR} ${CMAKE_BINARY_DIR}/lib ${include_${this}}
  )

