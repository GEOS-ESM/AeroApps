esma_set_this ()

foreach (base Mie_ MieObs_)

  esma_add_f2py_module (${base}
    SOURCES ${base}py.F90
    DESTINATION lib/Python/${this}
    LIBRARIES Chem_Base
    INCLUDEDIRS ${CMAKE_CURRENT_BINARY_DIR} ${CMAKE_BINARY_DIR}/lib
    ${include_${this}} ${esma_include}/Chem_Base
    )

endforeach()

esma_add_f2py_module (OMI_
  SOURCES OMI_py.F90
  DESTINATION lib/Python/${this}
  LIBRARIES Chem_Base VLIDORT
  INCLUDEDIRS ${CMAKE_CURRENT_BINARY_DIR} ${CMAKE_BINARY_DIR}/lib
  ${include_${this}} ${esma_include}/VLIDORT
  )
add_dependencies(OMI_ VLIDORT)

#esma_add_f2py_module (VLIDORT2_BRDF_
#  SOURCES VLIDORT_BRDF_py.F90
#  DESTINATION lib/Python/${this}
#  INCLUDEDIRS ${CMAKE_CURRENT_BINARY_DIR} ${CMAKE_BINARY_DIR}/lib ${include_${this}}
#  )

#esma_add_f2py_module (VLIDORT_OMI_
#  SOURCES VLIDORT_OMI_py.F90
#  DESTINATION lib/Python/${this}
#  LIBRARIES VLIDORT90
#  INCLUDEDIRS ${CMAKE_CURRENT_BINARY_DIR} ${CMAKE_BINARY_DIR}/lib
#  ${include_${this}}  ${esma_include}/VLIDORT90)

#esma_add_f2py_module (VLIDORT_cloud_
#  SOURCES VLIDORT_cloud_py.F90
#  DESTINATION lib/Python/${this}
#  INCLUDEDIRS ${CMAKE_CURRENT_BINARY_DIR} ${CMAKE_BINARY_DIR}/lib ${include_${this}}
#  )

esma_add_f2py_module (qsat_
  SOURCES qsat_py.F90
  DESTINATION lib/Python/${this}
  INCLUDEDIRS ${CMAKE_CURRENT_BINARY_DIR} ${CMAKE_BINARY_DIR}/lib ${include_${this}}
  )

esma_add_f2py_module (binObs_
  SOURCES binObs_py.F
  DESTINATION lib/Python/${this}
  INCLUDEDIRS ${CMAKE_CURRENT_BINARY_DIR} ${CMAKE_BINARY_DIR}/lib ${include_${this}}
  )

esma_add_f2py_module (test_mieobs
  SOURCES test_miobs.F90
  DESTINATION lib/Python/${this}
  INCLUDEDIRS ${CMAKE_CURRENT_BINARY_DIR} ${CMAKE_BINARY_DIR}/lib ${include_${this}}
  )

