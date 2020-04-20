esma_set_this ()

add_f2py_module (VLIDORT_BRDF_
  SOURCES VLIDORT_BRDF_py.F90
  DESTINATION lib/Python/${this}
  INCLUDEDIRS ${CMAKE_CURRENT_BINARY_DIR} ${CMAKE_BINARY_DIR}/lib ${include_${this}}
  )
  
