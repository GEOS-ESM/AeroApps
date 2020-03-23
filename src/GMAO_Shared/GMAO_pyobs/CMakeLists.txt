esma_set_this ()

set (SRCS
  VegType_Mod.F90
  VegType_io.c
  )

esma_add_library (${this}
  SRCS ${SRCS}
  DEPENDENCIES GMAO_gfio_r4
  INCLUDES ${INC_ESMF})

include(UseF2Py)

add_f2py_module(binObs_
   SOURCES binObs_py.F
   DESTINATION lib/Python/pyobs
   INCLUDEDIRS ${CMAKE_CURRENT_BINARY_DIR} ${CMAKE_BINARY_DIR}/lib ${include_${this}}
   )
add_dependencies(binObs_ ${this})

add_f2py_module(IGBP_
   SOURCES IGBP_py.F90
   DESTINATION lib/Python/pyobs
   ONLY getsimpleveg getdetailedveg
   LIBRARIES GMAO_pyobs
   INCLUDEDIRS ${CMAKE_CURRENT_BINARY_DIR} ${CMAKE_BINARY_DIR}/lib ${include_${this}}
   )
add_dependencies(IGBP_ ${this})

file (GLOB python_files pyobs/*.py)
install (FILES ${python_files} DESTINATION lib/Python/pyobs)
