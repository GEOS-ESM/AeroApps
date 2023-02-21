esma_set_this ()

set (SRCS
  VegType_Mod.F90
  VegType_io.c
  glint_mod.F90  
  sgp4_mod.F90  
  TLE_mod.F90
  cubed_sphere_utilities.F90
  )

esma_add_library (${this}
  SRCS ${SRCS}
  DEPENDENCIES GMAO_gfio_r4
  INCLUDES ${INC_ESMF})

find_package(F2PY2 REQUIRED)
esma_add_f2py2_module(binObs_
   SOURCES binObs_py.F
   DESTINATION lib/Python
   INCLUDEDIRS ${CMAKE_CURRENT_BINARY_DIR} ${CMAKE_BINARY_DIR}/lib ${include_${this}}
   )
add_dependencies(binObs_ ${this})

esma_add_f2py2_module(IGBP_
   SOURCES IGBP_py.F90
   DESTINATION lib/Python/pyobs
   ONLY getsimpleveg getdetailedveg
   LIBRARIES GMAO_pyobs
   INCLUDEDIRS ${CMAKE_CURRENT_BINARY_DIR} ${CMAKE_BINARY_DIR}/lib ${include_${this}}
   )
add_dependencies(IGBP_ ${this})

esma_add_f2py2_module(sgp4_
   SOURCES sgp4_py.F90
   DESTINATION lib/Python/pyobs
   LIBRARIES GMAO_pyobs
   INCLUDEDIRS ${CMAKE_CURRENT_BINARY_DIR} ${CMAKE_BINARY_DIR}/lib ${include_${this}}
   )
add_dependencies(sgp4_ ${this})

file (GLOB python_files pyobs/*.py)
install (FILES ${python_files} DESTINATION lib/Python/pyobs)

install(PROGRAMS active_aeronet.py  ext_sampler.py  stn_sampler.py  stn_xsect.py  trj_sampler.py  trj_xsect.py DESTINATION bin)

