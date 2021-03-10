esma_set_this()

set(srcs  PlumeRise_ConstMod.F90 PlumeRise_Mod.F90 LockPlume_Mod.F90 qsat_Mod.F90)

esma_add_library (${this}
  SRCS ${srcs}
  DEPENDENCIES GMAO_mpeu GMAO_hermes GMAO_pyobs # crtm?
  )

esma_add_f2py_module (ut_
  SOURCES ut_py.F90
  DESTINATION lib/Python/${this}
  INCLUDEDIRS ${CMAKE_CURRENT_BINARY_DIR} ${CMAKE_BINARY_DIR}/lib ${include_${this}}
  )

esma_add_f2py_module (PlumeRise_
  SOURCES qfed/PlumeRise_py.F90
  DESTINATION lib/Python/${this}
  INCLUDEDIRS ${CMAKE_CURRENT_BINARY_DIR} ${CMAKE_BINARY_DIR}/lib ${include_${this}}
  )

esma_add_f2py_module (LockPlume_
  SOURCES qfed/LockPlume_py.F90
  DESTINATION lib/Python/${this}
  INCLUDEDIRS ${CMAKE_CURRENT_BINARY_DIR} ${CMAKE_BINARY_DIR}/lib ${include_${this}}
  )

esma_add_f2py_module (tausfc_
  SOURCES qfed/tausfc_py.F90
  DESTINATION lib/Python/${this}
  INCLUDEDIRS ${CMAKE_CURRENT_BINARY_DIR} ${CMAKE_BINARY_DIR}/lib ${include_${this}}
  )

esma_add_f2py_module (crtmmodis_
  SOURCES qfed/crtmmodis_py.F90
  DESTINATION lib/Python/${this}
  INCLUDEDIRS ${CMAKE_CURRENT_BINARY_DIR} ${CMAKE_BINARY_DIR}/lib ${include_${this}}
  )

