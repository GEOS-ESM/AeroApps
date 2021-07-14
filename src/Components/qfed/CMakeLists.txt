esma_set_this()

set(srcs  PlumeRise_ConstMod.F90 PlumeRise_Mod.F90 LockPlume_Mod.F90 qsat_Mod.F90)

esma_add_library (${this}
  SRCS ${srcs}
  DEPENDENCIES GMAO_mpeu GMAO_hermes GMAO_pyobs # crtm?
  )

find_package(F2PY2 REQUIRED)

esma_add_f2py2_module (ut_
  SOURCES ut_py.F90
  DESTINATION lib/Python/${this}
  INCLUDEDIRS ${CMAKE_CURRENT_BINARY_DIR} ${CMAKE_BINARY_DIR}/lib ${include_${this}}
  )

esma_add_f2py2_module (PlumeRise_
  SOURCES qfed/PlumeRise_py.F90
  DESTINATION lib/Python/${this}
  INCLUDEDIRS ${CMAKE_CURRENT_BINARY_DIR} ${CMAKE_BINARY_DIR}/lib ${include_${this}}
  )

esma_add_f2py2_module (LockPlume_
  SOURCES qfed/LockPlume_py.F90
  DESTINATION lib/Python/${this}
  INCLUDEDIRS ${CMAKE_CURRENT_BINARY_DIR} ${CMAKE_BINARY_DIR}/lib ${include_${this}}
  )

esma_add_f2py2_module (tausfc_
  SOURCES qfed/tausfc_py.F90
  DESTINATION lib/Python/${this}
  INCLUDEDIRS ${CMAKE_CURRENT_BINARY_DIR} ${CMAKE_BINARY_DIR}/lib ${include_${this}}
  )

esma_add_f2py2_module (crtmmodis_
  SOURCES qfed/crtmmodis_py.F90
  DESTINATION lib/Python/${this}
  INCLUDEDIRS ${CMAKE_CURRENT_BINARY_DIR} ${CMAKE_BINARY_DIR}/lib ${include_${this}}
  )

