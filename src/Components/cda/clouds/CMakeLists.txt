esma_set_this ()

set (srcs
  mcsRegrid_Mod.F90
  mod_mkl.F90
  mod_StaticGridcolumnStatisticalModel.F90
  mod_GaussianCopula.F90
  mod_consts.F90
  mod_eswi.F90
  mod_icefrac.F90
  mod_reff.F90
  mod_simTauRe.F90
  mod_triangle.F90
  mod_utils.F90
  )

esma_add_library (${this}
  SRCS ${srcs}
  DEPENDENCIES Chem_Base MAPL GMAO_mpeu GMAO_ods gfio_r4
  PUBLIC_LIBS ${MKL_LIBRARIES}
  PUBLIC_INCLUDES ${MKL_INCLUDE_DIRS}
  )



add_f2py_module (Aero_
  SOURCES Aero_py.F90
  DESTINATION lib/Python/${this}
  INCLUDEDIRS ${CMAKE_CURRENT_BINARY_DIR} ${CMAKE_BINARY_DIR}/lib ${include_${this}}
  )

add_f2py_module (ICA_
  SOURCES ICA_py.F90
  DESTINATION lib/Python/${this}
  INCLUDEDIRS ${CMAKE_CURRENT_BINARY_DIR} ${CMAKE_BINARY_DIR}/lib ${include_${this}}
  )

add_f2py_module (mcsRegrid_
  SOURCES mcsRegrid_py.F90
  DESTINATION lib/Python/${this}
  INCLUDEDIRS ${CMAKE_CURRENT_BINARY_DIR} ${CMAKE_BINARY_DIR}/lib ${include_${this}}
  )
