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


find_package(F2PY2 REQUIRED)

esma_add_f2py2_module (Aero_
  SOURCES Aero_py.F90
  LIBRARIES MAPL.shared MAPL.base
  DESTINATION lib/Python/${this}
  INCLUDEDIRS ${CMAKE_CURRENT_BINARY_DIR} ${CMAKE_BINARY_DIR}/lib ${include_${this}} ${esma_include}/MAPL.shared ${esma_include}/MAPL.base
  )
add_dependencies(Aero_ MAPL.shared MAPL.base)

# MAT: Using add_f2py2 as I haven't figured out how to add MKL to python import test (needs vscdfnorm_)
add_f2py2_module (ICA_
  SOURCES ICA_py.F90
  DESTINATION lib/Python/${this}
  LIBRARIES ${this} MAPL.base
  INCLUDEDIRS ${CMAKE_CURRENT_BINARY_DIR} ${CMAKE_BINARY_DIR}/lib ${include_${this}} ${esma_include}/MAPL.base
  )
add_dependencies(ICA_ MAPL.base ${this})

esma_add_f2py2_module (mcsRegrid_
  SOURCES mcsRegrid_py.F90
  DESTINATION lib/Python/${this}
  LIBRARIES ${this}
  INCLUDEDIRS ${CMAKE_CURRENT_BINARY_DIR} ${CMAKE_BINARY_DIR}/lib ${include_${this}}
  )
add_dependencies(mcsRegrid_ ${this})

# these are executable scripts
set (PYSCRIPTS
    mcs_aero.py
    mcs.py
    renn.py
    )

install (
    PROGRAMS  ${PYSCRIPTS}
    DESTINATION bin/${this}
    )

# these are non-executaables
set (PYFILES
    eta.py
    ica.py
    )

set (RCFILE
    Aod_MODIS_Clouds.rc
    variables.rc
    )

set (OTHERFILES    
    edge-ica.ddf
    edge-noica.ddf
    README
    Re-Exploratory.ipynb
    )

install (
    FILES ${PYFILES} ${RCFILES} ${OTHERFILES}
    DESTINATION bin/${this}
    )
