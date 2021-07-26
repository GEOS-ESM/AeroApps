esma_set_this ()

set (srcs glint_mod.F90  sgp4_mod.F90  TLE_mod.F90)

esma_add_library (${this}
  SRCS ${srcs}
  DEPENDENCIES Chem_Base MAPL GMAO_mpeu GMAO_ods GMAO_gfio_r4
  )

ecbuild_add_executable(
    TARGET ut_sgp4.x
    SOURCES ut_sgp4.F90
    LIBS ${this}
    )

find_package(F2PY2 REQUIRED)
esma_add_f2py2_module (sgp4_
    SOURCES sgp4_py.F90
    DESTINATION lib/Python/${this}
    LIBRARIES ${this}
    INCLUDEDIRS ${CMAKE_CURRENT_BINARY_DIR} ${CMAKE_BINARY_DIR}/lib ${include_${this}}
    )
add_dependencies(sgp4_ ${this})

# these are executable scripts
set (PYSCRIPTS
    tle_sampler.py
    )

set (SCRIPTS
    wget_tle
    )

install (
    PROGRAMS  ${PYSCRIPTS} ${SCRIPTS}
    DESTINATION bin/${this}
    )

# these are non-executaables
set (RCFILES
    tle_sampler.rc
    )

install (
    FILES ${RCFILES}
    DESTINATION bin/${this}
    )

set (FILEDIRS
    TLE
    )

install (
    DIRECTORY ${FILEDIRS}
    DESTINATION bin/${this}
    USE_SOURCE_PERMISSIONS
    MESSAGE_NEVER
    ) 
