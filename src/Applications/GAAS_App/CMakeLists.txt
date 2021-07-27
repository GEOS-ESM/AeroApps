esma_set_this()

set(srcs aod_util.F90  m_ana.F90  m_duplicate.F  m_obs.F90  m_simul.F90  m_sqc.F90)
set_property(SOURCE ${srcs} PROPERTY COMPILE_FLAGS ${EXTENDED_SOURCE})

esma_add_library(${this}
  SRCS ${srcs}
  DEPENDENCIES GMAO_mpeu GMAO_ods GMAO_psas GMAO_hermes GMAO_gfio_r4
  Chem_Base MAPL
  INCLUDES ${INC_ESMF} ${INC_NETCDF}
  )
target_link_libraries(${this} PRIVATE OpenMP::OpenMP_Fortran)

set_property(SOURCE ana_aod.F PROPERTY COMPILE_FLAGS ${EXTENDED_SOURCE})
ecbuild_add_executable (
  TARGET ana_aod.x
  SOURCES ana_aod.F
  LIBS ${this}
  )

ecbuild_add_executable (
  TARGET mpana_aod.x
  SOURCES mpana_aod.F90
  LIBS ${this}
  )

ecbuild_add_executable (
  TARGET odsqc.x
  SOURCES odsqc.F90
  LIBS ${this}
  )

# these are executable scripts
set (PYSCRIPTS
    aerdt_l2a.py
    aeronet_all.py
    aod_data.py
    avhrr_all.py
    avhrr_l2a.py
    geo04_l2a.py
    misr_land.py
    modis_l2a.py
    mxd04_l2a.py
    patmosx_l2a.py
    )

set (SCRIPTS
    ana_aod.csh
    get_aero_obs.csh
    l2a_all.csh
    run_all.csh
    run_gaas_ana.csh
    ut_get_aero_obs.j
    )

install (
    PROGRAMS  ${PYSCRIPTS} ${SCRIPTS}
    DESTINATION bin/${this}
    )

# these are non-executaables
set (PYFILES
    avhrr_nnr.py
    geo04_nnr.py
    mxd04_nnr.py
    )

set (RCFILES
    ana_aodens.rc
    ana.rc
    mpana.rc
    obs.rc
    psas.rc
    sqc.rc
    )

set (PCFFILES
    aerdt_l2a.pcf
    avhrr_l2a.pcf
    modis_l2a.pcf
    )

install (
    FILES ${PYFILES} ${RCFILES} ${PCFFILES}
    DESTINATION bin/${this}
    )


set (FILEDIRS
    NCA
    )

install (
    DIRECTORY ${FILEDIRS}
    DESTINATION bin/${this}
    USE_SOURCE_PERMISSIONS
    MESSAGE_NEVER
    )
