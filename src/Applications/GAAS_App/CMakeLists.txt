esma_set_this()

set(srcs aod_util.F90  m_ana.F90  m_duplicate.F  m_obs.F90  m_simul.F90  m_sqc.F90)
set_property(SOURCE ${srcs} PROPERTY COMPILE_FLAGS ${EXTENDED_SOURCE})

esma_add_library(${this}
  SRCS ${srcs}
  DEPENDENCIES GMAO_mpeu GMAO_ods GMAO_psas GMAO_hermes GMAO_gfio_r4
  Chem_Base MAPL
  INCLUDES ${INC_ESMF} ${INC_NETCDF}
  )


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

