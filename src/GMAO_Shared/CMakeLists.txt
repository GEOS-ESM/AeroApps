# GSW (ecbuild project) includes ecbuild_system; disable it!
if( EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/@GSW )
   set (gsw_ECBUILD_SYSTEM_INCLUDED TRUE)
endif()

esma_add_subdirectories (
  @MAPL
  Chem_Base
  GMAO_etc
  GMAO_hermes
  GMAO_gfio
  GMAO_mpeu
  GMAO_ods
  GMAO_psas
  GMAO_pyobs
  )

