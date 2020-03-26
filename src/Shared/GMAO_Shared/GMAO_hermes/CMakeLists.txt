esma_set_this()

add_definitions (-DHERMES -DDEBUG_SHAVE -Dunix -D__unix__)
include_directories (${include_GEOS_Shared})
include_directories (${include_GMAO_gfio_r8})

string (REPLACE " " ";" flags ${FREAL8})
add_compile_options(${flags})

set (srcs
  m_ana2dyn.F m_const.f90 m_dyn.f90 m_dyn2dyn.f90 m_dynp.F90
  m_insitu.F m_interp.F90 m_lp.F m_maph.F90 m_maph_pert.F90 m_mapz.F m_mapz_pert.F90
  m_qsat.F m_remap.f90 m_set_eta.F90 m_vdc2vdc.F90 progsig.f90 m_topo_remap.F90
  m_interpack.F m_interpack_ad.F m_ec_set_eta.F90 m_ecdyn.f90 m_tick.f90
  m_dyn_util.F90 m_spline.f90
  )

esma_add_library(${this} SRCS ${srcs} DEPENDENCIES GMAO_mpeu GMAO_gfio_r8)
target_compile_options(${this} PRIVATE ${flags})

set (CMAKE_Fortran_FLAGS_RELEASE "-O3 ${common_Fortran_flags} ${ALIGNCOM}")

if (EXTENDED_SOURCE)
  set_source_files_properties (m_interpack.F    PROPERTIES COMPILE_FLAGS ${EXTENDED_SOURCE})
  set_source_files_properties (m_interpack_ad.F PROPERTIES COMPILE_FLAGS ${EXTENDED_SOURCE})
  set_source_files_properties (m_ana2dyn.F      PROPERTIES COMPILE_FLAGS ${EXTENDED_SOURCE})
  set_source_files_properties (m_insitu.F       PROPERTIES COMPILE_FLAGS ${EXTENDED_SOURCE})
  set_source_files_properties (m_set_eta.F90    PROPERTIES COMPILE_FLAGS ${EXTENDED_SOURCE})
  set_source_files_properties (ibc_upd.F        PROPERTIES COMPILE_FLAGS ${EXTENDED_SOURCE})
  set_source_files_properties (gcmbkg2ana.F     PROPERTIES COMPILE_FLAGS ${EXTENDED_SOURCE})
  set_source_files_properties (diag2dyn.F       PROPERTIES COMPILE_FLAGS ${EXTENDED_SOURCE})
endif ()
set_source_files_properties (GFIO_Shave.f90 PROPERTIES COMPILE_FLAGS ${PP})
set_source_files_properties (GFIO_remap.f90 PROPERTIES COMPILE_FLAGS ${PP})

ecbuild_add_executable (TARGET GFIO_remap.x     SOURCES GFIO_remap.f90     LIBS ${this} MAPL_cfio_r8)
ecbuild_add_executable (TARGET GFIO_Shave.x     SOURCES GFIO_Shave.f90     LIBS ${this} MAPL_cfio_r8)
ecbuild_add_executable (TARGET lcv2prs.x        SOURCES lcv2prs.F90        LIBS ${this} MAPL_cfio_r8 MAPL GEOS_Shared)

ecbuild_add_executable (TARGET ana5sfc.x        SOURCES ana5sfc.F90        LIBS ${this})
ecbuild_add_executable (TARGET diag2dyn.x       SOURCES diag2dyn.F         LIBS ${this})
ecbuild_add_executable (TARGET drs2dyn.x        SOURCES drs2dyn.f90        LIBS ${this})
ecbuild_add_executable (TARGET dyn2drs.x        SOURCES dyn2drs.f90        LIBS ${this})
ecbuild_add_executable (TARGET dyn2dyn.x        SOURCES dyn2dyn.f90        LIBS ${this})
ecbuild_add_executable (TARGET dyn_boot.x       SOURCES dyn_boot.f90       LIBS ${this})
ecbuild_add_executable (TARGET dyn2prs.x        SOURCES dyn2prs.f90        LIBS ${this})
ecbuild_add_executable (TARGET dyn2rs5.x        SOURCES dyn2rs5.f90        LIBS ${this})
ecbuild_add_executable (TARGET dyn52dyn.x       SOURCES dyn52dyn.f90       LIBS ${this})
ecbuild_add_executable (TARGET dyndiff.x        SOURCES dyndiff.f90        LIBS ${this})
ecbuild_add_executable (TARGET dyndims.x        SOURCES dyndims.f          LIBS ${this})
ecbuild_add_executable (TARGET dynp.x           SOURCES dynp.f90           LIBS ${this})
ecbuild_add_executable (TARGET echorc.x         SOURCES echorc.f90         LIBS ${this})
ecbuild_add_executable (TARGET fvsubset.x       SOURCES fvsubset.F90       LIBS ${this})
ecbuild_add_executable (TARGET getgfiodim.x     SOURCES getgfiodim.f90     LIBS ${this})
ecbuild_add_executable (TARGET ibc_upd.x        SOURCES ibc_upd.F          LIBS ${this})
ecbuild_add_executable (TARGET maph.x           SOURCES maph.f90           LIBS ${this})
ecbuild_add_executable (TARGET maph_pert.x      SOURCES maph_pert.f90      LIBS ${this})
ecbuild_add_executable (TARGET mapz.x           SOURCES mapz.f             LIBS ${this})
ecbuild_add_executable (TARGET recalcsfc.x      SOURCES recalcsfc.F90      LIBS ${this})
ecbuild_add_executable (TARGET reset_time.x     SOURCES reset_time.f       LIBS ${this})
ecbuild_add_executable (TARGET rs52dyn.x        SOURCES rs52dyn.f90        LIBS ${this})
ecbuild_add_executable (TARGET gcmbkg2ana.x     SOURCES gcmbkg2ana.F       LIBS ${this})
ecbuild_add_executable (TARGET dyndot.x         SOURCES dyndot.f90         LIBS ${this})
ecbuild_add_executable (TARGET dyn_pert_remap.x SOURCES dyn_pert_remap.F90 LIBS ${this})
ecbuild_add_executable (TARGET dyn_recenter.x   SOURCES dyn_recenter.f90   LIBS ${this})
ecbuild_add_executable (TARGET dyn_rankhist.x   SOURCES dyn_rankhist.F90   LIBS ${this})
ecbuild_add_executable (TARGET dyn2real_eta.x   SOURCES dyn2real_eta.f90   LIBS ${this})
ecbuild_add_executable (TARGET dyn_ncf2dyn.x    SOURCES dyn_ncf2dyn.f90    LIBS ${this})
ecbuild_add_executable (TARGET dyn_iupd.x       SOURCES dyn_iupd.F90       LIBS ${this})
ecbuild_add_executable (TARGET dyn_cov.x        SOURCES dyn_cov.f90        LIBS ${this})
ecbuild_add_executable (TARGET dyn_efsens.x     SOURCES dyn_efsens.f90     LIBS ${this})
ecbuild_add_executable (TARGET eta_echo.x       SOURCES eta_echo.f90       LIBS ${this})
ecbuild_add_executable (TARGET dyn_inflate.x    SOURCES dyn_inflate.f90    LIBS ${this})

install(PROGRAMS cnv2prs.pl ana5sfc DESTINATION bin)

file(GLOB resource_files *.rc)
install(FILES ${resource_files} DESTINATION etc)
