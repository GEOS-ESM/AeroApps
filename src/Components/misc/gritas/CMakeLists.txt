esma_set_this ()

set (srcs
  etc.f          gritas_gfio.f     gritas_grads.f
  gritas_kxkt.f  m_gritas_grids.f  rdel.f
  )

set_property(SOURCE ${srcs} PROPERTY COMPILE_FLAGS ${EXTENDED_SOURCE})

esma_add_library (${this}
  SRCS ${srcs}
  DEPENDENCIES GMAO_mpeu GMAO_ods GMAO_psas GMAO_hermes
  )

set_property(SOURCE grisas.f PROPERTY COMPILE_FLAGS ${EXTENDED_SOURCE})
ecbuild_add_executable (
  TARGET grisas.x
  SOURCES grisas.f
  LIBS ${this}
  )

set_property(SOURCE gritas.f PROPERTY COMPILE_FLAGS ${EXTENDED_SOURCE})
ecbuild_add_executable (
  TARGET gritas.x
  SOURCES gritas.f
  LIBS ${this}
  )


# these are executable scripts
set (SCRIPTS
    fdp
    )

install (
    PROGRAMS  ${SCRIPTS}
    DESTINATION bin/${this}
    )

# these are non-executaables
set (RCFILES
    gritas_amsua.rc
    gritas_amsub.rc
    gritas_aod.rc
    gritas_goes.rc
    gritas_goessond.rc
    gritas_hirs.rc
    gritas_laod.rc
    gritas_msu.rc
    gritas.rc
    gritas_ssu.rc
    gritas_upconv2.rc
    gritas_upconv.rc
    )    

install (
    FILES ${RCFILES}
    DESTINATION bin/${this}
    )
