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
