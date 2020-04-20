esma_set_this ()

set (srcs glint_mod.F90  sgp4_mod.F90  TLE_mod.F90)

esma_add_library (${this}
  SRCS ${srcs}
  DEPENDENCIES Chem_Base MAPL GMAO_mpeu GMAO_ods gfio_r4
  )
