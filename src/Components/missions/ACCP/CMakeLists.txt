esma_set_this()

ecbuild_add_executable (
  TARGET accp_polar_vlidort.x
  SOURCES accp_polar_vlidort.F90
  LIBS GMAO_mpeu GMAO_gfio_r4 MAPL_cfio_r4 Chem_Base MAPL VLIDORT2OS geo_vlidort ${NETCDF_LIBRARIES} MPI::MPI_Fortran

  )

#DESTINATION bin/missions/${this}
