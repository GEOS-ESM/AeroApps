esma_set_this()

ecbuild_add_executable (
    TARGET pace_vlidort.x
    SOURCES pace_vlidort.F90
  LIBS GMAO_mpeu GMAO_gfio_r4 MAPL_cfio_r4 Chem_Base MAPL VLIDORT2OS geo_vlidort ${NETCDF_LIBRARIES} MPI::MPI_Fortran
  NOINSTALL
  )

ecbuild_add_executable (
  TARGET pace_vlidort_multinode.x
  SOURCES pace_vlidort_multinode.F90
  LIBS GMAO_mpeu GMAO_gfio_r4 MAPL_cfio_r4 Chem_Base MAPL VLIDORT2OS geo_vlidort ${NETCDF_LIBRARIES} MPI::MPI_Fortran
  NOINSTALL
  )

#ecbuild_add_executable (
#  TARGET pace_vlidort_gasabs_multinode.x
#  SOURCES pace_vlidort_gasabs_multinode.F90
#  LIBS GMAO_mpeu GMAO_gfio_r4 MAPL_cfio_r4 Chem_Base MAPL VLIDORT2OS geo_vlidort ${NETCDF_LIBRARIES} MPI::MPI_Fortran
#  NOINSTALL
#  )

install (
    TARGETS pace_vlidort.x pace_vlidort_multinode.x
    DESTINATION bin/missions/${this}
    )

# these are executable scripts
#set (PYSCRIPTS
#    )

#set (SCRIPTS
#    )

#install (
#    PROGRAMS  ${PYSCRIPTS} ${SCRIPTS}
#    DESTINATION bin/missions/${this}
#    )

# these are non-executaables
#set (PYFILES
#    sdist.py
#    )

#set (RCFILES
#    )

#install (
#    FILES ${PYFILES} ${RCFILES} ${OTHERFILES}
#    DESTINATION bin/missions/${this}
#    )

#install (
#   DIRECTORY rc
#   DESTINATION bin/missions/${this}
#   USE_SOURCE_PERMISSIONS
#   MESSAGE_NEVER
#   )

#install (
#   DIRECTORY tle
#   DESTINATION bin/missions/${this}
#   USE_SOURCE_PERMISSIONS
#   MESSAGE_NEVER
#   )

#install (
#   DIRECTORY sampling
#   DESTINATION bin/missions/${this}
#   USE_SOURCE_PERMISSIONS
#   MESSAGE_NEVER
#   )

#install (
#   DIRECTORY sampling_M2R12K
#   DESTINATION bin/missions/${this}
#   USE_SOURCE_PERMISSIONS
#   MESSAGE_NEVER
#   )
