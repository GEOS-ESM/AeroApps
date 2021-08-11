esma_set_this()

ecbuild_add_executable (
  TARGET accp_polar_vlidort.x
  SOURCES accp_polar_vlidort.F90
  LIBS GMAO_mpeu GMAO_gfio_r4 MAPL_cfio_r4 Chem_Base MAPL VLIDORT2OS geo_vlidort ${NETCDF_LIBRARIES} MPI::MPI_Fortran
  NOINSTALL
  )

ecbuild_add_executable (
  TARGET accp_polar_vlidort_singlegeom.x
  SOURCES accp_polar_vlidort_singlegeom.F90
  LIBS GMAO_mpeu GMAO_gfio_r4 MAPL_cfio_r4 Chem_Base MAPL VLIDORT2OS geo_vlidort ${NETCDF_LIBRARIES} MPI::MPI_Fortran
  NOINSTALL
  )

install (
    TARGETS accp_polar_vlidort.x accp_polar_vlidort_singlegeom.x 
    DESTINATION bin/missions/${this}
    )

# these are executable scripts
set (PYSCRIPTS
    accp_polar_vlidort.py
    accp_sdist.py
    g5nr_ext_sampler.py
    g5nr_lidar_sampler.py
    merge_lidar.py
    polarimeter_swath.py
    run_accp_polar_vlidort_multinode.py
    run_accp_polar_vlidort.py
    run_accp_polar_vlidort_singlenode.py
    run_mp_lidar_vlidort.py
    run_polarimeter_swath.py
    )

set (SCRIPTS
    clean_mem.sh
    run_accp_polar_vlidort.j
    run_accp_polar_vlidort_multinode.j
    run_accp_polar_vlidort_singlenode.j
    run_accp_sdist.j
    run_ext_sampler.j
    run_lidar_sampler.j
    run_polarimeter_swath.j
    )

install (
    PROGRAMS  ${PYSCRIPTS} ${SCRIPTS}
    DESTINATION bin/missions/${this}
    )

# these are non-executaables
set (PYFILES
    sdist.py
    )

set (RCFILES
    accp_polar_vlidort.rc
    Chem_MieRegistry.rc
    brdf_mcd43c.pcf
    gpm.pcf
    lidar_ext.pcf
    lidar_files.pcf
    lidar_sampler_ler.pcf
    lidar_sampler_ndvi.pcf
    lidar_sampler.pcf
    polar07.pcf
    ss450.pcf
    )

install (
    FILES ${PYFILES} ${RCFILES} ${OTHERFILES}
    DESTINATION bin/missions/${this}
    )

install (
   DIRECTORY rc
   DESTINATION bin/missions/${this}
   USE_SOURCE_PERMISSIONS
   MESSAGE_NEVER
   )

install (
   DIRECTORY tle
   DESTINATION bin/missions/${this}
   USE_SOURCE_PERMISSIONS
   MESSAGE_NEVER
   )

install (
   DIRECTORY sampling
   DESTINATION bin/missions/${this}
   USE_SOURCE_PERMISSIONS
   MESSAGE_NEVER
   )
