esma_set_this()

ecbuild_add_executable (
  TARGET accp_polar_vlidort.x
  SOURCES accp_polar_vlidort.F90
  LIBS GMAO_mpeu GMAO_gfio_r4 MAPL_cfio_r4 Chem_Base MAPL ${VLIDORT_LIBRARY} geo_vlidort ${NETCDF_LIBRARIES} MPI::MPI_Fortran
  NOINSTALL
  )

ecbuild_add_executable (
  TARGET accp_polar_vlidort_singlegeom.x
  SOURCES accp_polar_vlidort_singlegeom.F90
  LIBS GMAO_mpeu GMAO_gfio_r4 MAPL_cfio_r4 Chem_Base MAPL ${VLIDORT_LIBRARY} geo_vlidort ${NETCDF_LIBRARIES} MPI::MPI_Fortran
  NOINSTALL
  )

ecbuild_add_executable (
  TARGET accp_polar_vlidort_cirrus.x
  SOURCES accp_polar_vlidort_cirrus.F90 cirrus_MieMod.F90
  LIBS GMAO_mpeu GMAO_gfio_r4 MAPL_cfio_r4 Chem_Base MAPL ${VLIDORT_LIBRARY} geo_vlidort ${NETCDF_LIBRARIES} MPI::MPI_Fortran
  NOINSTALL
  )

install (
    TARGETS accp_polar_vlidort.x accp_polar_vlidort_singlegeom.x accp_polar_vlidort_cirrus.x
    DESTINATION bin/missions/${this}
    )

# these are executable scripts
set (PYSCRIPTS
    accp_polar_vlidort_cirrus.py
    accp_polar_vlidort.py
    accp_sdist.py
    brdf_mcd43c_sampler.py
    fillvalue_fix.py
    fillvalue_fix_surface.py
    g5nr_ext_sampler.py
    g5nr_lidar_sampler.py
    lidar_sampler.py
    link_g5nr.py
    merge_lidar.py
    merge_random_mp.py
    merge_random.py
    plt_maxscat.py
    plt_nadirscat.py
    plt_nadirsolar.py
    plt_scat.py
    polarimeter_swath.py
    polarimeter_swath_single.py
    run_accp_polar_vlidort_cirrus.py
    run_accp_polar_vlidort_multinode.py
    run_accp_polar_vlidort.py
    run_accp_polar_vlidort_singlenode.py
    run_accp_sdist.py
    run_ext_sampler.py
    run_lidar_sampler.py
    run_mp_lidar_vlidort.py
    run_polarimeter_swath.py
    single_orbit.py
    )

set (SCRIPTS
    clean_mem.sh
    run_accp_polar_vlidort_cirrus.j
    run_accp_polar_vlidort.j
    run_accp_polar_vlidort_multinode.j
    run_accp_polar_vlidort_singlenode.j
    run_accp_sdist.j
    run_ext_sampler.j
    run_lidar_sampler.j
    run_lidar_sampler_m2r12k.j
    run_merge_random.g.j
    run_merge_random.g.m2r12k.j
    run_merge_random.ss.j
    run_merge_random.ss.m2r12k.j
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
    Aod_EOS.rc
    Chem_MieRegistry.rc
    brdf_mcd43c.pcf
    gpm045.pcf
    gpm050.pcf
    gpm055.pcf
    gpm060.pcf
    gpm.pcf
    landcover_mcd12c.pcf
    lidar_ext_gpmfine.m2r12k.pcf
    lidar_ext_gpmfine.pcf
    lidar_ext.m2r12k.pcf
    lidar_ext.pcf
    lidar_ext_ss450fine.m2r12k.pcf
    lidar_ext_ss450fine.pcf
    lidar_ext_ss450.pcf
    lidar_files_m2r12k.pcf
    lidar_files.pcf
    lidar_files_random.pcf
    lidar_sampler.gpm.pcf
    lidar_sampler_ler.gpm.pcf
    lidar_sampler_ler.pcf
    lidar_sampler.m2r12k.gpm.pcf
    lidar_sampler.m2r12k.ss450.pcf
    lidar_sampler_ndvi.pcf
    lidar_sampler_ndvi.ss450.pcf
    lidar_sampler.pcf
    merge_files_m2r12k.pcf
    merge_files.pcf
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

install (
   DIRECTORY sampling_M2R12K
   DESTINATION bin/missions/${this}
   USE_SOURCE_PERMISSIONS
   MESSAGE_NEVER
   )
