esma_add_subdirectories (solar_reference)
esma_set_this ()

set (srcs 
    sys_tracker.c 
    GeoAngles_py.F90 
    VLIDORT_LAMBERT.F90 
    VLIDORT_BRDF_MODIS.F90 
    VLIDORT_BRDF_MODIS_BPDF.F90
    LIDORT_LAMBERT.F90 
    LIDORT_BRDF_MODIS.F90 
    VLIDORT_SURFACE.F90 
    SURFACE.F90 
    VLIDORT_MieMod.F90
    VLIDORT_BRDF_CX.F90 
    VLIDORT_LAMBERT_BPDF.F90 
    VLIDORT_ROT.F90
    mp_netcdf_Mod.F90 
    netcdf_Mod.F90
    cloud_MieMod.F90)


esma_add_library (${this}
  SRCS ${srcs}
  DEPENDENCIES MAPL.shared MAPL.base VLIDORT2OS LIDORT3p7 Chem_Base ${NETCDF_LIBRARIES} MPI::MPI_Fortran

  )

ecbuild_add_executable (
  TARGET geo_vlidort.x
  SOURCES geo_vlidort.F90
  LIBS GMAO_mpeu GMAO_gfio_r4 MAPL_cfio_r4 Chem_Base MAPL VLIDORT2OS geo_vlidort ${NETCDF_LIBRARIES} MPI::MPI_Fortran
  )

ecbuild_add_executable (
  TARGET geo_vlidort_cloud.x
  SOURCES geo_vlidort_cloud.F90
  LIBS GMAO_mpeu GMAO_gfio_r4 MAPL_cfio_r4 Chem_Base MAPL VLIDORT2OS geo_vlidort ${NETCDF_LIBRARIES} MPI::MPI_Fortran
  )

# this needs de-bugging/updates to new driver interfaces
#ecbuild_add_executable (
#  TARGET geo_vlidort_vnncLUT.x
#  SOURCES geo_vlidort_vnncLUT.F90
#  LIBS GMAO_mpeu GMAO_gfio_r4 MAPL_cfio_r4 Chem_Base MAPL VLIDORT2OS geo_vlidort ${NETCDF_LIBRARIES} MPI::MPI_Fortran
#  )

ecbuild_add_executable (
  TARGET geo_lidort.x
  SOURCES geo_lidort.F90 netcdf_helper.F90
  LIBS GMAO_mpeu GMAO_gfio_r4 MAPL_cfio_r4 Chem_Base MAPL LIDORT3p7 geo_vlidort ${NETCDF_LIBRARIES} MPI::MPI_Fortran
  )

# this needs de-bugging/updates to new driver interfaces
#ecbuild_add_executable (
#  TARGET geo_vlidort_multinode.x
#  SOURCES geo_vlidort_multinode.F90 netcdf_helper.F90
#  LIBS GMAO_mpeu GMAO_gfio_r4 MAPL_cfio_r4 Chem_Base MAPL VLIDORT2OS geo_vlidort ${NETCDF_LIBRARIES} MPI::MPI_Fortran
#  )

ecbuild_add_executable (
  TARGET geo_vlidort_surface.x
  SOURCES geo_vlidort_surface.F90 netcdf_helper.F90
  LIBS GMAO_mpeu GMAO_gfio_r4 MAPL_cfio_r4 Chem_Base MAPL VLIDORT2OS geo_vlidort ${NETCDF_LIBRARIES} MPI::MPI_Fortran
  )

ecbuild_add_executable (
  TARGET geo_surface.x
  SOURCES geo_surface.F90 netcdf_helper.F90
  LIBS GMAO_mpeu GMAO_gfio_r4 MAPL_cfio_r4 Chem_Base MAPL geo_vlidort ${NETCDF_LIBRARIES} MPI::MPI_Fortran
  )

ecbuild_add_executable (
  TARGET geo_angles.x
  SOURCES geo_angles.F90
  LIBS GMAO_mpeu GMAO_gfio_r4 MAPL_cfio_r4 Chem_Base MAPL geo_vlidort ${NETCDF_LIBRARIES} MPI::MPI_Fortran
  )

# these are executable scripts
set (PYSCRIPTS
    geo_vlidort_AI_lc2.py
    geo_vlidort_blksurf_lc2.py
    geo_vlidort_cloud_lc2.py
    geo_vlidort_lc2.py
    geo_vlidort_ocean_lc2.py
    geo_warp.py
    run_geo_vlidort_AI_lc2.py
    run_geo_vlidort_blksurf_lc2.py
    run_geo_vlidort_cloud_lc2.py
    run_geo_vlidort_lc2.py
    run_geo_vlidort_ocean_lc2.py
    utAngles.py
    )

install (
    PROGRAMS  ${PYSCRIPTS}
    DESTINATION bin/${this}
    )


