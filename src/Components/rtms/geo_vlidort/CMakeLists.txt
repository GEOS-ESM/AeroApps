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
    netcdf_Mod.F90)


esma_add_library (${this}
  SRCS ${srcs}
  DEPENDENCIES MAPL.shared MAPL.base VLIDORT2OS LIDORT3p7 Chem_Base   DEPENDENCIES MAPL.shared MAPL.base VLIDORT2OS LIDORT3p7 Chem_Base ${NETCDF_LIBRARIES} MPI::MPI_Fortran

  )
