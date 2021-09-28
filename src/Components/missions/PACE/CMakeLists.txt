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

ecbuild_add_executable (
  TARGET pace_vlidort_gasabs_multinode.x
  SOURCES pace_vlidort_gasabs_multinode.F90
  LIBS GMAO_mpeu GMAO_gfio_r4 MAPL_cfio_r4 Chem_Base MAPL VLIDORT2OS geo_vlidort ${NETCDF_LIBRARIES} MPI::MPI_Fortran
  NOINSTALL
  )

install (
    TARGETS pace_vlidort.x pace_vlidort_multinode.x pace_vlidort_gasabs_multinode.x
    DESTINATION bin/missions/${this}
    )

# these are executable scripts
set (PYSCRIPTS
    alphaTable.py
    cloud_lc.py
    ext_calculator.py
    geos_to_standard.py
    hyperTest_ck_g5nr.py
    hyperTest_ck_o3_bremen.py
    hyperTest_ck_o3.py
    hyperTest_ck.py
    hyperTest_g5nr.py
    hyperTest_o3_bremen.py
    hyperTest_o3.py
    hyperTest.py
    hyperTest_uv_o3_bremen.py
    leo_sampler.py
    mcd12c_sampler.py
    mcd43c_sampler.py
    myd13c2_sampler.py
    nobm_sampler.py
    omi_ler_sampler.py
    pace_gasabs_lc.py
    pace_lb.py
    pace_lc.py
    pace_lc_to_l1b_mp.py
    pace_lc_to_l1b.py
    pace.py
    sleave_brdf.py
    )

set (SCRIPTS
    clean_mem.sh
    leo_vlidort_lc.j
    pace_lb.j
    pace_lc_array.j
    pace_lc_gas.j
    pace_lc_gas_multinode.j
    pace_lc.j
    pace_lc_multinode.j
    pace_lc_to_l1b.j
    pace_lc_to_l1b_mp.j
    )

install (
    PROGRAMS  ${PYSCRIPTS} ${SCRIPTS}
    DESTINATION bin/missions/${this}
    )

# these are non-executaables
set (PYFILES
    mcd12c.py
    mcd43c.py
    myd13c2.py
    nobm.py
    omi_ler.py
    )

set (RCFILES
    aer_Nv.rc
    Aod_EOS.rc
    asm_Nx.rc
    Chem_MieRegistry.rc
    chm_Nv.rc
    cld_nearest.rc
    cloud_calc.rc
    leo_sampler.rc
    leo_vlidort_cloud.rc
    met_Nv.rc
    met_Nv_uv.rc
    pace_vlidort_gas_470-nccs.rc
    pace_vlidort.rc
    variablesGCS.rc
    )

set (OTHERFILES
    g5nr.pcf
    g5nr_uv.pcf
    mcd12c.pcf
    mcd43c.pcf
    myd13c2.pcf
    nobm.pcf
    omi_ler.pcf
    pace_lb.pcf
    )

install (
    FILES ${PYFILES} ${RCFILES} ${OTHERFILES}
    DESTINATION bin/missions/${this}
    )

install (
   DIRECTORY correlated_k
   DESTINATION bin/missions/${this}
   USE_SOURCE_PERMISSIONS
   MESSAGE_NEVER
   )

install (
   DIRECTORY oci_tables
   DESTINATION bin/missions/${this}
   USE_SOURCE_PERMISSIONS
   MESSAGE_NEVER
   )

install (
   DIRECTORY hyperTest
   DESTINATION bin/missions/${this}
   USE_SOURCE_PERMISSIONS
   MESSAGE_NEVER
   )

