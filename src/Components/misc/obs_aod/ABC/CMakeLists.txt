esma_set_this ()

find_package(F2PY3 REQUIRED)
esma_add_f2py3_module (VLIDORT_BRDF_ABC_
  SOURCES VLIDORT_BRDF_ABC_py.F90
  DESTINATION lib/Python/${this}
  LIBRARIES VLIDORT90
  INCLUDEDIRS ${CMAKE_CURRENT_BINARY_DIR} ${CMAKE_BINARY_DIR}/lib ${include_${this}} ${esma_include}/VLIDORT90)
add_dependencies(VLIDORT_BRDF_ABC_ VLIDORT90)

# these are executable scripts
set (PYSCRIPTS
    mcd43c_download.py
    )

install (
    PROGRAMS  ${PYSCRIPTS}
    DESTINATION bin/${this}
    )

# these are non-executaables
set (PYFILES
    abc_avhrr.py
    abc_c6_aux.py
    abc_c6.py
    abc_misr.py
    abc_modis.py
    albedo.py
    anet.py
    aqc_deep.py
    aqc_mapss.py
    brdf.py
    colloc_deep.py
    colloc_misr.py
    evaluate_avhrr.py
    evaluate_cdr.py
    giant.py
    man.py
    mcd43c.py
    nn.py
    reduce_avhrr.py
    synoptic_avhrr.py
    )

set (CTLFILES
    albedo_clim.ctl
    MCD43C1.006.daily.ctl
    MOD13C2.006.monthly.ctl
    MYD13C2.006.monthly.ctl
    tavg1_2d_aer_Nx
    tavg1_2d_slv_Nx
    )

set (DDFFILES
    merra_flx-3hourly.ddf
    merra_flx-hourly.ddf
    merra_slv-hourly.ddf
    misr.ddf
    nnr_001.modo.ddf
    nnr_001.mydo.ddf
    )

install (
    FILES ${PYFILES} ${CTLFILES} ${DDFFILES}
    DESTINATION bin/${this}
    )
  
