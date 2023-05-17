esma_set_this ()

set (MODSRC
  m_ods.f90 m_roms.f90 m_odsmeta.f m_ods_structure.f m_ods_obsdiags.F90
  m_saber.f m_odsxsup.f90 m_Sndx.F90 m_ODS_Sample.f m_mraob.f90
  )

set (ODSSRC
  ods_cal2time.f  ods_cal2min.f
  ods_julian.f    ods_ncp_itoi4.f  ods_rget.f
  ods_caldat.f    ods_namemonth.f  ods_ncp_itor4.f  ods_scalei.f
  ods_case.f      ods_ncagtc.f     ods_ncp_rtoi1.f  ods_scaleirev.f
  ods_cget.f      ods_ncagti.f     ods_ncp_rtoi2.f  ods_scaler.f
  ods_checki.f    ods_ncagtr.f     ods_ncp_rtoi4.f  ods_scalerrev.f
  ods_checkr.f    ods_ncapti.f     ods_ncp_rtor4.f  ods_setparmc.f
  ods_coards.f    ods_ncaptr.f     ods_ncvgti.f     ods_setparmi.f
  ods_defwsp.f    ods_nccreate.f   ods_ncvgtr.f     ods_setparmr.f
  ods_file.f      ods_ncg_i1toi.f  ods_ncvpti.f     ods_strsearch.f
  ods_getattp.f   ods_ncg_i1tor.f  ods_ncvptr.f     ods_strsize.f
  ods_getc.f      ods_putc.f
  ods_ncvptc.f    ods_ncvgtc.f
  ods_geti.f      ods_ncg_i2toi.f  ods_ncvwsp.f     ods_time2cal.f
  ods_min2cal.f   ods_mraob.f
  ods_getlist.f   ods_ncg_i2tor.f  ods_nget.f       ods_type.f
  ods_getp.f      ods_ncg_i4toi.f  ods_nval.f       ods_updatep.f
  ods_getr.f      ods_ncg_i4tor.f  ods_parm.f       ods_valwsp.f
  ods_handle.f    ods_ncg_r4toi.f  ods_puti.f       ods_vartype.f
  ods_iget.f      ods_ncg_r4tor.f  ods_putlist.f    ods_verindex.f
  ods_imonth.f    ods_ncopen.f     ods_putr.f       ods_vertag.f
  ods_ncp_itoi1.f ods_readp.f      ods_writep.f
  ods_julhr.f     ods_ncp_itoi2.f  ods_resetp.f     odsnxtime.f
  ods_eget.f      ods_dcget.f      ods_dconv.f      ods_dcscan.f
  ods_drad.f      ods_dsbuv.f      ods_dpcp.f       ods_rmdupl.f
  ods_do3lev.f 
  )

# This is equivalent to FOPT=$(FOPT3) in GNU Make
if (CMAKE_Fortran_COMPILER_ID MATCHES Intel)
   string (REPLACE "${OPTREPORT0}" "" CMAKE_Fortran_FLAGS_RELEASE ${CMAKE_Fortran_FLAGS_RELEASE})
   string (REPLACE "${FTZ}" "" CMAKE_Fortran_FLAGS_RELEASE ${CMAKE_Fortran_FLAGS_RELEASE})
   string (REPLACE "${ALIGN_ALL}" "" CMAKE_Fortran_FLAGS_RELEASE ${CMAKE_Fortran_FLAGS_RELEASE})
   string (REPLACE "${NO_ALIAS}" "" CMAKE_Fortran_FLAGS_RELEASE ${CMAKE_Fortran_FLAGS_RELEASE})
endif ()
esma_add_library (${this}
  SRCS ${MODSRC} ${ODSSRC}
  DEPENDENCIES GMAO_mpeu GMAO_pilgrim
  INCLUDES ${INC_NETCDF}
  )

ecbuild_add_executable (TARGET ods_scan.x SOURCES ods_scan.f LIBS ${this} ${NETCDF_LIBRARIES})
ecbuild_add_executable (TARGET dconv2ods.x SOURCES dconv2ods.f LIBS ${this} ${NETCDF_LIBRARIES})
ecbuild_add_executable (TARGET ods_blank.x SOURCES ods_blank.f90 LIBS ${this} ${NETCDF_LIBRARIES})
ecbuild_add_executable (TARGET ods_maker.x SOURCES ods_maker.f90 LIBS ${this} ${NETCDF_LIBRARIES})
ecbuild_add_executable (TARGET ods_sample.x SOURCES ods_sample.f LIBS ${this} ${NETCDF_LIBRARIES})
ecbuild_add_executable (TARGET ods_synX2Y.x SOURCES ods_synX2Y.f90 LIBS ${this} ${NETCDF_LIBRARIES})
ecbuild_add_executable (TARGET ods_tally.x SOURCES ods_tally.f90 LIBS ${this} ${NETCDF_LIBRARIES})
ecbuild_add_executable (TARGET saber2ods.x SOURCES saber2ods.f LIBS ${this} ${NETCDF_LIBRARIES})
ecbuild_add_executable (TARGET odslist SOURCES odslist.f LIBS ${this} ${NETCDF_LIBRARIES})
ecbuild_add_executable (TARGET odsmatch SOURCES odsmatch.f LIBS ${this} ${NETCDF_LIBRARIES})
ecbuild_add_executable (TARGET odsselect SOURCES odsselect.f LIBS ${this} ${NETCDF_LIBRARIES})
ecbuild_add_executable (TARGET odsshuffle SOURCES odsshuffle.f LIBS ${this} ${NETCDF_LIBRARIES})
ecbuild_add_executable (TARGET odsstats SOURCES odsstats.f LIBS ${this} ${NETCDF_LIBRARIES})
ecbuild_add_executable (TARGET odsemean SOURCES odsemean.f LIBS ${this} ${NETCDF_LIBRARIES})

if (EXTENDED_SOURCE)
  set_target_properties (${this} PROPERTIES COMPILE_FLAGS ${EXTENDED_SOURCE})
  set_target_properties (saber2ods.x PROPERTIES COMPILE_FLAGS ${EXTENDED_SOURCE})
  set_target_properties (dconv2ods.x PROPERTIES COMPILE_FLAGS ${EXTENDED_SOURCE})
  set_target_properties (odsmatch  PROPERTIES COMPILE_FLAGS ${EXTENDED_SOURCE})
  set_target_properties (odslist  PROPERTIES COMPILE_FLAGS ${EXTENDED_SOURCE})
  set_target_properties (odsselect  PROPERTIES COMPILE_FLAGS ${EXTENDED_SOURCE})
  set_target_properties (odsshuffle  PROPERTIES COMPILE_FLAGS ${EXTENDED_SOURCE})
  set_target_properties (odsstats  PROPERTIES COMPILE_FLAGS ${EXTENDED_SOURCE})
  set_target_properties (odsemean  PROPERTIES COMPILE_FLAGS ${EXTENDED_SOURCE})
endif()

install(
   PROGRAMS unique_instruments.csh obimp_summary.pl diag2ods
   DESTINATION bin
   )

file (GLOB rc_files *.rc)
install (
   FILES ${rc_files} odslist.odt
   DESTINATION etc
   )

install (
   PROGRAMS pyods/__init__.py
   DESTINATION lib/Python/pyods
   )

find_package(F2PY3)
if (F2PY3_FOUND)
  # Using add_ here as the import test threw an SSL error on GMAO desktop
  add_f2py3_module(pyods_
    SOURCES pyods_.F90
    LIBRARIES GMAO_ods GMAO_mpeu ${NETCDF_LIBRARIES}
    INCLUDEDIRS ${CMAKE_CURRENT_BINARY_DIR} ${CMAKE_BINARY_DIR}/lib ${BASEDIR}/lib ${include_${this}} ${INC_NETCDF}
    DESTINATION lib/Python/pyods
    USE_MPI
    )
  add_dependencies(pyods_ ${this})
endif ()

separate_arguments(ODSMETA_ARGS NATIVE_COMMAND "sed -e 's/integer, parameter//g' -e 's/:://g' -e 's/\!/   #/g' -e 's/^[ ]*//g'")
add_custom_command (
   OUTPUT odsmeta.py
   WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
   SOURCE ${CMAKE_CURRENT_SOURCE_DIR}/odsmeta.h
   COMMAND grep 'integer, parameter' ${CMAKE_CURRENT_SOURCE_DIR}/odsmeta.h | ${ODSMETA_ARGS} | grep -v idsats > odsmeta.py
   )

add_custom_target(make_odsmeta ALL DEPENDS odsmeta.py)

install(
   PROGRAMS ${CMAKE_CURRENT_BINARY_DIR}/odsmeta.py
   DESTINATION lib/Python/pyods
   )

