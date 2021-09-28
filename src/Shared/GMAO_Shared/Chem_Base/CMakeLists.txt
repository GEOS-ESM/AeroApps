esma_set_this ()

set (srcs
  Chem_RegistryMod.F90
  Chem_ArrayMod.F90
  Chem_BundleMod.F90
  Chem_Mod.F90
  Chem_InitMod.F90
  Chem_MieTableMod.F90
  Chem_MieMod.F90
  Chem_AodMod.F90
  Chem_SimpleBundleMod.F90
  MieObs_py.F90
  )
esma_add_library(${this} SRCS ${srcs} DEPENDENCIES GMAO_mpeu MAPL GMAO_gfio_r4)
target_include_directories (${this} PUBLIC ${INC_ESMF} ${INC_NETCDF})
target_compile_definitions (${this} PRIVATE GEOS5)

ecbuild_add_executable(TARGET gogo.x SOURCES gogo.F90 LIBS ${this})

find_package(F2PY2 REQUIRED)
# MAT: Note uses add_f2py2 because I can't figure out the import test for netCDF f2py code
add_f2py2_module(MieObs_ 
   SOURCES MieObs_py.F90 
   DESTINATION lib/Python
   LIBRARIES Chem_Base MAPL GMAO_mpeu ${NETCDF_LIBRARIES} ${ESMF_LIBRARY}
   INCLUDEDIRS ${CMAKE_CURRENT_BINARY_DIR} ${CMAKE_BINARY_DIR}/lib ${esma_include}/Chem_Base ${BASEDIR}/lib ${include_${this}} ${INC_NETCDF}
   USE_MPI
   )
add_dependencies(MieObs_ Chem_Base GMAO_mpeu ${this})


foreach (exe Chem_Aod.x Chem_Aod3d.x ctl_crst.x Chem_BundleToG5rs.x reff_calculator.xx ext_calculator.xx)
  get_filename_component (basename ${exe} NAME_WE) # base file name
  ecbuild_add_executable (
    TARGET ${exe}
    SOURCES ${basename}.F90
    LIBS ${this} ${OpenMP_Fortran_LIBRARIES})
  # Do not know how to propagate the flags for OpenMP link step from lower libraries.
  set_target_properties(${exe} PROPERTIES LINK_FLAGS "${OpenMP_Fortran_FLAGS}")
endforeach ()

# Copy RC files to build tree
file (GLOB_RECURSE rc_files RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} *.rc)
foreach ( file ${rc_files} )
   get_filename_component( dir ${file} DIRECTORY )
   install( FILES ${file} DESTINATION etc/${dir} )
endforeach()

install(PROGRAMS ext_calculator.py reff.py extinction_calculator.csh DESTINATION bin)
install(PROGRAMS mieobs.py DESTINATION lib/Python)
