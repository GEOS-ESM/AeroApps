set (CMAKE_Fortran_FLAGS_RELEASE "-O3 ${GEOS_Fortran_Release_FPE_Flags} ${ALIGNCOM}")

ecbuild_add_executable (
  TARGET argopt
  SOURCES argopt.c
  )

ecbuild_add_executable (
  TARGET fresolv.x
  SOURCES fresolv.f90
  LIBS GMAO_mpeu
  FFLAGS ${FREAL8})

ecbuild_add_executable (
  TARGET gethdftinc.x
  SOURCES gethdftinc.f
  LIBS GMAO_mpeu ${NETCDF_LIBRARIES} Threads::Threads
  INCLUDES ${INC_NETCDF}
  FFLAGS ${FREAL8})

ecbuild_add_executable (
  TARGET zeit_ci.x
  SOURCES zeit_ci.f zeit_usage.f
  LIBS GMAO_mpeu
  FFLAGS ${FREAL8})

ecbuild_add_executable (
  TARGET zeit_co.x
  SOURCES zeit_co.f zeit_usage.f
  LIBS GMAO_mpeu
  FFLAGS ${FREAL8})

ecbuild_add_executable (
  TARGET zeit_pr.x
  SOURCES zeit_pr.f zeit_usage.f zeit_read.f zeit_jday.f
  LIBS GMAO_mpeu
  FFLAGS ${FREAL8})

ecbuild_add_executable (
  TARGET mkdrstdate.x
  SOURCES mkdrstdate.f
  FFLAGS ${FREAL8})

ecbuild_add_executable (
  TARGET rst_date
  SOURCES rst_date.f
  FFLAGS ${FREAL8})

##### Installation Section #####

## bin ##

file(GLOB perl_files *.pl)
# We must remove two files because of DASPERL
list(REMOVE_ITEM perl_files
   ${CMAKE_CURRENT_SOURCE_DIR}/Err_Log.pl 
   ${CMAKE_CURRENT_SOURCE_DIR}/fndate.pl
   )

file(GLOB perl_packages *.pm)
file(GLOB python_files *.py)

set(perlscripts
   acquire
   acquire_obsys
   esma_mpirun
   g5_modules_perl_wrapper
   jobIDfilter
   lsfpL
   mkdir_remote
   n4zip
   n4zip_ods
   now
   num_days_in_month
   pasta
   pesto
   prund.pl
   r_dist
   require_obsys
   rndasfiles
   tick
   token_resolve
   vED
   stripname
   )

set(otherbin
   Functions.csh
   pyrob
   RmShmKeys_sshmpi.csh
   rmshmkeyhere.sh
   )

install (
   PROGRAMS ${perl_files} ${perl_packages} ${python_files} ${perlscripts} ${otherbin}
   DESTINATION bin
   )

set(perlscripts_with_dasperl
   Err_Log.pl
   bacon
   fvgetf
   fndate.pl
   grads2grib
   rget
   rst4ncepf
   rst2rcyc
   rput
   )

set (DASPERL /usr/bin/perl)
foreach (script ${perlscripts_with_dasperl})
   configure_file(${script} ${script} @ONLY)
   install(PROGRAMS ${CMAKE_CURRENT_BINARY_DIR}/${script} DESTINATION bin)
endforeach ()

## etc ##

file(GLOB obsys_rc obsys*.rc)
file(GLOB templates *.tmpl)

install (
   FILES ${obsys_rc} ${templates} pesto.arc PUBLICTAG
   DESTINATION etc
   )

## grads ##

set (grads_files
   lats4d.gs
   lats4d
   info.gs
   )

install (PROGRAMS ${grads_files}
   DESTINATION lib/grads
   )
