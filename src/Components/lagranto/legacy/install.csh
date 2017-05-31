#!/bin/csh

# -----------------------------------------------------------------------------
# Set some general parameters
# -----------------------------------------------------------------------------

set path = ( $BASEDIR/`uname`/bin $path )
    
# Usage
if ( $#argv == 0 ) then
  echo "install.sh [lib|core|goodies|links|all|docu|clean|tag] "
  exit 0
endif

# Set the mode
set mode = $1

# Set path to SVN repository
set svnpath=https://svn.iac.ethz.ch/pub/lagranto.ecmwf/

# Set paths for development and for synchronisation (operational)
#set path_devel = "${DYN_TOOLS}//lagranto.ecmwf/"
#set path_sync  = "${DYN_TOOLS}/lagranto/"

# Set Fortran compiler
setenv FORTRAN  ifort

# Set netcdf format (ive, cdo, mch) 
set ncdf = cdo

# Init netCDF library depending on the Fortran compiler
if ( "${FORTRAN}" == "pgf90" ) then
  module load netcdf/4.2.1-pgf90
  module list

else if ( "${FORTRAN}" == "gfortran" ) then
  module load gfortran
  module load netcdf/4.2.1

else if ( "${FORTRAN}" == "ifort" ) then
  module load ifort/10.1.017
  module load netcdf/4.1.1-ifort

else
  echo "Fortran Compiler ${FORTRAN} not supported... Exit"
  exit 1

endif

# -----------------------------------------------------------------------------
# Create a new tag in SVN repository
# -----------------------------------------------------------------------------

if ( "${mode}" == "tag" ) then
   svn info
   if ( "${#argv}" != 2 ) then
     echo "Usage: install.csh tag id <id=tag number>"
   else
     set tagnr = $2
   endif
   svn copy ${svnpath}/trunk ${svnpath}/tags/${tagnr} -m "Release ${tagnr}"
   exit 0
endif

# -----------------------------------------------------------------------------
# Set internal parameters and detailed installation mode
# -----------------------------------------------------------------------------

# Set LAGRANTO environment variable
#setenv LAGRANTO ${path_devel}

# Set netCDF paths
setenv NETCDF_LIB "`nc-config --flibs` -L/opt/mpi/ifort/12.1/bin//../lib -lmpi -lmpi_cxx -lmpi_f77 -limf -lm" 
setenv NETCDF_INC `nc-config --fflags`/netcdf

# Set list of core programs
set core  = "create_startf caltra trace select density lidar"

# Set list of goodies
set tools = "traj2num lsl2rdf changet extract getmima gettidiff getvars list2lsl lsl2list mergetra newtime reformat timeres trainfo difference datelist tracal" 

# Set list of libraries
set libs  = "iotra ioinp inter times libcdfio libcdfplus"

# Core programs
foreach prog ( $core )
   if ( "${prog}" == "${mode}" ) then
      set core  = ${prog}
      set mode  = "core"
   endif
end

# Goodies
foreach tool ( $tools )
   if ( "${tool}" == "${mode}" ) then
      set tools = ${tool}
      set mode  = "goodies"
   endif
end

# Libraries
foreach lib ( $libs )
   if ( "${lib}" == "${mode}" ) then
      set libs  = ${lib}
      set mode  = "lib"
   endif
end

# Check that the mode is ok 
if ( "${mode}" == "all"     ) goto modeok
if ( "${mode}" == "lib"     ) goto modeok
if ( "${mode}" == "core"    ) goto modeok
if ( "${mode}" == "goodies" ) goto modeok
if ( "${mode}" == "links"   ) goto modeok
if ( "${mode}" == "clean"   ) goto modeok
if ( "${mode}" == "docu"    ) goto modeok
if ( "${mode}" == "sync"    ) goto modeok
echo "Unsupported mode ${mode} ... Stop"
exit 1

modeok:

# -----------------------------------------------------------------------------
# Make clean 
# -----------------------------------------------------------------------------

if ( "${mode}" == "clean" ) then

cd ${LAGRANTO}/

foreach prog ( $core )
   \rm -f ${prog}/${prog} ${prog}/${prog}.o
end
\rm -f trace/calvar.o select/special.o

foreach tool ( $tools )
  \rm -f goodies/${tool} goodies/${tool}.o 
end
\rm -f goodies/*.mod goodies/*.o

\rm lib/*.a lib/*.o

foreach prog ( $core )
   \rm -f bin/${prog} bin/${prog}.sh  bin/${prog}.ecmwf
end
\rm -f bin/seltra bin/seltra.sh bin/seltra.ecmwf
foreach tool ( $tools )
  \rm -f bin/${tool} bin/${tool}.sh bin/${tool}.ecmwf
end
\rm -f bin/lagrantohelp.sh bin/lagrantohelp.ecmwf
\rm -f bin/startf bin/startf.sh bin/startf.ecmwf
\rm -f bin/lagranto.sh bin/lagranto.ecmwf

\rm ${LAGRANTO}/startf

exit 0

endif

# -----------------------------------------------------------------------------
# Install reference 
# -----------------------------------------------------------------------------

if ( ("${mode}" == "docu") | ("${mode}" == "all" ) ) then

echo "Installing documentation"
echo "-----------------------------------------------------------------"
echo

cd ${LAGRANTO}/docu/reference/

\rm -f reference.ps
\rm -f reference2.ps
groff -man ../man/*.0 > reference2.ps
ps2pdf reference2.ps

latex title
dvips title
ps2pdf title.ps

gs -dBATCH -dNOPAUSE -q -sDEVICE=pdfwrite -sOutputFile=reference.pdf title.pdf reference2.pdf

\rm -f *.aux *.log *.dvi
\rm -f title.ps  reference2.ps
\rm -f title.pdf reference2.pdf

endif

if ( "${mode}" == "docu" ) exit 0



# -----------------------------------------------------------------------------
# Install libraries 
# -----------------------------------------------------------------------------

if ( ("${mode}" == "lib") | ("${mode}" == "all" ) ) then

echo "-----------------------------------------------------------------"
echo "Installing libraries"
echo "-----------------------------------------------------------------"
echo

# Change to library directory
cd ${LAGRANTO}/lib

# Set the correct netCDF interface
echo
if ( "${ncdf}" == "ive" ) then
   echo " ioinp_ive.f -> ioinp.f"
   \cp ioinp_ive.f ioinp.f
endif
if  ( "${ncdf}" == "nc" ) then
   echo " ioinp_cdf.f -> ioinp.f"
   \cp ioinp_nc.f ioinp.f
endif
if  ( "${ncdf}" == "cdo" ) then
   echo " ioinp_cdo.f -> ioinp.f"
   \cp ioinp_cdo.f ioinp.f
endif
if  ( "${ncdf}" == "mch" ) then
   echo " ioinp_mch.f -> ioinp.f"
   \cp ioinp_mch.f ioinp.f
endif
echo

# Loop over all libraries - compile and make library
foreach lib ( $libs )

\rm -f ${lib}.a
\rm -f ${lib}.o
echo ${FORTRAN} -c -O ${lib}.f
${FORTRAN} -c -O ${NETCDF_INC} ${lib}.f
ar r ${lib}.a ${lib}.o
\rm -f ${lib}.l ${lib}.o
ranlib ${lib}.a
if ( ! -f ${lib}.a ) then
  echo "Problem in compiling ${lib} ... Stop"
  exit 1
endif

end

endif

if ( "${mode}" == "lib" ) exit 0

# -----------------------------------------------------------------------------
# Check that libraries are ok
# -----------------------------------------------------------------------------

echo
echo "-----------------------------------------------------------------"
echo "Check that all libraries are available"
echo "-----------------------------------------------------------------"
echo

# Change to library directory
cd ${LAGRANTO}/lib

# Check whether all libraries are available
foreach lib ( $libs )

if ( ! -f ${lib}.a ) then
  echo "Library ${lib} missing... Stop"
  exit 1
else
  ls -l ${lib}.a    
endif

end

# Exit if only libraries shoudl be installed
if ( "${mode}" == "lib" ) exit 0

# -----------------------------------------------------------------------------
# Compile Lagrango core programs
# -----------------------------------------------------------------------------

if ( ("${mode}" == "all" ) | ("${mode}" == "core" ) ) then

echo
echo "-----------------------------------------------------------------"
echo "Installing Lagranto core programs"
echo "-----------------------------------------------------------------"

foreach prog ( $core )

echo
echo "----- $prog"
echo
cd ${LAGRANTO}/${prog}
\rm -f ${prog}.o 
\rm -f ${prog}
if ( "${prog}" == "trace"  ) \rm calvar.o
if ( "${prog}" == "select" ) \rm special.o
\rm -f ${prog}
make -f ${prog}.make
if ( ! -f ${prog} ) then
  echo "Problem in compiling ${prog} ... Stop"
  exit 1
endif

end

endif

if ( "${mode}" == "core" ) exit 0

# -----------------------------------------------------------------------------
# Check that all Lagranto core programs are available
# -----------------------------------------------------------------------------

echo
echo "-----------------------------------------------------------------"
echo "Check that all Lagranto core programs are available"
echo "-----------------------------------------------------------------"
echo

foreach prog ( $core )

  cd ${LAGRANTO}/${prog}
  if ( ! -f ${prog} ) then
    echo "${prog} is missing... Stop"
    exit 1
  else
    ls -l ${prog}    
  endif

end

# Exit if only core programs shoudl be installed
if ( "${mode}" == "core" ) exit 0

# -----------------------------------------------------------------------------
# Compile Lagrango goodies
# -----------------------------------------------------------------------------

if ( ("${mode}" == "all" ) | ("${mode}" == "goodies" ) ) then

echo
echo "-----------------------------------------------------------------"
echo "Installing Lagranto goodies"
echo "-----------------------------------------------------------------"

cd ${LAGRANTO}/goodies

foreach tool ( $tools )

echo
echo "----- ${tool}"
echo
\rm -f ${tool}.o
\rm -f ${tool}
if ( -f ${tool}.make ) then
   make -f ${tool}.make
else if ( -f ${tool}.install ) then
   ./${tool}.install
endif

if ( ! -f ${tool} ) then
  echo "Problem in compiling ${tool} ... Stop"
  exit 1
endif

end

endif

if ( "${mode}" == "goodies" ) exit 0

# -----------------------------------------------------------------------------
# Check that all Lagranto goodies are available
# -----------------------------------------------------------------------------

echo
echo "-----------------------------------------------------------------"
echo "Check that all Lagranto goodies are available"
echo "-----------------------------------------------------------------"
echo

cd ${LAGRANTO}/goodies

foreach tool ( $tools )

if ( ! -f ${tool} ) then
  echo "${tool} is missing... Stop"
  exit 1
else
  ls -l ${tool} 
endif

end

endif

# Exit if only goodies should be installed
if ( "${mode}" == "goodies" ) exit 0

# -----------------------------------------------------------------------------
# Create links to programs
# -----------------------------------------------------------------------------

if ( ("${mode}" == "all" ) | ("${mode}" == "links" )  ) then

echo
echo "-----------------------------------------------------------------"
echo "Create links in ${LAGRANTO}/bin/"
echo "-----------------------------------------------------------------"
echo

if ( ! -d ${LAGRANTO}/bin ) mkdir ${LAGRANTO}/bin
cd ${LAGRANTO}/bin

ln -svf ${LAGRANTO}/bin/lagranto            lagranto.sh
ln -svf ${LAGRANTO}/bin/lagranto            lagranto.ecmwf  
ln -svf ${LAGRANTO}/bin/lagrantohelp        lagrantohelp.sh
ln -svf ${LAGRANTO}/bin/lagrantohelp        lagrantohelp.ecmwf

ln -svf ${LAGRANTO}/caltra/caltra.sh        caltra.sh
ln -svf ${LAGRANTO}/startf/create_startf.sh create_startf.sh
ln -svf ${LAGRANTO}/select/select.sh        select.sh
ln -svf ${LAGRANTO}/select/select.sh        seltra.sh
ln -svf ${LAGRANTO}/trace/trace.sh          trace.sh
ln -svf ${LAGRANTO}/density/density.sh      density.sh
ln -svf ${LAGRANTO}/startf/create_startf.sh startf.sh
ln -svf ${LAGRANTO}/lidar/lidar.sh          lidar.sh

ln -svf ${LAGRANTO}/caltra/caltra.sh        caltra.ecmwf
ln -svf ${LAGRANTO}/startf/create_startf.sh create_startf.ecmwf
ln -svf ${LAGRANTO}/select/select.sh        select.ecmwf
ln -svf ${LAGRANTO}/select/select.sh        seltra.ecmwf
ln -svf ${LAGRANTO}/trace/trace.sh          trace.ecmwf
ln -svf ${LAGRANTO}/density/density.sh      density.ecmwf
ln -svf ${LAGRANTO}/startf/create_startf.sh startf.ecmwf
ln -svf ${LAGRANTO}/lidar/lidar.sh          lidar.ecmwf

ln -svf ${LAGRANTO}/caltra/caltra.sh        caltra
ln -svf ${LAGRANTO}/startf/create_startf.sh create_startf
ln -svf ${LAGRANTO}/select/select.sh        select
ln -svf ${LAGRANTO}/select/select.sh        seltra
ln -svf ${LAGRANTO}/trace/trace.sh          trace
ln -svf ${LAGRANTO}/density/density.sh      density
ln -svf ${LAGRANTO}/startf/create_startf.sh startf
ln -svf ${LAGRANTO}/lidar/lidar.sh          lidar

foreach tool ( $tools )

ln -svf ${LAGRANTO}/goodies/${tool}.sh     ${tool}.sh 
ln -svf ${LAGRANTO}/goodies/${tool}.sh     ${tool}.ecmwf
ln -svf ${LAGRANTO}/goodies/${tool}.sh     ${tool} 

end

# Set link for create_startf / startf
\rm -f ${LAGRANTO}/startf
ln -svf ${LAGRANTO}/create_startf ${LAGRANTO}/startf

# Set extra name for <select> to avoid conflict in BASH
ln -svf ${LAGRANTO}/select/select.sh        seltra
ln -svf ${LAGRANTO}/select/select.sh        seltra.sh
ln -svf ${LAGRANTO}/select/select.sh        seltra.ecmwf

endif


# -----------------------------------------------------------------------------
# Synchronise ( development -> operational ) 
# -----------------------------------------------------------------------------

if ( ("${mode}" == "all" ) | ("${mode}" == "sync" )  ) then

echo
echo "-----------------------------------------------------------------"
echo "Sync ( lagranto.ecmwf -> lagranto )"
echo "-----------------------------------------------------------------"
echo

cd ${path_sync}/bin/

ln -svf ${path_devel}/bin/lagranto.sh            lagranto.ecmwf
ln -svf ${path_devel}/bin/lagrantohelp.sh        lagrantohelp.ecmwf

ln -svf ${path_devel}/caltra/caltra.sh           caltra.ecmwf
ln -svf ${path_devel}/startf/create_startf.sh    create_startf.ecmwf
ln -svf ${path_devel}/select/select.sh           select.ecmwf
ln -svf ${path_devel}/trace/trace.sh             trace.ecmwf
ln -svf ${path_devel}/density/density.sh         density.ecmwf
ln -svf ${path_devel}/startf/create_startf.sh    startf.ecmwf
ln -svf ${path_devel}/lidar/lidar.sh             lidar.ecmwf
ln -svf ${path_devel}/lidar/seltra.sh            seltra.ecmwf

foreach tool ( $tools )

ln -svf ${path_devel}/goodies/${tool}.sh         ${tool}.ecmwf

end
 
# Set all permissions
chmod -R og+rx ${path_sync}/bin/

endif

# -----------------------------------------------------------------------------
# Final tasks
# -----------------------------------------------------------------------------

echo
echo "-----------------------------------------------------------------"
echo "Installation complete"
echo "-----------------------------------------------------------------"
echo
echo "Please set the environmental variable LAGRANTO"
echo
echo "    setenv LAGRANTO ${LAGRANTO}"
echo
 




