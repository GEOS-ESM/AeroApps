#!/bin/csh

# ---------------------------------------------------------------------
# Usage, parameter settings
# ---------------------------------------------------------------------

# Write usage information
if ( (${#argv} == 0) | (${#argv} < 2) ) then
  echo 
  ${LAGRANTO}/bin/lagrantohelp density short
  echo
  exit 0
endif

# Write title
echo 
echo '========================================================='
echo '       *** START OF PREPROCESSOR DENSITY ***              '
echo

# ---------------------------------------------------------------------
# Handle input parameters
# ---------------------------------------------------------------------

# Get input parameters
set inpfile = $1
set outfile = $2

# Set default values
set radius     = 100
set unit       = 'km'
set grid       = "360 180 -180. -90. 1. 1."
set mode       = 'keep'
set tratime    = 'all'
set param      = 0
set gridtype   = 'latlon'
set sel_file   = 'nil'
set sel_format = 'nil'
set field      = 'nil'    
set crefile    = 0

# Handle optional arguments
while ( $#argv > 0 )

  switch ( $argv[1] )

   case -radius
     set radius = $argv[2]
     set unit   = $argv[3]
     echo "Flag '-radius'     -> ${radius} ${unit} (user defined)"
     shift;
     shift;
   breaksw

   case -time
     set tratime=$argv[2]
     echo "Flag '-tratime'    -> ${tratime} (user defined)"
     shift;
   breaksw

   case -interp
     set param=$argv[2]
     
     if ( "$argv[3]" == "h"   ) set mode = "time"
     if ( "$argv[3]" == "km"  ) set mode = "space"
     if ( "$argv[3]" == "deg" ) set mode = "grid"
     
     echo "Flag '-interp'      -> ${mode} ${param} (user defined)"
     shift;
     shift;
   breaksw

   case -create
     set crefile = 1
     echo "Flag '-create'      -> true (user defined)"
   breaksw

   case -latlon
     set gridtype = 'latlon'
     if ( "$argv[2]" == "dynamic" ) then
	    set grid = "0 0 0 0 0 0"
	    echo "Flag '-latlon'       -> dynamic  (user defined)"
	    shift;
     else
        set nlon   = $argv[2]
        set nlat   = $argv[3]
        set lonmin = $argv[4]
        set latmin = $argv[5]
        set dlon   = $argv[6]
        set dlat   = $argv[7]
        set grid   = "${nlon} ${nlat} ${lonmin} ${latmin} ${dlon} ${dlat}"
        echo "Flag '-latlon      ->  ${grid}  (user defined)"
        shift;
        shift;
        shift;
        shift;
        shift;
        shift;
     endif
   breaksw

   case -rotated
     set gridtype = 'rotated'
     set clon     =  $argv[2]
     set clat     =  $argv[3]
     set nlonlat = $argv[4]
     set dlonlat = $argv[5]
     set grid   = "${clon} ${clat} ${nlonlat} ${dlonlat}"
     echo "Flag '-rotated     ->  ${clon}, ${clat}, ${nlonlat}, ${dlonlat}  (user defined)"
     shift;
     shift;
     shift;
     shift;
    breaksw

    case -index
       set sel_file   = $argv[2]
       set sel_format = 'index'
       echo "Flag '-index'     -> ${sel_file} (user defined)"
       shift;
   breaksw

   case -boolean
       set sel_file   = $argv[2]
       set sel_format = 'boolean'
       echo "Flag '-boolean'   -> ${sel_file} (user defined)"
       shift;
   breaksw

   case -field
       set field   = $argv[2]
       echo "Flag '-field'   -> ${field} (user defined)"
       shift;
   breaksw

  endsw
 
  shift;

end

# ---------------------------------------------------------------------
# Do some checks and preparation, then run the program
# ---------------------------------------------------------------------

# Rename field <time> to <TIME> to avoid conflict with time coordinate
# on the netCDF file
if ( "${field}" == "time" ) set field = "TIME"    

# Determine the time step 
if ( "${tratime}" == "all" ) then
    set step = 0
else
    set timelist = (`${LAGRANTO}/bin/trainfo.sh $inpfile times`) 
    set step     = 0
    set found    = 0
    foreach val ( ${timelist} )
       @ step = ${step} + 1
       if ( "${tratime}" == "${val}" ) then
          set found   = ${step}
       endif
     end
     if ( ${found} == 0 ) then
       echo "Invalid time ${tratime} for gridding"
       echo "${timelist}"
       exit 1
      endif
      set step = ${found}
endif

# Check consistency of arguments
if ( ( "${mode}" == "time" ) & ( ${step} != 0 ) ) then
    echo " ERROR: Options 'interp -time' and 'step' incompatible' "
    exit 1
endif
if ( ( "${mode}" == "space" ) & ( ${step} != 0 ) ) then
    echo " ERROR: Options 'interp -space' and 'step' incompatible' "
    exit 1
endif

# Get trajectory info
set ntra   = (`${LAGRANTO}/bin/trainfo.sh $inpfile ntra`)
set ntime  = (`${LAGRANTO}/bin/trainfo.sh $inpfile ntim`)
set nfield = (`${LAGRANTO}/bin/trainfo.sh $inpfile ncol`)

# Check whether selection file is available
if ( "${sel_file}" != "nil" ) then
   if ( ! -f ${sel_file} ) then
     echo " ERROR: selection file ${sel_file} is missing... Stop"
     exit 1
   endif
endif

# Check whether output file exists - set the <crefile> flag
if ( "${crefile}" == "0" ) then
   if ( ! -f ${outfile} ) then
       set crefile = 1
       echo
       echo "${outfile} will be created... "
   else
       echo
       echo "${outfile} will be modified ..."
   endif
else
   echo
   echo "${outfile} will be created... "
endif

# Chewck whether the variable exists - set the <crevar> flag
if ( "${crefile}" == "0" ) then    
    set varlist = ` ${LAGRANTO}/goodies/getvars ${outfile}` 
    set crevar  = 1
    foreach var ( ${varlist} )
       if ( "${var}" == "${field}" ) set crevar = 0
    end
else
   set crevar = 1
endif

# Prepare parameter file and run program
\rm -f density.param
touch density.param
echo ${inpfile}                 >> density.param
echo ${outfile}                 >> density.param
echo \"${field}\"               >> density.param
echo ${ntime} ${nfield} ${ntra} >> density.param
echo ${gridtype}                >> density.param
echo ${grid}                    >> density.param
echo ${radius} ${unit}          >> density.param
echo ${mode}                    >> density.param
echo ${param}                   >> density.param
echo ${step}                    >> density.param
echo \"${sel_file}\"            >> density.param
echo \"${sel_format}\"          >> density.param
echo ${crefile}                 >> density.param
echo ${crevar}                  >> density.param

# Write status info
echo 
echo '       *** END OF PREPROCESSOR DENSITY ***              '
echo '========================================================='
echo

# Run density
${LAGRANTO}/density/density

# Make clean
\rm -f density.param

exit 0




