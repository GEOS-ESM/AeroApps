#!/bin/csh

# Set some parameters
set crifile=${PWD}/select.parsed

# Set base directories (run+prog)
set cdfdir=${PWD}
set tradir=${PWD}

# Write usage information
if ( ${#argv} ==  0 ) then
  echo 
  ${LAGRANTO}/bin/lagrantohelp select short
  echo  
  exit 0
endif

# Write list of special selection criteria
if ( "$1" == "-special" ) then
    grep '%)' ${LAGRANTO}/select/special.f
    exit 0
endif

# Write title
echo 
echo '========================================================='
echo '       *** START OF PREPROCESSOR SELECT ***              '
echo

# Save input arguments
set inpfile=$1
set outfile=$2
set crit="$3"
shift
shift
shift

# Handle optional arguments
echo
echo '---- OPTIONAL FLAGS -------------------------------------'
echo  

set noclean = 'false'
set format  = 'trajectory'
set regionf = 'regionf'
set trigger = 'nil'

while ( $#argv > 0 )

  switch ( $argv[1] )

   case -noclean
     set noclean = 'true'
     echo "noclean       -> true (user defined)"
   breaksw

   case -trigger
     set trigger = '-trigger'
     echo "trigger      ->  true (user defined)"
   breaksw

   case -boolean
     set format = 'boolean'
     echo "format       ->  boolean (user defined)"
   breaksw

   case -index
     set format = 'index'
     echo "format       ->  index (user defined)"
   breaksw

   case -count
     set format = 'count'
     echo "format       ->  count (user defined)"
   breaksw

   case -startf
     set format = 'startf'
     echo "format       ->  startf (user defined)"
   breaksw

   case -regionf
     set regionf = $argv[2]
     echo "regionf                -> ${regionf} (user defined)"
     shift;
   breaksw

  endsw
 
  shift;

end

# Decide whether <select> is a file or an explicit criterion
set flag_select = 'criterion'
set test = `echo ${crit} | grep ':' | wc -c`
if ( "${test}" == "0" ) then 
  set flag_select     = 'file'
  set flag_selectfile = $crit
  if ( -f ${flag_selectfile} ) then 
     set crit             = `cat ${flag_selectfile}`
  else
     echo " ERROR: criterion file ${flag_selectfile} is missing... Stop"
     exit 1
  endif
endif

# Get the start, end and reference date for the tracing
set ntra      =  `${LAGRANTO}/goodies/trainfo.sh ${inpfile} ntra`
set ntim      =  `${LAGRANTO}/goodies/trainfo.sh ${inpfile} ntim`
set ncol      =  `${LAGRANTO}/goodies/trainfo.sh ${inpfile} ncol`
set times     =  `${LAGRANTO}/goodies/trainfo.sh ${inpfile} times`
set vars      =  `${LAGRANTO}/goodies/trainfo.sh ${inpfile} vars`

# Split the criterion into subunits (with Perl)
${LAGRANTO}/select/select.perl "${crit}" >! ${crifile} 
if ( "${status}" != "0" ) then
     echo "ERROR:  Parser <select> failed"
     exit 1
endif

# Write some status information
echo
echo '---- DIRECTORIES AND PROGRAMS ---------------------------'
echo    
echo "CDF directory         : ${cdfdir}"
echo "TRA directory         : ${tradir}"
echo "PROGRAM SELECT        : ${LAGRANTO}/select/select"
echo "PARSER                : ${LAGRANTO}/select/select.perl"
echo "CRITERION FILE        : ${crifile}"
echo
echo '---- INPUT PARAMETERS -----------------------------------'
echo    
echo "Input file            : ${inpfile}"
echo "Output file           : ${outfile}"
echo "Output format         : ${format}"
if ( "${flag_select}" == "criterion" ) then
    echo "Criterion             : ${crit}"
else
    echo "Criterion             : ${crit} (from file ${flag_selectfile})"
endif
echo
echo '---- INPUT FILE -----------------------------------------'
echo    
echo "# TRA                 : ${ntra}"
echo "# TIMES               : ${ntim}"
echo "# COLUMNS             : ${ncol}"
echo "Times                 : ${times}"
echo "Variables             : ${vars}"
echo
echo '---- PARSED CRITERION ------------------------------------'
echo 
cat ${crifile}

# Finish the preprocessor
echo 
echo '       *** END OF PREPROCESSOR SELECT ***              '
echo '========================================================='
echo

# Run the selection programme
echo \"${inpfile}\" >! select.param
echo \"${outfile}\" >> select.param
echo \"${format}\"  >> select.param
echo \"${crifile}\" >> select.param
echo ${ntra}        >> select.param
echo ${ntim}        >> select.param
echo ${ncol}        >> select.param
echo \"${regionf}\" >> select.param
echo \"${trigger}\" >> select.param

${LAGRANTO}/select/select

if ( "${status}" != "0" ) then
  echo "ERROR:  Program <select> failed"
  exit 1
endif

# Make clean
if ( "${noclean}" == "false" ) then
  \rm -f select.param
  \rm -f ${crifile}
endif

exit 0
