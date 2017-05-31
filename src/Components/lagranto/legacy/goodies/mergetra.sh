#!/bin/csh

# -----------------------------------------------------------------------------
# Set some parameters
# -----------------------------------------------------------------------------

# Write usage information
if ( ${#argv} == 0) then
  echo 
  ${LAGRANTO}/bin/lagrantohelp mergetra short 
  echo  
  exit 0
endif

# Set input file
set inpfile1 = $1
set inpfile2 = $2
set outfile  = $3
shift; shift; shift

# Handle optional arguments 
set datecheck = "datecheck"
while ( $#argv > 0 )
   if ( "$1" == "-nodatecheck" ) set datecheck="nodatecheck"
   shift
end

# Set Fortran program
set prog=${LAGRANTO}/goodies/mergetra

# -----------------------------------------------------------------------------
# Run program 
# -----------------------------------------------------------------------------

set dims1=`${LAGRANTO}/goodies/trainfo.sh ${inpfile1} dim` 
set dims2=`${LAGRANTO}/goodies/trainfo.sh ${inpfile2} dim` 

\rm -f mergetra.param
echo \"${inpfile1}\"  >! mergetra.param 
echo \"${inpfile2}\"  >> mergetra.param
echo \"${outfile}\"   >> mergetra.param
echo ${dims1}         >> mergetra.param
echo ${dims2}         >> mergetra.param
echo \"${datecheck}\" >> mergetra.param

${prog}

exit 0

