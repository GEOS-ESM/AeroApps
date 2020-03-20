#!/bin/csh

# Write usage information
if ( ${#argv} == 0) then
  echo 
  ${LAGRANTO}/bin/lagrantohelp gettidiff short
  echo  
  exit 0
endif

${LAGRANTO}/goodies/gettidiff $*

exit 0

