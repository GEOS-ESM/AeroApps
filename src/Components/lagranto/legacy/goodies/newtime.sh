#!/bin/csh

# Write usage information
if ( ${#argv} == 0) then
  echo 
  ${LAGRANTO}/bin/lagrantohelp newtime short
  echo  
  exit 0
endif

${LAGRANTO}/goodies/newtime $*

exit 0

