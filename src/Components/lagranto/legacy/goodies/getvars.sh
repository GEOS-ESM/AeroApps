#!/bin/csh

# Write usage information
if ( ${#argv} == 0) then
  echo 
  ${LAGRANTO}/bin/lagrantohelp getvars short
  echo  
  exit 0
endif

${LAGRANTO}/goodies/getvars $*

exit 0

