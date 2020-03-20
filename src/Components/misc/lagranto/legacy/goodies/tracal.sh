#!/bin/csh

# Write usage information
if ( ${#argv} == 0) then
  echo 
  ${LAGRANTO}/bin/lagrantohelp tracal short
  echo  
  exit 0
endif

set inpfile = $1
set outfile = $2
set expr    = $3

set dim=`${LAGRANTO}/goodies/trainfo.sh  ${inpfile} dim` 

\rm -f tracal.param
echo \"${inpfile}\"  >! tracal.param
echo \"${outfile}\"  >> tracal.param
echo \"${expr}\"     >> tracal.param
echo ${dim}          >> tracal.param

${LAGRANTO}/goodies/tracal

\rm -f tracal.param

exit 0

