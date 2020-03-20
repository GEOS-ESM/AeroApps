#!/bin/csh


# Write usage information
if ( ${#argv} == 0) then
  echo 
  ${LAGRANTO}/bin/lagrantohelp lsl2list short
  echo  
  exit 0
endif

set inpfile=$1
set outfile=$2

set prog1=${LAGRANTO}/goodies/lsl2list

set dim=`${LAGRANTO}/goodies/trainfo.sh ${inpfile} dim` 

\rm -f lsl2list.param
echo \"${inpfile}\"       >! lsl2list.param
echo \"${outfile}\"       >> lsl2list.param
echo ${dim}               >> lsl2list.param

${prog1}

#\rm -f lsl2list.param

exit 0

