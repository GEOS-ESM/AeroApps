#!/bin/csh


# Write usage information
if ( ${#argv} == 0) then
  echo 
  ${LAGRANTO}/bin/lagrantohelp lsl2rdf short
  echo  
  exit 0
endif

set inpfile=$1
set outfile=$2

set prog1=${LAGRANTO}/goodies/lsl2rdf

set dim=`${LAGRANTO}/goodies/trainfo.sh ${inpfile} dim` 

\rm -f lsl2rdf.param
echo \"${inpfile}\"       >! lsl2rdf.param
echo \"${outfile}\"       >> lsl2rdf.param
echo ${dim}               >> lsl2rdf.param

${prog1}

\rm -f lsl2rdf.param

exit 0

