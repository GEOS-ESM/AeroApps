#!/bin/csh


# Write usage information
if ( ${#argv} == 0) then
  echo 
  ${LAGRANTO}/bin/lagrantohelp reformat short
  echo  
  exit 0
endif

set inpfile=$1
set outfile=$2

set prog1=${LAGRANTO}/goodies/reformat
set prog2=${LAGRANTO}/goodies/trainfo.sh

set dim=`${prog2} ${inpfile} dim` 

\rm -f reformat.param
echo \"${inpfile}\"  >! reformat.param
echo \"${outfile}\"  >> reformat.param
echo ${dim}          >> reformat.param

${prog1}

\rm -f reformat.param

exit 0

