#!/bin/csh


# Write usage information
if ( ${#argv} == 0) then
  echo 
  ${LAGRANTO}/bin/lagrantohelp traj2num short
  echo  
  exit 0
endif

set inpfile = $1
set outfile = $2
set mode    = $3

set dim=`${LAGRANTO}/goodies/trainfo.sh ${inpfile} dim` 

\rm -f traj2num.param
echo \"${inpfile}\"  >! traj2num.param
echo \"${outfile}\"  >> traj2num.param
echo ${dim}          >> traj2num.param
echo \"${mode}\"     >> traj2num.param

${LAGRANTO}/goodies/traj2num

\rm -f traj2num.param

exit 0

