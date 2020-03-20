#!/bin/csh

# -----------------------------------------------------------------------------
# Set some parameters
# -----------------------------------------------------------------------------

# Set input file
set inpfile=$1
if ( ${#argv} == 2 ) then
  set mode=$2
else
  set mode='all'
endif

# Write usage information
if ( ${#argv} == 0) then
  echo 
  ${LAGRANTO}/bin/lagrantohelp trainfo short
  echo  
  exit 0
endif

# Set Fortran program
set prog=${LAGRANTO}/goodies/trainfo

# -----------------------------------------------------------------------------
# Run program 
# -----------------------------------------------------------------------------

\rm -f trainfo.param
echo \"${inpfile}\" >! trainfo.param 
echo \"${mode}\"    >> trainfo.param

${prog}

\rm -f trainfo.param

exit 0

# -----------------------------------------------------------------------------
# Old code: shell script extraction of ntra,ntim,ncol
# -----------------------------------------------------------------------------

# Get line numbers of first trajectory block (separated by empty line) 
set first=4
loop1:
  @ first = ${first} + 1
  set line=`sed -ne "${first},${first}p" ${inpfile}` 
if ( "${line}" != "" ) goto loop1
@ final = ${first} + 1
loop2:
  @ final = ${final} + 1
  set line=`sed -ne "${final},${final}p" ${inpfile}`
if ( "${line}" != "" ) goto loop2
@ first = ${first} + 1
@ final = ${final} - 1

# Set the number of fields, of times and of trajectories
set ntime=`echo "1 + ${final} - ${first}" | bc`  
set line=`sed -ne "${first},${first}p" ${inpfile}`
set ncol=`echo ${line} | awk '{print NF}'` 
set nlines=`wc -l ${inpfile} | awk '{print $1}'`
set ntra=`echo "(${nlines}-4)/(${ntime}+1)" | bc`

# Write info
echo ${ntra} ${ntime} ${ncol}
