#!/bin/csh


# Write usage information
if ( ${#argv} == 0) then
  echo 
  ${LAGRANTO}/bin/lagrantohelp timeres short
  echo  
  exit 0
endif

# Get input and output trajectory file
set inpfile = $1
set outfile = $2

# Handle optional arguments

set mode = "cubic"
set shift = 0
set value = 0
set unit  = "min"

while ( $#argv > 0 )

  switch ( $argv[1] )

   case -h
     set unit  = "h"
     set value = $argv[2]
     shift;
   breaksw

   case -min
     set unit  = "min"
     set value = $argv[2]
     shift;
   breaksw

   case -linear
     set mode  = "linear"
   breaksw

   case -cubic
     set mode  = "cubic"
   breaksw
   
   case -shifth
     set shift = $argv[2]
     shift;
   breaksw

  endsw

  shift;

end

# Get the dimensions of the trajectory file
set dim=`${LAGRANTO}/goodies/trainfo.sh ${inpfile} dim` 

# Prepare parameter file and run program
\rm -f timeres.param
echo \"${inpfile}\"  >! timeres.param
echo \"${outfile}\"  >> timeres.param
echo ${dim}          >> timeres.param
echo ${value}        >> timeres.param
echo \"${unit}\"     >> timeres.param
echo \"${mode}\"     >> timeres.param
echo ${shift}        >> timeres.param

${LAGRANTO}/goodies/timeres

# Make clean
#\rm -f timeres.param

exit 0

