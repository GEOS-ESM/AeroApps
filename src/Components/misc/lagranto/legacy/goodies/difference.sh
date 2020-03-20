#!/bin/csh


# Write usage information
if ( ${#argv} == 0) then
  echo 
  ${LAGRANTO}/bin/lagrantohelp difference short
  echo  
  exit 0
endif

# Get input and output trajectory file
set inpfile1  = $1
set inpfile2  = $2
set outfile   = $3
set fieldname = $4

# Handle optional arguments

set mode = "single"

while ( $#argv > 0 )

  switch ( $argv[1] )

   case -single
     set mode  = "single"
   breaksw

   case -max
     set mode  = "max"
   breaksw

  endsw

  shift;

end

# Get the dimensions of the trajectory files
set dim1=`${LAGRANTO}/goodies/trainfo.sh ${inpfile1} dim` 
set dim2=`${LAGRANTO}/goodies/trainfo.sh ${inpfile2} dim` 

# Prepare parameter file and run program
\rm -f difference.param
echo \"${inpfile1}\"  >! difference.param
echo \"${inpfile2}\"  >> difference.param
echo \"${outfile}\"   >> difference.param
echo ${dim1}          >> difference.param
echo ${dim2}          >> difference.param
echo \"${mode}\"      >> difference.param
echo \"${fieldname}\" >> difference.param

${LAGRANTO}/goodies/difference

# Make clean
#\rm -f difference.param

exit 0

