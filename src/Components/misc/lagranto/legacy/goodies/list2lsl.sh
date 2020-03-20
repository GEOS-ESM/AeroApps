#!/bin/csh


# Write usage information
if ( ${#argv} == 0) then
 echo 
  ${LAGRANTO}/bin/lagrantohelp list2lsl short
  echo  
  exit 0
endif

set inpfile=$1
set outfile=$2

set refdate   = `date +'%Y%m%d_%H%M'`
set timevalue = 0

while ( $#argv > 0 )

   switch ( $argv[1] )

   case -ref
     set refdate=$argv[2]
     shift;
   breaksw

   case -time
     set timevalue=$argv[2]
     shift;
   breaksw

  endsw
 
  shift;

end

set ntra=`wc -l ${inpfile} | awk '{print $1}'`

# Split the reference date
set yyyy=`echo ${refdate}   | cut -c 1-4`
set   mm=`echo ${refdate}   | cut -c 5-6`
set   dd=`echo ${refdate}   | cut -c 7-8`
set   hh=`echo ${refdate}   | cut -c 10-11`
set  min=`echo ${refdate}00 | cut -c 12-13`

\rm -f list2lsl.param
echo \"${inpfile}\"                      >! list2lsl.param
echo \"${outfile}\"                      >> list2lsl.param
echo ${ntra}                             >> list2lsl.param
echo ${yyyy} ${mm} ${dd} ${hh} ${min} 00 >> list2lsl.param 
echo ${timevalue}                        >> list2lsl.param 

${LAGRANTO}/goodies/list2lsl

\rm -f list2lsl.param

exit 0

