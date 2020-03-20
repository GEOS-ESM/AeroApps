#!/bin/csh

# Write usage information
if ( ${#argv} == 0) then
  echo 
  ${LAGRANTO}/bin/lagrantohelp extract short
  echo  
  exit 0
endif

# Set the input and output filenames
set inpfile = $1
set outfile = $2
shift; shift;

# Handle optional arguments (get a list, remove , and / separators)
set mode=$1
shift;
set str=
while ( $#argv > 0 )
   set str = "${str} $1"
   shift
end

set str=`echo $str | sed -e "s/,//g"`  
set str=`echo $str | sed -e "s/\///g"`  

# Check that the mode is valid
if ( "${mode}" == "-var"      )  goto next
if ( "${mode}" == "-time"     )  goto next
if ( "${mode}" == "-tra"      )  goto next
if ( "${mode}" == "-startf"   )  goto next
if ( "${mode}" == "-index"    )  goto next
if ( "${mode}" == "-boolean"  )  goto next
if ( "${mode}" == "-pattern"  )  goto next
if ( "${mode}" == "-leaving"  )  goto next
if ( "${mode}" == "-staying"  )  goto next

echo " Invalid mode ${mode}..."
exit 1

next:

# Set program names
set prog1=${LAGRANTO}/goodies/extract
set prog2=${LAGRANTO}/goodies/trainfo.sh

# Get trajectory dimensions
set dim=`${prog2} ${inpfile} dim` 

# Run program
\rm -f extract.param
echo \"${inpfile}\"  >! extract.param
echo \"${outfile}\"  >> extract.param
echo \"${mode}\"     >> extract.param
echo ${dim}          >> extract.param
echo \"${str}\"      >> extract.param

${prog1}

#\rm -f extract.param

exit 0

