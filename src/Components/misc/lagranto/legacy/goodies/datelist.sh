#!/bin/csh

# Write usage information
if ( ${#argv} == 0) then
  echo 
  ${LAGRANTO}/bin/lagrantohelp datelist short
  echo  
  exit 0
endif

# Handle fixed arguments
set filename = $1
set mode     = $2

# Redirect output to screen if requested
if ( "${filename}" == "stdout" ) set filename="/dev/stdout"
if ( "${filename}" == "screen" ) set filename="/dev/stdout"

# Handle optional arguments
set interval = "6"

while ( $#argv > 0 )

  switch ( $argv[1] )

   case -create
     set startdate = $argv[2] 
     set finaldate = $argv[3] 
     shift;
     shift;
   breaksw

   case -indir
     set dirname = $argv[2] 
     shift;
   breaksw

   case -next
     set date = $argv[2] 
     shift;
   breaksw

   case -prev
     set date = $argv[2] 
     shift;
   breaksw

   case -isin
     set date = $argv[2] 
     shift;
   breaksw

   case -interval
     set interval = $argv[2] 
     shift;
   breaksw

   case -overlap
     set file1 = $argv[2] 
     set file2 = $argv[3] 
     shift;
     shift;
   breaksw

   case -onlyin1
     set file1 = $argv[2] 
     set file2 = $argv[3] 
     shift;
     shift;
   breaksw

   case -onlyin2
     set file1 = $argv[2] 
     set file2 = $argv[3] 
     shift;
     shift;
   breaksw

   case -totime
     set refdate = $argv[2] 
     shift;
   breaksw

   case -todate
     set refdate = $argv[2] 
     shift;
   breaksw

  endsw

  shift;

end

# Mode: -create startdate finaldate
if ( "${mode}" == "-create" ) then
    
   \rm -f datelist.param
   echo \"${filename}\"    >  datelist.param
   echo \"${mode}\"        >> datelist.param
   echo \"${startdate}00\" >> datelist.param
   echo \"${finaldate}00\" >> datelist.param
   echo ${interval}        >> datelist.param
   
   ${LAGRANTO}/goodies/datelist

endif

# Mode: -first 
if ( "${mode}" == "-first" ) then
    
   head -1 ${filename} 

endif

# Mode: -last 
if ( "${mode}" == "-last" ) then
    
   tail -1 ${filename} 

endif

# Mode: -ndates 
if ( "${mode}" == "-ndates" ) then
    
   wc -l ${filename} | awk '{ print $1}'

endif

# Mode: -timerange 
if ( "${mode}" == "-timerange" ) then
  
   set firstdate = `head -1 ${filename}`
   set finaldate = `tail -1 ${filename}`  

   ${LAGRANTO}/goodies/gettidiff ${finaldate} ${firstdate}

endif

# Mode: -indir {directory name}
if ( "${mode}" == "-indir" ) then

    ls -1 ${dirname} | perl -ne 'print  if s/.*([0-9]{8}_[0-9]{2}).*/\1/' | sort | uniq >! ${filename}

endif

# Mode: -next
if ( "${mode}" == "-next" ) then

  set last = `tail -1 ${filename}` 
  if ( "${date}" != "${last}" ) then
     set next = `sed -n "/${date}/{n;p;}" ${filename}`
     echo ${next}
  else
     echo "nil"
  endif

endif

# Mode: -prev
if ( "${mode}" == "-prev" ) then

  set first = `head -1 ${filename}` 
  if ( "${date}" != "${first}" ) then
     set prev = `sed -n "/${date}/{g;p;};h" ${filename}`
     echo ${prev}
  else
     echo "nil"
  endif

endif 

# Mode: -isin
if ( "${mode}" == "-isin" ) then

  set flag = `sed -n "/${date}/p" ${filename}`
    
  if ( "${flag}" != "" ) then
    echo "1"
  else
    echo "0"
  endif

endif 

# Mode: -overlap
if ( "${mode}" == "-overlap" ) then
    if ( "${filename}" != "/dev/stdout" ) then 
       set outfile = ${filename}.$$
       \rm -f ${outfile}
       grep -f ${file1} ${file2} > ${outfile}
       \mv ${outfile} ${filename}
    else
       grep -f ${file1} ${file2}
    endif
endif

# Mode: -onlyin1
if ( "${mode}" == "-onlyin1" ) then
    
   \rm -f datelist.param
   echo \"${filename}\"    >  datelist.param
   echo \"${mode}\"        >> datelist.param
   echo \"${file1}\"       >> datelist.param
   echo \"${file2}\"       >> datelist.param

   ${LAGRANTO}/goodies/datelist
    
endif

# Mode: -onlyin2
if ( "${mode}" == "-onlyin2" ) then
    set outfile1 = "tmp1.$$"
    set outfile2 = "tmp2.$$"
    \rm -f ${outfile1}
    \rm -f ${outfile2}
    grep -f  ${file2}    ${file1} >! ${outfile1}
    grep -vf ${outfile1} ${file2} >! ${outfile2}
    if ( "${filename}" != "/dev/stdout" ) then
      \mv ${outfile2} ${filename}
    else
      cat  ${outfile2}
    endif
    \rm -f ${outfile1}
    \rm -f ${outfile2}
endif

# Mode: -totime 
if ( "${mode}" == "-totime" ) then
       
   \rm -f datelist.param
   echo \"${filename}\"    >  datelist.param
   echo \"${mode}\"        >> datelist.param
   echo \"${refdate}00\"   >> datelist.param

   ${LAGRANTO}/goodies/datelist

endif

# Mode: -todate
if ( "${mode}" == "-todate" ) then
       
   \rm -f datelist.param
   echo \"${filename}\"    >  datelist.param
   echo \"${mode}\"        >> datelist.param
   echo \"${refdate}00\"   >> datelist.param

   ${LAGRANTO}/goodies/datelist

endif

# Mode: -randsample
if ( "${mode}" == "-randsample" ) then
       
   \rm -f datelist.param
   echo \"${filename}\"    >  datelist.param
   echo \"${mode}\"        >> datelist.param
   echo ${randsample}      >> datelist.param

   ${LAGRANTO}/goodies/datelist

endif

exit 0

