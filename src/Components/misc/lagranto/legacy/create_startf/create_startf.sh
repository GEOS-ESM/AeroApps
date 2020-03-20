#!/bin/csh -f

# -----------------------------------------------------------------------------------------------------------------
# Usage and parameter handling
# -----------------------------------------------------------------------------------------------------------------

# Write usage information
if ( ${#argv} == 0 ) then
  echo 
  ${LAGRANTO}/bin/lagrantohelp create_startf short
  echo
  exit 1
endif

# -----------------------------------------------------------------------------------------------------------------
# Handle input parameters
# -----------------------------------------------------------------------------------------------------------------

# Prefix of primary and secondary file; and set the filenames
set charp='P'
set chars='S'

# Write title
echo 
echo '========================================================='
echo '       *** START OF PREPROCESSOR CREATE_STARTF ***              '
echo

# Extract arguments and split reference date
set refdate   = $1
set ofile     = $2
set specifier = `echo $3` 
shift
shift
shift

# Check whether specifier is a file or whether it is explicitly written
set flag_criterion = 'file'
set test = `echo ${specifier} | grep '@' | wc -c`
if ( "${test}" != "0" ) then 
  set flag_criterion = 'criterion'
endif

# Get the criterion from the file
if ( "${flag_criterion}" == "file" ) then
    if ( -f ${specifier} ) then
       set filename  = ${specifier}
       set specifier = `cat ${specifier}`
    else
       echo " ERROR: cannot read criterion from file ${specifier}... Stop"
       exit 1
    endif
endif

echo "---- INPUT PARAMETERS ----------------------------------"
echo 
echo "Reference date        : ${refdate}"
if ( "${flag_criterion}" == "criterion" ) then 
   echo "Specifier             : ${specifier}"
else
   echo "Specifier             : ${specifier} [from file ${filename}]"
endif
echo "Output file           : ${ofile}"
echo

# Handle optional arguments
set tvfile    = 'tracevars'
set changet   = 'false'
set noclean   = 'false'
set regionf   = 'regionf'
set timecheck = 'no' 

while ( $#argv > 0 )

  switch ( $argv[1] )

   case -t
     set tvfile = $argv[2]
     echo "tvfile                -> ${tvfile} (user defined)"
     shift;
   breaksw

   case -changet
     set changet = 'true'
     echo "changet               -> true (user defined)"
   breaksw

   case -noclean
     set noclean = 'true'
     echo "noclean               -> true (user defined)"
   breaksw

   case -timecheck
     set timecheck = 'yes'
     echo "timecheck               -> yes (user defined)"
   breaksw

   case -regionf
     set regionf = $argv[2]
     echo "regionf                -> ${regionf} (user defined)"
     shift;
   breaksw

   endsw
 
   shift;

end

# No change of times necessary if no check requested
if ( "${timecheck}" == "no" ) then
   set  changet = 'false'
endif

# Split the reference date
set yyyy=`echo ${refdate}   | cut -c 1-4` 
set   mm=`echo ${refdate}   | cut -c 5-6` 
set   dd=`echo ${refdate}   | cut -c 7-8` 
set   hh=`echo ${refdate}   | cut -c 10-11` 
set  min=`echo ${refdate}00 | cut -c 12-13` 

# Set base directories (run+prog)
set tradir=${PWD}

# Set program paths and filenames 
set parfile=${tradir}/create_startf.param 
set crifile=${tradir}/create_startf.criterion

# Decide whether a tracing and selection is necessary
# If so, some intermediate file are writen (in Fortran binary)
set flag = `${LAGRANTO}/startf/create_startf.perl "${specifier}" | tail -1` 
if ( "${flag}" != "nil" ) then
    set format1 = ".3"
    set format2 = ".3"
else
    set format1 = ""
    set format2 = ""
endif

# Write status information
echo
echo '---- DIRECTORIES AND PROGRAMS ---------------------------'
echo    
echo "PROGRAM CREATE_STARTF : ${LAGRANTO}/startf/create_startf"
echo "PARAMETER file        : ${parfile}"
echo "CRITERION file        : ${crifile}"
echo "RUN directory         : ${tradir}"
echo

# -----------------------------------------------------------------------------------------------------------------
# Set the primary and scecondary data files (necessary for interpolation if intermediate reference date)
# -----------------------------------------------------------------------------------------------------------------

# Find a first data file (if possible corresponding to reference date)
set file=${charp}${yyyy}${mm}${dd}_${hh}
if ( ( -f ${file} ) && ( ${min} == 0 ) ) then
   set timeshift=0
   set date0=${yyyy}${mm}${dd}_${hh}
   set date1=${yyyy}${mm}${dd}_${hh}
   set pfile0=${file1}
   set pfile1=${file}
   goto label3
else
  set file=`ls ${charp}[0-9_]*[0-9] | head -1 | sed -e 's/@//'`
endif

# Determine time increment (in hours) between data files
set date1=`echo $file | cut -c 2-12`
set n=`ls ${charp}[0-9_]*[0-9] | grep -n $date1 | awk -F: '{print $1}'`
@ n ++
set date2=`ls ${charp}[0-9_]*[0-9] | head -$n | tail -1 | cut -c 2-12`
set timeinc=`${LAGRANTO}/goodies/gettidiff $date2 $date1`
endif
if ( $timeinc == 0 ) then
    echo " ERROR: cannot set the time increment between input files ***"
    exit 1
endif

# Search the first file to use
set timeshift=
foreach i ( ${charp}????????_?? )

  set date0 = `echo $i | cut -c 2-12`
  set td1  = `${LAGRANTO}/goodies/gettidiff ${yyyy}${mm}${dd}_${hh} ${date0}`

  if ( ( $td1 >= 0 ) && ( $td1 < $timeinc ) )  then
      set timeshift=$td1
      set pfile0=${charp}${date0}
      goto label2
  endif

end

# Check if no P files are available for the specified time period
if ( $timeshift == "" ) then
  echo " ERROR: no data files available for the specified reference date"
  exit 1
endif

# Everything is fine so far: proceed
label2:

# Set the next date and check its availability
if ( ( ${timeshift} != 0 ) || ( ${min} > 0 ) ) then
    set date1=`${LAGRANTO}/goodies/newtime $date0 $timeinc`
    if ( ! -f ${charp}${date1} ) then
       echo " ERROR: file with primary data is missing for $date1"
       exit 1
    else
       set pfile1=${charp}${date1}
    endif
else
    set date1=${date0}
    set pfile1=${pfile0}
endif

# Set the final timeshift 
if ( ${min} != 00 ) then
    set timeshift=${timeshift}.${min}
endif

# Everything is fine!
label3:

# Set secondary files and check their availability
set sfile0=${chars}${date0}
set sfile1=${chars}${date1}

# Write status information
echo '---- DATA FILES -----------------------------------------'
echo    
echo "Primary files                         : ${pfile0}"
echo "                                      : ${pfile1}"
echo "Secondary files                       : ${sfile0}"
echo "                                      : ${sfile1}"
echo "Timeshift to first data file (hh.mm)  : ${timeshift}"
echo "Time increment of data files          : ${timeinc}"
echo

# --------------------------------------------------------------------------------------------------------------
# Create the start positions (without selection)
# -----------------------------------------------------------------------------------------------------------------

# Set times relative to the reference date
if ( "${changet}" == "true" ) then
  echo '---- CHANGE TIMES ON DATA FILES  ------------------------'
  echo   
  ${LAGRANTO}/goodies/changet.sh ${refdate} ${pfile0}
  ${LAGRANTO}/goodies/changet.sh ${refdate} ${pfile1}
  if ( -f  ${sfile0} ) then
    ${LAGRANTO}/goodies/changet.sh ${refdate} ${sfile0}
  endif
  if ( -f  ${sfile1} ) then
    ${LAGRANTO}/goodies/changet.sh ${refdate} ${sfile1}
  endif
endif

# Write parameters to parameter file and create the starting positions
\rm -f ${parfile}
echo \"${pfile0}\" \"${pfile1}\"    >! ${parfile}
echo \"${sfile0}\" \"${sfile1}\"    >> ${parfile}
echo \"${ofile}${format1}\"         >> ${parfile}
echo \"${regionf}\"                 >> ${parfile}
echo ${yyyy}                        >> ${parfile}
echo ${mm}                          >> ${parfile}
echo ${dd}                          >> ${parfile}
echo ${hh}                          >> ${parfile}
echo ${min}                         >> ${parfile}
echo 00                             >> ${parfile}
echo ${timeshift}                   >> ${parfile}
echo ${timeinc}                     >> ${parfile}

# Analyse the specifier and append to parameter file
${LAGRANTO}/startf/create_startf.perl "${specifier}" >> ${parfile}

if ( "${status}" != "0" ) then
  echo "ERROR:  Preprocessor <create_startf> failed"
  exit 1
endif

# Write selection criterion to file
\rm -f ${crifile}
tail -1 ${parfile} >! ${crifile} 

# Write flag for no time check
echo \"${timecheck}\"               >> ${parfile}

# Write title
echo 
echo '       *** END OF PREPROCESSOR CREATE_STARTF ***'        
echo '========================================================='
echo

# Create the startf
    cd ${tradir}
    pwd
${LAGRANTO}/startf/create_startf

if ( "${status}" != "0" ) then
  echo "ERROR:  Program <create_startf> failed"
  exit 1
endif

# --------------------------------------------------------------------------------------------------------------
# Apply selection (first tracing then selection)
# --------------------------------------------------------------------------------------------------------------

# Stop if no tracing and selection is necessary
if ( "${flag}" == "nil" ) goto finish

# Tracing of extra variables
if ( -f  ${ofile}${format2} ) then
  \rm -f ${ofile}${format2}
endif
if ( "${timecheck}" == "no" ) then
  ${LAGRANTO}/trace/trace.sh ${ofile}${format1} ${ofile}${format2} -v ${tvfile} -notimecheck
else
  ${LAGRANTO}/trace/trace.sh ${ofile}${format1} ${ofile}${format2} -v ${tvfile}
endif
\rm -f ${ofile}${format1}

# Selection
if ( -f  ${ofile} ) then
  \rm -f ${ofile}
endif
${LAGRANTO}/select/select.sh ${ofile}${format2} ${ofile} `cat ${crifile}` 
\rm -f ${ofile}${format2}

# --------------------------------------------------------------------------------------------------------------
# Final tasks (make clean)
# --------------------------------------------------------------------------------------------------------------

finish:

echo $noclean

if ( "${noclean}" == "false" ) then
  \rm -f ${crifile}
  \rm -f ${parfile}
 endif

exit 0 
  


