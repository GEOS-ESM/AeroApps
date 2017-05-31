#!/bin/csh

if ( $#argv < 1 ) then
  echo 
  ${LAGRANTO}/bin/lagrantohelp changet short
  echo  
  exit 0
endif

set stdate=$1
set fort=fort.9

if ( $#argv == 2 ) then
  set string=$2
else
  set string='P[_0-9]*[0-9] S[_0-9]*[0-9]'
endif

foreach i ( $string )

  \rm -f $fort
  touch $fort

# set date = `echo $i | sed -e 's/[A-Za-z_]*//' | cut -c 1-9`
  set date = `echo $i | sed -e 's/[A-Za-z_]*//'`
  echo $i >> $fort

  set date=`echo ${date}00 | cut -c 1-13` 
  set stdate=`echo ${stdate}00 | cut -c 1-13`

  set tim = `${LAGRANTO}/goodies/gettidiff $date $stdate`
  echo ${tim} >> $fort

  ${LAGRANTO}/goodies/changet

end

#\rm -f $fort
