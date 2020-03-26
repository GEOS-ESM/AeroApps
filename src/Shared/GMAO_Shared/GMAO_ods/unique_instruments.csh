#!/bin/csh 
#
# 20Mar2009 - Todling - Initial code (there has to be a smart why to do this in perl
#

if ( !($?FVWORK)  ) then
  echo "Error (unique_diag.csh): should define env var FVWORK"
  exit 1
endif
unalias cd
cd $FVWORK
set lst = `/bin/ls *diag_*_anl.*bin`
set nlst = ""
foreach fn ( $lst )
  set chop = `echo $fn | cut -d. -f2 | cut -c6-`
  set len = `echo $chop | wc -c` 
  @ lene = $len - 5
  set chop = `echo $chop | cut -c1-$lene`
  set nlst = ( $nlst $chop )
end
echo $nlst
