#!/bin/csh -f

if ( $#argv < 1 ) then
    echo "pr_instance.csh  diag_fn"
    exit 1
else

endif

if ( $#argv > 1 ) then
   set fsize = $2
else
   set fsize = 20
endif

set fdir = ${fsize}ha

set diag_fn = $1  # includes dir
set in_dir = `dirname $diag_fn`/../../..
set phis_fn = $in_dir/phis.hdf

set out_fn = `echo $diag_fn | \
              sed -e "s@ceres_b55/diag@pr_v1/$fdir@" \
                  -e 's/c403_cer_01.diag.eta/pr_v1.hghts.sfc/' `

set out_dir = `dirname $out_fn` 
#set basen = `basename $out_fn`
#set basen = $basen:r
#set timestamp = $basen:e
#set nhms = `echo $timestamp | awk '{print substr($1,10,4)"00"}'`

   mkdir -p $out_dir

   if ( -e $out_fn ) then
    echo "pr_instance: output file $out_fn exists, skipping it..."
   else
    PlumeRise_sa.x -fsize $fsize -o $out_fn $phis_fn $diag_fn
   endif

