#!/bin/csh

setenv EXPID  f513a_rt
setenv FVHOME /discover/nobackup/$user/$EXPID
setenv FVROOT /discover/nobackup/projects/gmao/advda/$user/4OPS/g5151/GEOSadas/Linux
setenv FVWORK /discover/nobackup/$user/AODOBS
setenv GID g0613
setenv group_list "SBATCH -A $GID"

setenv NCPUS_AOD   8
setenv MPIRUN_AOD  "mpirun "

set path = ( . $FVROOT/bin $path )

if (! -d $FVWORK ) mkdir -p $FVWORK
cd $FVWORK

/bin/cp $FVHOME/run/obsys-gaas.rc .

get_aero_obs.csh  20151128 180000 1
get_aero_obs.csh  20151128 210000 1

run_gaas_ana.csh  $EXPID 20151128 180000 1 

