#!/bin/csh

setenv EXPID  prePP_rt
setenv FVHOME /discover/nobackup/projects/gmao/advda/$user/$EXPID
setenv FVROOT `cat $FVHOME/.FVROOT`
setenv FVWORK /discover/nobackup/$user/AODWORK
setenv GID g0613
setenv MODIS_L2_HDF 1
setenv group_list "SBATCH -A $GID"

setenv NCPUS_AOD   8
setenv MPIRUN_AOD  "mpiexec_mpt "

set path = ( . $FVROOT/bin $path )
source $FVROOT/bin/g5_modules

if (! -d $FVWORK ) mkdir -p $FVWORK
cd $FVWORK

/bin/cp $FVHOME/run/obsys-gaas.rc .

get_aero_obs.csh  20160604 180000 1
get_aero_obs.csh  20160604 210000 1

run_gaas_ana.csh  $EXPID 20160604 180000 1  $FVWORK

