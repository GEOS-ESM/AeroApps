#!/bin/csh -f

setenv MYNAME RmShmKeys_sshmpi

                                    setenv FAILED 0
if( (! $?FVROOT ) & (! $?GEOSBIN) ) setenv FAILED 1

if ( $FAILED ) then
  env
  echo " ${MYNAME}: not all required env vars defined"
  exit 1
endif

if( $?FVROOT ) then
     set pathname = $FVROOT/bin
endif
if( $?GEOSBIN ) then
     set pathname = $GEOSBIN
endif

if( $?PBS_NODEFILE ) then
   sleep 10

   set nodes = `cat $PBS_NODEFILE | uniq`
   foreach node ($nodes)
      echo sshmpi $node $pathname/rmshmkeyhere.sh
           sshmpi $node $pathname/rmshmkeyhere.sh &
   end

   wait
endif
