#!/bin/csh

# Set base location
#  cd /home/pcolarco/geos_ae

# setup environment
  setenv AEROAPPS /home/pcolarco/geos_aerosols/pcolarco/AeroApps
  source $AEROAPPS/env@/g5_modules
  setenv PYTHONPATH ./:$AEROAPPS/install/lib/Python/scat:$AEROAPPS/install/lib/Python/:$AEROAPPS/install/lib/Python/pyobs:$AEROAPPS/install/lib/Python/pyods
  set path = ( . $AEROAPPS/install/bin $path )
# ? Why do I need this next line?
  setenv LD_LIBRARY_PATH $AEROAPPS/install/lib/:$LD_LIBRARY_PATH

  ln -s /discover/nobackup/projects/gmao/share/dao_ops/fvInput_nc3/ ExtData
  
# run sampling scripts
  set tracks = `\ls -1 /discover/nobackup/projects/gmao/iesa/aerosol/campaigns/CPEXCV/*.ict`
  echo $tracks

  prund.pl -H `hostname` -d `echo $tracks` &
  mpirun -np 48 prund.pl -H `hostname` ./plot_trj_ext_curtain.py %s MERRA2

#  foreach track (`\ls -1 ../data/*R0.ict`)
#    echo $track
#    ./plot_ext_curtain.py $track >> /dev/null
#   end

