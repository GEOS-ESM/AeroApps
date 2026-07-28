#!/bin/csh

# Set base location
  cd /home/pcolarco/ARCSIX/eval

# setup environment
  setenv AEROAPPS $NOBACKUP/AeroApps
  source $AEROAPPS/env@/g5_modules
  setenv PYTHONPATH ./:$AEROAPPS/install/lib/Python/scat:$AEROAPPS/install/lib/Python/:$AEROAPPS/install/lib/Python/pyobs:$AEROAPPS/install/lib/Python/pyods
  set path = ( . $AEROAPPS/install/bin $path )
# ? Why do I need this next line?
  setenv LD_LIBRARY_PATH $AEROAPPS/install/lib/:$LD_LIBRARY_PATH

# run sampling scripts
  set tracks = `\ls -1 ../data/*R0.ict`
  prund.pl -H `hostname` -d `echo $tracks` &
  mpirun -np 48 prund.pl -H `hostname` ./plot_trj_mass_along.py %s fp

#  foreach track (`\ls -1 ../data/*R0.ict`)
#    echo $track
#    ./plot_mass_curtain.py $track >> /dev/null
#   end

