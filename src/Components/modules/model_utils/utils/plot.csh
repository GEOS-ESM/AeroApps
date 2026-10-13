#!/bin/csh

# Set base location

# setup environment
  setenv AEROAPPS /home/pcolarco/geos_aerosols/pcolarco/AeroApps
  source $AEROAPPS/env@/g5_modules
  setenv PYTHONPATH ./:$AEROAPPS/install/lib/Python/scat:$AEROAPPS/install/lib/Python/:$AEROAPPS/install/lib/Python/pyobs:$AEROAPPS/install/lib/Python/pyods
  set path = ( . $AEROAPPS/install/bin $path )
# ? Why do I need this next line?
  setenv LD_LIBRARY_PATH $AEROAPPS/install/lib/:$LD_LIBRARY_PATH

# run sampling scripts
  set files = `\ls -1 /home/pcolarco/fp/chem/Y2026/M0[6,7,8]/f5430_fp.inst1_2d_hwl_Nx.2026*.nc4`

  prund.pl -H `hostname` -d `echo $files` &
  mpirun -np 48 prund.pl -H `hostname` ./make_nifty_plots.py %s

#  foreach track (`\ls -1 ../data/*R0.ict`)
#    echo $track
#    ./plot_ext_curtain.py $track >> /dev/null
#   end

