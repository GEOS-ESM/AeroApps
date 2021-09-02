#!/bin/bash

# This script runs the two VLIDORT package solar tests which have been set up to
# run using OpenMP (see http://openmp.org for further information about OpenMP)

# Sample command:
# vlidort_run_OMP.bash gfortran

# Note: gfortran is the only OpenMP 3.0 compatible compiler set up in the makefile
#       at this time (9.28.14)

# Note: (10/16/14) IFORT Also now tested with OpenMP in the makefile

# Do some set up

ulimit -s unlimited         # Maximum stack size of the main thread
ulimit -c unlimited         # Maximum size of core file created if a problem occurs

#export OMP_STACKSIZE=50M    # Maximum stack size of OpenMP-spawned threads
#export OMP_STACKSIZE=100M
#export OMP_STACKSIZE=150M
#export OMP_STACKSIZE=200M
export OMP_STACKSIZE=512M


# Run VLIDORT tests using OpenMP
if [ -e makefile ] ; then
  echo
  #check to see if active makefile is current
  if [ $(diff -q vlidort_v_test/makefile makefile | wc -l) != "0" ] ; then
    echo 'makefile has changed - copying new version up to script directory ...'
    cp vlidort_v_test/makefile .
  else
    echo 'makefile in script directory is up to date ...'
  fi
else
  echo 'copying current test makefile up to script directory ...'
  cp vlidort_v_test/makefile .
fi

# Replace active "vlidort_pars.f90" file with the standard one for these tests if needed
if [ $(diff -q vlidort_def/vlidort_pars.f90_save vlidort_def/vlidort_pars.f90 | wc -l) != "0" ] ; then
  make clean
  cp vlidort_def/vlidort_pars.f90_save vlidort_def/vlidort_pars.f90
fi

echo
echo "making OpenMP tests ..."
echo
make openmp FC=$1 OPENMP=t $2
echo

fn_prefix='v2p8p3_'
fn_suffix='.exe'

for filename in v2p8p3_*_OMP_tester.exe ; do
  #chop off filename frontend
  temp_fn=${filename#$fn_prefix}
  #chop off filename backend
  testname=${temp_fn%$fn_suffix}

  echo
  echo "running $testname ..."
  echo
  ./$filename
  echo
done

rm v2p8p3_*_OMP_tester.exe
#rm makefile

echo
echo 'done'
echo

exit 0
