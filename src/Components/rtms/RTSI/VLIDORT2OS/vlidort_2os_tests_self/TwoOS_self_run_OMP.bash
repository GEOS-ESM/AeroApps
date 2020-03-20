#!/bin/bash

# This script runs the TwoOS package tests which have been set up to
# run using OpenMP (see http://openmp.org for further information about OpenMP)

# Sample command:
# TwoOS_self_run_OMP.bash gfortran


# Do some set up

ulimit -s unlimited         # Maximum stack size of the main thread
ulimit -c unlimited         # Maximum size of core file created if a problem occurs

#export OMP_STACKSIZE=50M    # Maximum stack size of OpenMP-spawned threads
export OMP_STACKSIZE=100M
#export OMP_STACKSIZE=150M
#export OMP_STACKSIZE=200M
#export OMP_STACKSIZE=512M


# Run tests using OpenMP

echo
echo "making OpenMP tests ..."
echo
make openmp FC=$1 OPENMP=t $2
echo

fn_prefix='TwoOS_'
fn_suffix='.exe'

for filename in TwoOS_*_OMPtester.exe ; do
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

#rm makefile

echo
echo 'done'
echo

exit 0
