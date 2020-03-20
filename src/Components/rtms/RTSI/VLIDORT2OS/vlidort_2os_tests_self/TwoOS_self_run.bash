#!/bin/bash

# Sample commands:
# TwoOS_self_run.bash gfortran
# TwoOS_self_run.bash ifort


# Run TwoOS_SELF tests

echo
echo "making main tests ..."
echo
make main FC=$1 $2

fn_prefix="TwoOS_"
fn_suffix='.exe'
for filename in TwoOS_*_tester.exe ; do
  #chop off filename frontend
  temp_fn=${filename#$fn_prefix}
  #chop off filename backend
  testname=${temp_fn%$fn_suffix}
  echo
  echo "running $testname ..."
  ./$filename
done

#rm makefile

echo
echo 'done'
echo

exit 0
