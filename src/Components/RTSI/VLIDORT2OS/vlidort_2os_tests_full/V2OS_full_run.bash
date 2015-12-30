#!/bin/bash

# Sample commands:
# V2OS_full_run.bash gfortran
# V2OS_full_run.bash ifort


# Run V2OS FULL tests

echo
echo "making main tests ..."
echo
make main FC=$1 $2
echo

fn_prefix="V2OS_"
fn_suffix='.exe'
for filename in V2OS_*_tester.exe ; do
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
