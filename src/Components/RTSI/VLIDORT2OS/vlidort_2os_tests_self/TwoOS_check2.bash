#!/bin/bash

# Sample command:
# TwoOS_check2.bash gfortran
# TwoOS_check2.bash  ifort


# Run VLIDORT 2OS checks

file_count=0
echo
echo 'checking results of tests...'
echo
for filename in TwoOS_*.all* ; do
  ref_filename=saved_results/$2/$filename
  if [ -e $ref_filename ]; then
    TwoOS_diff $ref_filename $filename
    echo $((++file_count)) files processed
  fi
done

echo
echo 'done'
echo

exit 0
