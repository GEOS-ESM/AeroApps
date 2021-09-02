#!/bin/bash

# Sample commands:
# vlidort_check2.bash s gfortran
# vlidort_check2.bash v ifort
# vlidort_check2.bash sc gfortran

echo
if   [ "$1" = "s" ] ; then
  echo "doing scalar checks"
elif [ "$1" = "v" ] ; then
  echo "doing vector checks"
elif [ "$1" = "sc" ] ; then
  echo "doing spherical correction checks"
fi

# Run VLIDORT checks

cd vlidort_$1_test

file_count=0
echo
echo 'checking results of tests...'
echo
for filename in results_*.* ; do
  ref_filename=saved_results/$2/$filename
  if [ -e $ref_filename ]; then
    ../vlidort_diff $ref_filename $filename
    echo $((++file_count)) files processed
  fi
done
#for filename in results_brdf_*check*.* ; do
#  if [ "$filename" != 'results_brdf_*check*.*' ] ; then
#    ref_filename=saved_results/$2/$filename
#    if [ -e $ref_filename ]; then
#      ../vlidort_diff $ref_filename $filename
#      echo $((++file_count)) files processed
#    fi
#  fi
#done

echo
echo 'done'
echo

exit
