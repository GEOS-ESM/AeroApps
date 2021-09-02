#!/bin/bash

# New spherical correction tests for VLIDORT 2.8.3. 31 March 2021

# Sample commands:
# vlidort_run_SpherCorr.bash gfortran
# vlidort_run_SpherCorr.bash ifort


# Select desired VLIDORT test(s) to run
# Note to user: just set to 1 to activate desired test(s),
#               then execute the script

           #Test(s)
           #-------

# (1) Individual tests

# main:

test[1]=1  #SpherCorr
test[2]=0  #LCS_SpherCorr
test[3]=0  #LPS_SpherCorr


#### Run the selected test(s) ####

# Define test names:

# (1) Individual tests

testname[1]='SpherCorr_tester_V3'
testname[2]='LCS_SpherCorr_tester_V2'
testname[3]='LPS_SpherCorr_tester_V2'


echo
echo "doing spherical correction test(s)"

# Define executable file prefix and suffix
fn_prefix='V2p8p3_'
fn_suffix='.exe'

# Run desired VLIDORT test(s)
if [ -e makefile ] ; then
  echo
  #check to see if active makefile is current
  if [ $(diff -q vlidort_sc_test/makefile makefile | wc -l) != "0" ] ; then
    echo 'makefile has changed - copying new version up to script directory ...'
    cp vlidort_sc_test/makefile .
  else
    echo 'makefile in script directory is up to date ...'
  fi
else
  echo 'copying current test makefile up to script directory ...'
  cp vlidort_sc_test/makefile .
fi

for ((i=1 ; i<=3 ; i++)) ; do
  if [ "${test[i]}" = "1" ] ; then
    #replace active "vlidort_pars.f90" file with special one for these tests if needed
    if [ $(diff -q vlidort_def/vlidort_pars.f90_Sphericity_Correction_Tests vlidort_def/vlidort_pars.f90  | wc -l) != "0" ] ; then
      make clean
      cp vlidort_def/vlidort_pars.f90_Sphericity_Correction_Tests vlidort_def/vlidort_pars.f90
    fi

    #do selected test
    exec="${fn_prefix}${testname[i]}${fn_suffix}"
    echo $exec
    echo
    echo "making ${testname[i]} ..."
    echo
    make $exec FC=$1 $2
    echo
    echo "running ${testname[i]} ..."
    echo
    ./$exec
    echo
  fi  
done

echo
echo 'done'
echo

exit 0
