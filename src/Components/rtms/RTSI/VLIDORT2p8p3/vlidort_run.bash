#!/bin/bash

# Sample commands:
# vlidort_run.bash s gfortran
# vlidort_run.bash v ifort


echo
if [ "$1" = "s" ] ; then
  echo "doing scalar tests"
elif [ "$1" = "v" ] ; then
  echo "doing vector tests"
fi


# Run VLIDORT tests
if [ -e makefile ] ; then
  echo
  #check to see if active makefile is current
  if [ $(diff -q vlidort_$1_test/makefile makefile | wc -l) != "0" ] ; then
    echo 'makefile has changed - copying new version up to script directory ...'
    cp vlidort_$1_test/makefile .
  else
    echo 'makefile in script directory is up to date ...'
  fi
else
  echo 'copying current test makefile up to script directory ...'
  cp vlidort_$1_test/makefile .
fi

# Replace active "vlidort_pars.f90" file with the standard one for these tests if needed
if [ $(diff -q vlidort_def/vlidort_pars.f90_save vlidort_def/vlidort_pars.f90 | wc -l) != "0" ] ; then
  make clean
  cp vlidort_def/vlidort_pars.f90_save vlidort_def/vlidort_pars.f90
fi

#do main tests:
echo
echo "making main tests ..."
echo
make main FC=$2 $3
echo

fn_prefix="$12p8p3_"
fn_suffix='.exe'
for filename in $12p8p3_*_tester.exe ; do
  #chop off filename frontend
  temp_fn=${filename#$fn_prefix}
  #chop off filename backend
  testname=${temp_fn%$fn_suffix}

  echo "running $testname ..."
  echo
  ./$filename
  echo
done

#do additional special tests (common to both scalar & vector):

echo
make clean

#Planetary test

#replace standard "vlidort_pars.f90" file with special one for this test
cp vlidort_def/vlidort_pars.f90_Planetary_test vlidort_def/vlidort_pars.f90

echo
echo 'making Planetary_tester ...'
echo
make $12p8p3_Planetary_tester.exe FC=$2 $3
echo
echo 'running Planetary_tester ...'
echo
./$12p8p3_Planetary_tester.exe

echo
make clean

#LWCoupling test

#replace "vlidort_pars.f90" file with special one for this test
cp vlidort_def/vlidort_pars.f90_LWCoupling vlidort_def/vlidort_pars.f90

echo
echo 'making LWCoupling_tester ...'
echo
make $12p8p3_LWCoupling_tester.exe FC=$2 $3
echo
echo 'running LWCoupling_tester ...'
echo
./$12p8p3_LWCoupling_tester.exe

echo
make clean

if [ "$1" = "s" ] ; then
  #do additional scalar-only test(s):

  #(1) F-matrix / Z-matrix test

  #return standard "vlidort_pars.f90" file to its place
  cp vlidort_def/vlidort_pars.f90_save vlidort_def/vlidort_pars.f90

  echo
  echo 'making vfzmat_tester ...'
  echo
  make s2p8p3_vfzmat_tester.exe FC=$2 $3
  echo
  echo 'running vfzmat_tester ...'
  echo
  ./s2p8p3_vfzmat_tester.exe

  echo
  make clean

  #(2) Scalar surface-leaving (VSLEAVE) self test

  #replace standard "vlidort_pars.f90" file with special one for this test
  cp vlidort_def/vlidort_pars.f90_vsleave_test vlidort_def/vlidort_pars.f90

  echo
  echo 'making vsleave_self_tester ...'
  echo
  make s2p8p3_vsleave_self_tester.exe FC=$2 $3
  echo
  echo 'running vsleave_self_tester ...'
  echo
  ./s2p8p3_vsleave_self_tester.exe
  
  echo
  make clean

  #return standard "vlidort_pars.f90" file to its place
  cp vlidort_def/vlidort_pars.f90_save vlidort_def/vlidort_pars.f90
fi

if [ "$1" = "v" ] ; then
  #do additional vector-only test(s):

  #(1) Vector surface-leaving (VSLEAVE) test

  #replace standard "vlidort_pars.f90" file with special one for this test
  cp vlidort_def/vlidort_pars.f90_vsleave_test vlidort_def/vlidort_pars.f90

  echo
  echo 'making vsleaveplus_tester ...'
  echo
  make v2p8p3_vsleaveplus_tester.exe FC=$2 $3
  echo
  echo 'running vsleaveplus_tester ...'
  echo
  ./v2p8p3_vsleaveplus_tester.exe

  echo
  make clean

  #(2) Siewert (2000) validation

  #replace "vlidort_pars.f90" file with special one for this test
  cp vlidort_def/vlidort_pars.f90_Siewert2000 vlidort_def/vlidort_pars.f90

  echo
  echo 'making Siewert2000_tester ...'
  echo
  make v2p8p3_Siewert2000_tester.exe FC=$2 $3
  echo
  echo 'running Siewert2000_tester ...'
  echo
  ./v2p8p3_Siewert2000_tester.exe
  
  echo
  make clean

  #return standard "vlidort_pars.f90" file to its place
  cp vlidort_def/vlidort_pars.f90_save vlidort_def/vlidort_pars.f90
fi

#rm makefile

echo
echo 'done'
echo

exit 0
