#!/bin/bash

# Sample commands:
# vlidort_run_extras.bash v ifort nstokes3
# vlidort_run_extras.bash v gfortran obsgeo
# vlidort_run_extras.bash v ifort doublet

echo
if [ "$1" = "v" ] ; then
  echo "doing vector $3 extra tests"
else
  echo 'the 1st argument of vlidort_run_extras.bash must be "v"'
  exit 1
fi

# Check to see if active makefile exists and is current
echo
if [ -e makefile ] ; then
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

# Prepare the tests
echo
rm -f *.exe

# Define test-specific config files for the desired tests
config_list="V2p8p3_VLIDORT_ReadInput VBRDF_ReadInput"
if [ "$3" = "nstokes3" ] ; then
  #config_list="V2p8p3_VLIDORT_ReadInput VBRDF_ReadInput"
  fn_suffix='ns3'
elif [ "$3" = "obsgeo" ] || [ "$3" = "doublet" ] ; then
  #config_list="V2p8p3_VLIDORT_ReadInput VBRDF_ReadInput V2p8p3_vsleaveplus_tester VSLEAVE_ReadInput"
  fn_suffix=$3
else
  echo 'the 3rd argument of vlidort_run_extras.bash must be "nstokes3", "obsgeo", or "doublet"'
  exit 1
fi

# Copy those test-specific config files up to the "vlidort_v_test" subdir
for fn in $config_list ; do
  #echo "copying $fn config file up to vlidort_v_test"
  cp vlidort_v_test/config/$fn_suffix/$fn.cfg_$fn_suffix vlidort_v_test/$fn.cfg
done

# Make some test executables
echo "making $3 extra tests ..."
echo
make $3 FC=$2 $4

# Run VLIDORT extra tests
echo
fn_prefix="v2p8p3_"
fn_suffix='.exe'
for filename in v2p8p3_*_tester.exe ; do
  #chop off filename frontend
  temp_fn=${filename#$fn_prefix}
  #chop off filename backend
  testname=${temp_fn%$fn_suffix}

  echo "running $testname ..."
  echo
  ./$filename
  echo
done

#if [ "$3" = "obsgeo" ] || [ "$3" = "doublet" ] ; then
##if [ "$3" = "dummy" ] ; then
#  #do additional vector test(s) requiring special "vlidort.pars.f90" settings:
#
#  echo
#  make clean
#
#  #Vector surface-leaving (VSLEAVE) test
#
#  #replace standard "vlidort_pars.f90" file with special one for this test
#  cp vlidort_def/vlidort_pars.f90_vsleave_test vlidort_def/vlidort_pars.f90
#
#  echo
#  echo 'making vsleaveplus_tester ...'
#  echo
#  make v2p8p3_vsleaveplus_tester.exe FC=$2 $4
#  echo
#  echo 'running vsleaveplus_tester ...'
#  echo
#  ./v2p8p3_vsleaveplus_tester.exe
#
#  #return standard "vlidort_pars.f90" file to its place
#  cp vlidort_def/vlidort_pars.f90_save vlidort_def/vlidort_pars.f90
#
#  echo
#  make clean
#fi

# Return default ("base") config files back to the "vlidort_v_test" subdir
for fn in $config_list ; do
  #echo "returning $fn default ('base') config file to vlidort_v_test"
  cp vlidort_v_test/config/base/$fn.cfg_base vlidort_v_test/$fn.cfg
done

rm -f *.exe
#rm makefile

echo
echo 'done'
echo

exit 0
