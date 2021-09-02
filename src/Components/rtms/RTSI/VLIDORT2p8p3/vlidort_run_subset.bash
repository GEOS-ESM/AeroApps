#!/bin/bash

# Sample commands:
# vlidort_run_subset.bash s gfortran
# vlidort_run_subset.bash v ifort


# Select desired VLIDORT test(s) to run
# Note to user: just set to 1 to activate desired test(s),
#               then execute the script

           #Test(s)
           #-------

# (1) Individual tests

# main (both scalar & vector):
test[1]=1  #solar
test[2]=0  #solar_lpcs
test[3]=0  #brdf_self
test[4]=0  #thermal
test[5]=0  #thermal_lpcs
test[6]=0  #brdfplus

# special (both scalar & vector):
test[7]=0  #Planetary        ---> NEW Version 2.8.1
test[8]=0  #LWCoupling       ---> NEW Version 2.8.1

# special (scalar-only tests):
test[9]=0  #vfzmat           ---> NEW Version 2.8
test[10]=0 #vsleave_self     ---> NEW Version 2.8

# special (vector-only tests):
test[11]=0 #vsleaveplus
test[12]=0 #Siewert2000

# (2) Set of tests
test[13]=0 #all solar
test[14]=0 #all thermal
test[15]=0 #all main


#### Run the selected test(s) ####

# Define test names:

# (1) Individual tests 

# main (both scalar & vector):
testname[1]='solar'
testname[2]='solar_lpcs'
testname[3]='brdf_self'
testname[4]='thermal'
testname[5]='thermal_lpcs'
testname[6]='brdfplus'

# special (both scalar & vector):
testname[7]='Planetary'
testname[8]='LWCoupling'

# special (scalar-only tests):
testname[9]='vfzmat'
testname[10]='vsleave_self'

# special (vector-only tests):
testname[11]='vsleaveplus'
testname[12]='Siewert2000'

# (2) Sets of tests
testname[13]='solar'
testname[14]='thermal'
testname[15]='main'


echo
if [ "$1" = "s" ] ; then
  echo "doing scalar test(s)"
elif [ "$1" = "v" ] ; then
  echo "doing vector test(s)"
fi

# Define executable file prefix and suffix
fn_prefix="$12p8p3_"
fn_suffix='.exe'

# Run desired VLIDORT test(s)
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

for ((i=1 ; i<=15 ; i++)) ; do
  #define file executable(s) to run
  #echo 'i = '$i
  #echo 'test[i] = '${test[i]}

  if [ $i -le 12 ] ; then
    #do individual test
    #echo 'inside individual test section'

    if [ "${test[i]}" = "1" ] ; then
      if [ $i -eq 7 ] ; then
        #replace active "vlidort_pars.f90" file with special one for
        #  test #7 (Planetary test) if needed
        #echo 'dud line before'
        if [ $(diff -q vlidort_def/vlidort_pars.f90_Planetary_test vlidort_def/vlidort_pars.f90 | wc -l) != "0" ] ; then
          make clean
          cp vlidort_def/vlidort_pars.f90_Planetary_test vlidort_def/vlidort_pars.f90
        fi
      elif [ $i -eq 8 ]  ; then
        #replace active "vlidort_pars.f90" file with special one for
        #  test #8 (LWCoupling test) if needed
        #echo 'dud line before'
        if [ $(diff -q vlidort_def/vlidort_pars.f90_LWCoupling vlidort_def/vlidort_pars.f90 | wc -l) != "0" ] ; then
          make clean
          cp vlidort_def/vlidort_pars.f90_LWCoupling vlidort_def/vlidort_pars.f90
        fi
      elif [ $i -eq 10 ] && [ "$1" = "s" ] ; then
        #replace active "vlidort_pars.f90" file with special one for
        #  scalar test #10 (vsleave self test) if needed
        #echo 'dud line before'
        if [ $(diff -q vlidort_def/vlidort_pars.f90_vsleave_test vlidort_def/vlidort_pars.f90 | wc -l) != "0" ] ; then
          make clean
          cp vlidort_def/vlidort_pars.f90_vsleave_test vlidort_def/vlidort_pars.f90
        fi
      elif [ $i -eq 11 ] && [ "$1" = "v" ] ; then
        #replace active "vlidort_pars.f90" file with special one for
        #  vector test #11 (vsleaveplus test) if needed
        #echo 'dud line before'
        if [ $(diff -q vlidort_def/vlidort_pars.f90_vsleave_test vlidort_def/vlidort_pars.f90 | wc -l) != "0" ] ; then
          make clean
          cp vlidort_def/vlidort_pars.f90_vsleave_test vlidort_def/vlidort_pars.f90
        fi
      elif [ $i -eq 12 ] && [ "$1" = "v" ] ; then
        #replace active "vlidort_pars.f90" file with special one for
        #  vector test #12 (Siewert2000 test) if needed
        #echo 'dud line before'
        if [ $(diff -q vlidort_def/vlidort_pars.f90_Siewert2000 vlidort_def/vlidort_pars.f90 | wc -l) != "0" ] ; then
          make clean
          cp vlidort_def/vlidort_pars.f90_Siewert2000 vlidort_def/vlidort_pars.f90
        fi
      else
        #replace active "vlidort_pars.f90" file with the standard one
        #  "vlidort_pars.f90_save" if needed
        #echo 'dud line before'
        if [ $(diff -q vlidort_def/vlidort_pars.f90_save vlidort_def/vlidort_pars.f90 | wc -l) != "0" ] ; then
          make clean
          cp vlidort_def/vlidort_pars.f90_save vlidort_def/vlidort_pars.f90
        fi
      fi

      #do selected test
      exec="${fn_prefix}${testname[i]}_tester${fn_suffix}"
      #echo $exec
      echo
      echo "making ${testname[i]}_tester ..."
      echo
      make $exec FC=$2 $3
      echo
      #exit 0
      echo "running ${testname[i]}_tester ..."
      echo
      ./$exec
      echo
    fi

  elif [ $i -ge 13 ] && [ $i -le 14 ] ; then
    #do a set of tests (a subset of "main")
    #echo 'inside set of tests section'
    if [ "${test[i]}" = "1" ] ; then
      #replace active "vlidort_pars.f90" file with the standard one if needed
      #echo 'dud line before'
      if [ $(diff -q vlidort_def/vlidort_pars.f90_save vlidort_def/vlidort_pars.f90 | wc -l) != "0" ] ; then
        make clean
        cp vlidort_def/vlidort_pars.f90_save vlidort_def/vlidort_pars.f90
      fi

      #do selected set
      echo
      echo "making ${testname[i]} tests ..."
      echo
      make ${testname[i]} FC=$2 $3
      echo
      filename_template=${fn_prefix}${testname[i]}*tester${fn_suffix}
      for filename in ${filename_template} ; do
        echo "filename = $filename"
        #chop off filename frontend
        temp_fn=${filename#$fn_prefix}
        #chop off filename backend
        fn=${temp_fn%$fn_suffix}

        echo "running $fn ..."
        echo
        ./$filename
        echo
      done
    fi
  elif [ $i -eq 15 ] ; then
    #do all "main" tests
    #echo 'inside set of tests section'
    if [ "${test[i]}" = "1" ] ; then
      #replace active "vlidort_pars.f90" file with the standard one if needed
      #echo 'dud line before'
      if [ $(diff -q vlidort_def/vlidort_pars.f90_save vlidort_def/vlidort_pars.f90 | wc -l) != "0" ] ; then
        make clean
        cp vlidort_def/vlidort_pars.f90_save vlidort_def/vlidort_pars.f90
      fi

      #do selected set
      echo
      echo "making ${testname[i]} tests ..."
      echo
      make ${testname[i]} FC=$2 $3
      echo
      for ((j=1 ; j<=6 ; j++)) ; do
        exec="${fn_prefix}${testname[j]}_tester${fn_suffix}"

        echo "running ${testname[j]}_tester ..."
        echo
        ./$exec
        echo
      done
    fi
  fi
done

#rm makefile

echo
echo 'done'
echo

exit 0
