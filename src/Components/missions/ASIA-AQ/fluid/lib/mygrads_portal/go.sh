#!/bin/sh

SRCPATH=/discover/nobackup/aconaty/DEV_TEST_FLUID_20170711/dev/lib

while read file; do

  local_file=`basename $file`

  if [ ! -f $local_file ]; then
    echo "$local_file is missing"
  else
    diff -q $file $local_file
  fi

done

exit 0
