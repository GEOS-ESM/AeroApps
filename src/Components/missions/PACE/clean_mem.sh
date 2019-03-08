#!/bin/bash

segments=$( ipcs -m | grep $USER | grep '^0x1b2e' | awk '{print $1}' )

if [[ -z $segments ]]
then
  echo "No shared memory segments to remove."
else
  echo "Removing segments..."
  for seg in $segments
  do
     ipcrm -M $seg
  done
fi
