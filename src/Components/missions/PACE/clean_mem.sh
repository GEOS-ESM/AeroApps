#!/bin/bash

segments=$( ipcs -m | grep '^0x' | awk '{print $1}' )

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