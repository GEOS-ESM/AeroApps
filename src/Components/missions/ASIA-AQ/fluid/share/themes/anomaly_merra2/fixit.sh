#!/bin/sh

while read fname; do
    ncatted -a Title,global,m,c,'No Title\n' $fname
done

exit 0
