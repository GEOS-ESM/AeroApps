#!/bin/sh

fname=$1
ncks=$BASEDIR/Linux/bin/ncks

$ncks -O -d latitude,15.0,70.0 -d longitude,190.0,310.0 $fname $fname.tmp
#/bin/mv -f $fname.tmp $fname

exit 0
