#!/bin/sh

if [ $# -ge 1 ]; then
  idate=$1
else
  idate=`date "+%Y%m%d"`
fi

echo $idate
