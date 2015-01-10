#!/bin/tcsh -f

set pid = $1 

sleep 5
kill -9 $pid
