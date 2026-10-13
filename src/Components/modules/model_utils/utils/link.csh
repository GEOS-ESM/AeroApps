#!/bin/tcsh

set i = 0
foreach file (`\ls -1 f5430*png`)
@ i = $i + 1

set istr = $i
if($i<1000) set istr = "0"$istr
if($i<100)  set istr = "0"$istr
if($i<10)   set istr = "0"$istr

echo f5430_fp.inst1_2d_hwl_Nx.$istr.png
ln -s $file f5430_fp.inst1_2d_hwl_Nx.$istr.png

end
