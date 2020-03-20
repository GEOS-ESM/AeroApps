#!/bin/csh -f
#
#PBS -l walltime=8:00:00
#PBS -l select=1:ncpus=8:proc=neha
#PBS -N VLIDORT
#PBS -q general
#PBS -W group_list=s0980

#               ---------------
#               Set Environment
#               ---------------
#

umask 022
limit stacksize unlimited

set ARCH = `uname -s`
set MACH = `uname -m`

setenv ESMADIR __ESMADIR__ 
setenv BASEDIR __BASEDIR__

set path = ( $ESMADIR/$ARCH/bin $BASEDIR/$ARCH/bin $path )

if ( $?PYTHONPATH ) then
     setenv PYTHONPATH $ESMADIR/$ARCH/lib/Python:$PYTHONPATH
else
     setenv PYTHONPATH $ESMADIR/$ARCH/lib/Python
endif

#               --------------------
#                Run the Simulations
#               --------------------
#

foreach nymd ( 20080701 )
    foreach nhms ( 120000 )

       omi_l3a.py $nymd nhms

    end
end


