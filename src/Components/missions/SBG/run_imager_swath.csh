#!/bin/csh -f


umask 022

limit stacksize unlimited

#######################################################################
#           Architecture Specific Environment Variables
#######################################################################
source $HOME/.cshrc
source ./setup_env


##################################################################
######
######         Do Sampling
######
##################################################################
python3 -u imager_swath.py -v 2024-01-16T18:00 2024-01-16T19:00 -d 200 --dt_units milliseconds -D 1 --out_year 2006 --no_ss imager_files.pcf ssd650.pcf sbg.pcf 
