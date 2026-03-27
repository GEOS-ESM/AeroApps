#!/bin/csh -f
#
# This script runs the gt_generate_parcels utility
# to create a netcdf file of initial parcel positions
# for the gigatraj model.
#
# We will create a set of parcels randomly distributed 
# within a disc of radius $R centered at ($Lat, $Lon),
# between potential temperature of 380.0 and 390.0 K.

# Call is
# ./initialize_profile.csh outstring lon lat z
# where outstring is added to the file name "profile_init.{$outstring}.nc"
# lon and lat are the central longitude and latitude for a disk of radius
# R (set below) and z is the midpoint (in km) of the parcel positions (+/- 1 km)

set outstring = $1

# define the disc
set R   = 100.0
set Lat = $3
set Lon = $2

# define the vertical coordinates and range
set VCoord = "PAlt"
set VUnits = "km"
set z      = $4

#Found this here:
#https://earthsci.stanford.edu/computing/unix/programming/shell/expressions.php
alias MATH 'set \!:1 = `echo "\!:3-$" | bc -l`'
MATH VMin = $z - 1
MATH VMax = $z + 1

# The name of the output text file
set OutputFile="parcel_init.$outstring.nc"
echo $outstring

# Use this format in the output:
#  lOngitude   lAtitude   Vertical
set Format="%o %a %v"

# Generate the parcels
${ROOT}/bin/gt_generate_parcels \
    --number 2500 \
    --random \
    --clat ${Lat} \
    --clon ${Lon} \
    --radius ${R} \
    --vertical ${VCoord} \
    --vunits ${VUnits} \
    --zlow ${VMin} \
    --zhigh ${VMax} \
    --format "${Format}" \
    --netcdf "${OutputFile}"



