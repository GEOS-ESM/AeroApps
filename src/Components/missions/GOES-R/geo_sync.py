#!/usr/bin/env python3

"""
	Utility to make one scantime file from decomposed
	domain files.
"""
import os

from optparse        import OptionParser
import numpy as np
from netCDF4 import Dataset
import math as math
#---
def getCoords(options, inst='GOES-R'):
  """
  Get coordinates from instrument geo-location file.
  Assumes the file has albed been tighten in the E-W domain
  """
  nc = Dataset(options.geoFile)

  lon = nc.variables['clon'][:,:]
  lat = nc.variables['clat'][:,:]
  missing = nc.variables['clon'].missing_value
  return (nc,lon,lat,missing)
#------------------------------------ M A I N ------------------------------------

if __name__ == "__main__":
	outFile   = 'goes-r.lg1.scantime.nc4'
	geoFile   = '/nobackup/GOES-R/LevelG/invariant/goes-r.lg1.invariant.nc4'
	npet      = 1200

	#   Parse command line options
	#   --------------------------
	parser = OptionParser(usage="Usage: %prog [OPTIONS] ", version='1.0.0' )

	parser.add_option("-o", "--output", dest="outFile", default=outFile,
			              help="Output NetCDF file (default=%s)"\
	                        %outFile )
	parser.add_option("-n", "--np", dest="npet", default=npet,
			              help="number of processors (default=%i)"\
	                        %npet )
	parser.add_option("-g", "--geolocation", dest="geoFile", default=geoFile,
              help='GOES-R Level G1 geo-location file (default=%s)'%geoFile)

	(options, args) = parser.parse_args()
	npet      = options.npet
	outFile   = options.outFile

	# Read in Data
	# ------------
	ncGeo, clon, clat, MISSING = getCoords(options)
	ym, xm = clon.shape

	nclr    = np.zeros(npet)
	if (npet >= xm):
		nclr[0:xm] = 1
	else: 
	  nclr[:] = xm/npet
	  left = xm % npet
	  if (npet >= left):
	    nclr[0:left] = nclr[0:left] + 1
	  else: 
	    while (left > 0): 
	      if (npet >= left):
	        nclr[0:left] = nclr[0:left] + 1
	        left = 0
	      else:
	        nclr = nclr + 1
	        left = left - npet

	counti = nclr
	starti = np.zeros(npet)
	starti[0] = 1
	for i in np.arange(npet-1):
		starti[i+1] = starti[i] + counti[i] 
	endi   = starti + counti -1

	starti = starti.astype(int)
	endi   = endi.astype(int)


	# -- Create outfile
	ncout = Dataset(outFile,'w',format='NETCDF4')

	# -- Global attributes
	ncout.institution = "NASA/Goddard Space Flight Center"
	ncout.source      = "Global Model and Assimilation Office"
	ncout.history     = "Created from geo_sync.F90 output"
	ncout.comment     = "This file contains scanTimes for the GOES-R domain synchronized with TEMPO"
	ncout.contact     = "P.Castellanos<patricia.castellanos@nasa.gov>"

	# -- Define dimensions and variables
	ns    = ncout.createDimension('ns', ym)
	ew    = ncout.createDimension('ew', xm)
	scanTime = ncout.createVariable('scanTime','f8',('ns','ew',))
	scanTime.units = 'seconds since 0001-01-01 00:00:00.0'
	scanTime.missing_value = 1e15

	# -- Read in data and write to variable
	for t,start in enumerate(starti):
		end   = endi[t]
		ncin  = Dataset('sync/goes-r.lg1.scantime.'+str(start).zfill(4)+'-'+str(end).zfill(4)+'.nc4')
		scanTime[:,start-1:end]  = ncin.variables['scanTime'][:,:]
		ncin.close()

	ncout.close()

