#!/bin/python
# -*- coding: utf-8 -*-
""" Runscript for geo_sampler"""

import os
from datetime import datetime, timedelta 
from dateutil.parser import parse
from distutils.dir_util import mkpath
import numpy as np


#########################################################

if __name__ == "__main__":
    
    startdate = '2005-01'
    enddate   = '2005-12'
    dataname  = 'SurfLER'
    instname  = 'omi'

    calculon  = '/nobackup/TEMPO'
    nccs      = '/discover/nobackup'

    dt = 1
    startdate = parse(startdate+'-01T00:00:00')
    enddate   = parse(enddate+'-01T00:00:00')
    
    if (os.path.exists(calculon)):
        indir     = calculon + '/' + dataname
        outdir    = indir
    else:
        indir     = nccs + '/projects/gmao/osse2/pub/c1440_NR/OBS/TEMPO/DATA/' + dataname
        outdir    = indir

    while (startdate <= enddate):
    	date = str(startdate.month).zfill(2)
    	command = 'geo_sampler.py -C -N -v -r geo_sampler_omi_ler.rc -o ' + outdir + \
    					'/tempo-'+ instname + '.' + dataname + '.'+ date + '.nc4 ' +\
    					str(startdate.date()) + 'T' + str(startdate.time())
    	print command
    	os.system(command)    					
    	
    	if (startdate.month < 12):
    		startdate = str(startdate.year) + '-' + str(startdate.month+dt).zfill(2) + '-01T00:00:00'
    	else:
    		startdate = str(startdate.year + 1) + '-01-01T00:00:00'

    	startdate = parse(startdate)