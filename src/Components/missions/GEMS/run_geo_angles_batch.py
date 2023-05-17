#!/bin/python3
# -*- coding: utf-8 -*-
""" Runscript for geo_angles"""
from datetime import datetime, timedelta 
from dateutil.parser import parse
import os
import subprocess
from distutils.dir_util import mkpath
import numpy as np
import math



def make_rcfile(indir,outdir,date,instname,layout=None):

    dirname = outdir+'/LevelB/Y'+str(date.year)+'/M'+str(date.month).zfill(2)+'/D'+str(date.day).zfill(2)
    if not os.path.exists(dirname):
        os.makedirs(dirname)

    rcfile = open('geo_angles.rc','w')
    rcfile.write('INDIR: ' + indir + '\n')
    rcfile.write('OUTDIR: ' + outdir + '\n')
    rcfile.write('DATE: '+str(date.year)+str(date.month).zfill(2)+str(date.day).zfill(2)+'\n')
    rcfile.write('TIME: '+str(date.hour).zfill(2)+'\n')
    rcfile.write('INSTNAME: ' + instname + '\n')
    if layout is not None:
        rcfile.write('LAYOUT: ' + layout + '\n')
    rcfile.close()    
        
#########################################################

if __name__ == "__main__":
    
    startdate   = '2005-12-31T00:00:00'
    enddate     = '2005-12-31T23:00:00'
    layout      = None
    episode     = None
    instname    = 'gems'
    calculon    = '/nobackup/GEMS'

    if type(episode) is int: episode = [episode]
    if episode is None: episode = [None]

    for e in episode:
        if e == 1:
            startdate  = '2005-12-31T00:00:00Z'
            enddate    = '2006-01-01T23:00:00Z'
        elif e == 2:
            startdate  = '2006-07-27T00:00:00Z'
            enddate    = '2006-08-09T23:00:00Z'
        elif e == 3:
            startdate  = '2007-03-28T00:00:00Z'
            enddate    = '2007-03-29T23:00:00Z'
        elif e == 4:
            startdate  = '2007-04-10T00:00:00Z'
            enddate    = '2007-04-11T23:00:00Z'
    
        indir     = calculon 
        outdir    = calculon
        
        dt = timedelta(hours=1)
        startdate = parse(startdate)
        enddate   = parse(enddate)

        # Create RC-files

        while (startdate <= enddate):
            make_rcfile(indir,outdir,startdate,instname,layout=layout)

            runstring = 'mpirun -np 8 ./geo_angles.x geo_angles.rc'
            os.system(runstring)

            startdate = startdate + dt

    


