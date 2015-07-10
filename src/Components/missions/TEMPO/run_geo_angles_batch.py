#!/bin/python
# -*- coding: utf-8 -*-
""" Runscript for geo_angles"""
from datetime import datetime, timedelta 
from dateutil.parser import parse
import os
import subprocess
from distutils.dir_util import mkpath
import numpy as np
import math



def make_rcfile(indir,outdir,date,instname,):

    rcfile = open('geo_angles.rc','w')
    rcfile.write('INDIR: ' + indir + '\n')
    rcfile.write('OUTDIR: ' + outdir + '\n')
    rcfile.write('DATE: '+str(date.year)+str(date.month).zfill(2)+str(date.day).zfill(2)+'\n')
    rcfile.write('TIME: '+str(date.hour).zfill(2)+'\n')
    rcfile.write('INSTNAME: ' + instname + '\n')
    rcfile.close()    
        
#########################################################

if __name__ == "__main__":
    
    startdate   = '2007-04-10T00:00:00'
    enddate     = '2007-04-11T23:00:00'
    instname    = 'tempo'
    calculon    = '/nobackup/TEMPO'
    
    indir     = calculon 
    outdir    = calculon
    
    dt = timedelta(hours=1)
    startdate = parse(startdate)
    enddate   = parse(enddate)

    # Create RC-files

    while (startdate <= enddate):
        make_rcfile(indir,outdir,startdate,instname)

        runstring = './geo_angles.x geo_angles.rc'
        os.system(runstring)

        startdate = startdate + dt

    


