#!/bin/python
# -*- coding: utf-8 -*-
""" Runscript to clean up directories for geo_vlidort runs"""
from datetime import datetime, timedelta 
from dateutil.parser import parse
import os
import shutil
from distutils.dir_util import mkpath
import numpy as np
import math

def destroy_workspace(date,ch,code):
    dirname = str(date.date())+'T'+str(date.hour).zfill(2)+'.'+ch+'.'+code
    if os.path.exists(dirname):
        shutil.rmtree(dirname)

    slurm = dirname + '*.out'
    if os.path.isfile(slurm):
        os.remove(slurm)

            
#########################################################

if __name__ == "__main__":
    
    startdate = '2005-12-31T00:00:00'
    enddate   = '2005-12-31T23:00:00'
    channels  = '470','550','670'
    surface   = 'MAIACRTLS'
    interp    = 'interpolate'
    code      = 'vector'
        
    dt        = timedelta(hours=1)
    startdate = parse(startdate)
    enddate   = parse(enddate)

    if (interp.lower() == 'interpolate'):
        code = code + '.i' 

    code = code + surface

    # destroy working directories 
    while (startdate <= enddate):
        for i, ch in enumerate(channels):

            destroy_workspace(startdate,ch,code)

       
        startdate = startdate + dt

