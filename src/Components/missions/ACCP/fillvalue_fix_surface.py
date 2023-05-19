#!/usr/bin/env python3

"""
    Wrapper to submit jobs to sbatch
"""

import os
import subprocess
from   datetime        import datetime, timedelta
from   dateutil.parser import parse         as isoparser
import argparse
import numpy as np
from netCDF4 import Dataset

ler = ['SRFLER388','SRFLER354']

ndvi = ['NDVI']


COLS = {'ler': ler, 
        'ndvi': ndvi} 


if __name__ == '__main__':

    #Defaults
    DT_hours = 1
    rootDir     = '/discover/nobackup/projects/gmao/osse2/pub/c1440_NR/OBS/A-CCP/'
    rootDir     = '/nobackup/3/pcastell/A-CCP'

    parser = argparse.ArgumentParser()
    parser.add_argument("iso_t1",help='starting iso time')
    parser.add_argument("iso_t2",help='ending iso time')
    parser.add_argument("orbit",help='orbit name')

    args = parser.parse_args()


    sdate = isoparser(args.iso_t1)
    edate = isoparser(args.iso_t2)
    orbit = args.orbit
    dt    = timedelta(hours=DT_hours)
    while sdate < edate:
        Y = str(sdate.year)
        M = str(sdate.month).zfill(2)
        D = str(sdate.day).zfill(2)
        ymd = Y+M+D
        hh = str(sdate.hour).zfill(2)
        for col in COLS:
            if col == 'ndvi':
                inFile = '{}/{}/LevelB/surface/NDVI/MYD13C2/006/Y{}/M{}/D{}/{}-g5nr.lb2.{}.{}_{}00z.nc4'.format(rootDir,orbit.upper(),Y,M,D,orbit,col,ymd,hh)
            elif col == 'ler':
                inFile = '{}/{}/LevelB/surface/LER/OMI/Y{}/M{}/D{}/{}-g5nr.lb2.{}.{}_{}00z.nc4'.format(rootDir,orbit.upper(),Y,M,D,orbit,col,ymd,hh)

            if os.path.exists(inFile):
                nc = Dataset(inFile,'r+')

                VARS = COLS[col]

                for var in VARS:
                    varobj = nc.variables[var]
                    undef = varobj.missing_value
                    uu    = undef.astype(np.float32)
                    varobj.missing_value = uu

                nc.close()

        sdate += dt
