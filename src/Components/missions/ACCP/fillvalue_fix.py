#!/usr/bin/env python

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

aer_Nv = ['PS','DELP','RH',
          'DU001','DU002','DU003','DU004','DU005',
          'SS001','SS002','SS003','SS004','SS005',
          'SO2','SO4',
          'BCPHOBIC','BCPHILIC',
          'OCPHOBIC','OCPHILIC']

aer_Nv = ['AIRDENS']

met_Nv = ['TROPPB','TROPT','TROPQ','IWP','LWP',
          'U10M','V10M','PS','TS','PBLH',
          'CLDLOW','CLDMID','CLDHGH','CLDTOT',
          'TAULOW','TAUMID','TAUHGH','TAUTOT',
          'FRSNO','FRSEAICE','T','U','V','QV']

asm_Nx = ['AREA','FRLAKE','FRLAND','FRLANDICE',
          'FROCEAN','PHIS','SGH']

chm_Nv = ['CO','CO2','O3','SO2']

COLS = {'aer_Nv': aer_Nv}
#        'met_Nv': met_Nv, 
#        'asm_Nx': asm_Nx, 
#        'chm_Nv': chm_Nv}


if __name__ == '__main__':

    #Defaults
    DT_hours = 1
    rootDir     = '/discover/nobackup/projects/gmao/osse2/pub/c1440_NR/OBS/A-CCP/'
#    rootDir     = '/nobackup/3/pcastell/A-CCP'

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
            inFile = '{}/{}/LevelB/Y{}/M{}/D{}/{}-g5nr.lb2.{}.{}_{}00z.nc4'.format(rootDir,orbit.upper(),Y,M,D,orbit,col,ymd,hh)
            nc = Dataset(inFile,'r+')

            VARS = COLS[col]

            for var in VARS:
                varobj = nc.variables[var]
                undef = varobj.missing_value
                uu    = undef.astype(np.float32)
                varobj.missing_value = uu

            nc.close()

        sdate += dt
