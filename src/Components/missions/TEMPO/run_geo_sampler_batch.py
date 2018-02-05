#!/bin/python
# -*- coding: utf-8 -*-
""" Runscript for geo_sampler"""

import os
from datetime import datetime, timedelta 
from dateutil.parser import parse
from distutils.dir_util import mkpath
import numpy as np
from dateutil.relativedelta import relativedelta

instname = 'TEMPO'
nccs     = '/discover/nobackup/projects/gmao/osse2/pub/c1440_NR/OBS/TEMPO/DATA/'
calculon = '/nobackup/{}'.format(instname)

class SCENARIO(object):
    def __init__(self,startdate,enddate,dt,layout,rootDir,varlist):
        self.instname = instname

        self.startdate = parse(startdate)
        self.enddate   = parse(enddate)
        self.dt        = dt
        self.rootDir   = rootDir
        self.layout    = layout
        self.rootDir   = rootDir

        if type(varlist) is str:
            varlist = [varlist]

        self.varlist = varlist

    def sample(self):
        for dataname in self.varlist:
            if dataname == 'SurfLER':
                self.do_SurfLER()
            if dataname == 'aer_Nv':
                self.do_aerNv()
            if dataname == 'MCD43D':
                self.do_MCD43D()


    def do_MCD43D(self):
        if self.layout is not None:
            ntiles = int(self.layout[0])*int(self.layout[1])
            for tile in range(ntiles):
                laycode = self.layout + str(tile)
                
                invariant = '{}/LevelG/invariant/tempo.lg1.invariant.{}.nc4'.format(self.rootDir,laycode)
                outdir    = '{}/BRDF/vMCD43D'.format(self.rootDir)
                date = self.startdate
                while date <= self.enddate: 
                    outfile   = '{}/MCD43D.{}.{}.nc4'.format(outdir,date.strftime('%Y%j'),laycode)
                    landdir   = '{}/LevelB/{}/'.format(self.rootDir,date.strftime('Y%Y/M%m/D%d'))
                    landfile  = '{}/tempo-g5nr-icacl-TOTWPDF-GCOP-SKEWT.{}_{}z.{}.nc4'.format(landdir,date.strftime('%Y%m%d'),date.strftime('%H'),laycode)
                    climdir   = '{}/BRDF/vMCD43C'.format(self.rootDir)
                    climfile  = '{}/tempo.mcd43c1_climatology.{}.nc4'.format(climdir,laycode)
                    command = './mcd43d_sampler.py -v -C -g {} -o {} -l {} -i TEMPO -c {} {}'.format(invariant,outfile,landfile,climfile,date.strftime('%Y%j'))

                    print command
                    os.system(command)

                    date += self.dt
        else:            
            invariant = '{}/LevelG/invariant/tempo.lg1.invariant.nc4'.format(self.rootDir)
            outdir    = '{}/BRDF/vMCD43D'.format(self.rootDir)
            date = self.startdate
            while date <= self.enddate:                
                outfile   = '{}/MCD43D.{}.nc4'.format(outdir,date.strftime('%Y%j'))
                landdir   = '{}/LevelB/{}/'.format(self.rootDir,date.strftime('Y%Y/M%m/D%d'))
                landfile  = '{}/tempo-g5nr-icacl-TOTWPDF-GCOP-SKEWT.{}_{}z.nc4'.format(landdir,date.strftime('%Y%m%d'),date.strftime('%H'))
                climdir   = '{}/BRDF/vMCD43C'.format(self.rootDir)
                climfile  = '{}/tempo.mcd43c1_climatology.nc4'.format(climdir)
                command = './mcd43d_sampler.py -v -C -g {} -o {} -l {} -i TEMPO -c {} {}'.format(invariant,outfile,landfile,climfile,date.strftime('%Y%j'))

                print command
                os.system(command)

                date += self.dt


    def do_aerNv(self):

        invariant = '{}/LevelG/invariant/tempo.lg1.invariant.nc4'.format(self.rootDir)
        rcfile    = '/discover/nobackup/rgovinda/TEMPO/aer_Nv.rc'
        date = self.startdate
        while (date <= self.enddate):
            outdir   = '{}/LevelB/{}/'.format(self.rootDir,date.strftime('Y%Y/M%m/D%d'))
            outfile   = '{}/tempo-g5nr.lb2.aer_Nv.{}_{}z.nc4'.format(outdir,date.strftime('%Y%m%d'),date.strftime('%H'))

            command = './geo_sampler.py -v -C -g {} -o {} -r {} {} {}'.format(invariant,outfile,rcfile,date.isoformat(),date.isoformat())

            print command
            os.system(command)              

            date += self.dt


    def do_SurfLER(self):
        indir = self.rootDir + '/SurfLER'
        outdir = indir
        dt = 1
        date = self.startdate
        while (date <= self.enddate):
            MM = str(date.month).zfill(2)
            command = 'geo_sampler.py -C -N -v -r geo_sampler_omi_ler.rc -o ' + outdir + \
                            '/tempo-omi.SurfLER.' + MM + '.nc4 ' +\
                            date.isoformat()
            print command
            os.system(command)                      
            
            if (date.month < 12):
                date = str(date.year) + '-' + str(date.month+dt).zfill(2) + '-01T00:00:00'
            else:
                date = str(date.year + 1) + '-01-01T00:00:00'

            date = parse(date)    


#########################################################

if __name__ == "__main__":
    
    startdate = '2006-02-15T21:00:00'
    enddate   = '2006-12-15T21:00:00'
    dt        = relativedelta(months=1)
    layout    = '41'
    rootDir   = '/nobackup/3/pcastell/TEMPO/CLD_DATA/'
    varlist   = ['MCD43D']

    scen = SCENARIO(startdate,enddate,dt,layout,rootDir,varlist)
    scen.sample()


    
