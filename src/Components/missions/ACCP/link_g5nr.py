#!/usr/bin/env python3

"""
    Link surface files for M2R12K runs
"""

import os
import sys
from glob import glob
from datetime import datetime, timedelta


if __name__ == "__main__":

    orbit = 'SS450'
    startdate = datetime(2013,8,1,0)
    enddate   = datetime(2013,9,1,0)
    dt        = timedelta(hours=1)


    rootDir = "/discover/nobackup/projects/gmao/osse2/pub/c1440_NR/OBS/A-CCP/"
 
    # change dirs
    os.chdir(rootDir)


    # surface files
    surface = orbit+'/LevelB/surface/'
    ler     = surface + 'LER/OMI'
    brdf    = surface + 'BRDF/MCD43C1/006/'
    angles  = orbit+'/LevelB/'

    surfaceo = orbit+'-M2R12K/LevelB/surface/'
    lero    = surfaceo + 'LER/OMI'
    brdfo   = surfaceo + 'BRDF/MCD43C1/006/'

    if not os.path.exists(lero):
        os.makedirs(lero)
    if not os.path.exists(brdfo):
        os.makedirs(brdfo)

    dirs = glob(ler+startdate.strftime('/Y2006/M%m/*'))
    for dd in dirs:
        if os.path.isdir(dd):
            files = glob(dd+'/*')
            for ff in files:
                of = ff.replace('2006','2013').replace(orbit,orbit+'-M2R12K')
                od = os.path.dirname(of)
                if not os.path.exists(od):
                    os.makedirs(od)
                if os.path.islink(of):
                    os.unlink(of)
                os.symlink(rootDir+ff,of)
        if os.path.isfile(dd):
            of = dd.replace('2006','2013').replace(orbit,orbit+'-M2R12K')
            od = os.path.dirname(of)
            if not os.path.exists(od):
                os.makedirs(od)
            if os.path.islink(of):
                os.unlink(of)
            os.symlink(rootDir+ff,of)

    dirs = glob(brdf+startdate.strftime('/Y2006/M%m/*'))
    for dd in dirs:
        if os.path.isdir(dd):
            files = glob(dd+'/*')
            for ff in files:
                of = ff.replace('2006','2013').replace(orbit,orbit+'-M2R12K')
                od = os.path.dirname(of)
                if not os.path.exists(od):
                    os.makedirs(od)
                if os.path.islink(of):
                    os.unlink(of)
                os.symlink(rootDir+ff,of)   
        if os.path.isfile(dd):
            of = dd.replace('2006','2013').replace(orbit,orbit+'-M2R12K')
            od = os.path.dirname(of)
            if not os.path.exists(od):
                os.makedirs(od)
            if os.path.islink(of):
                os.unlink(of)
            os.symlink(rootDir+ff,of)

    dirs = glob(angles+startdate.strftime('/Y2006/M%m/*'))
    for dd in dirs:
        if os.path.isdir(dd):
            files = glob(dd+'/*polar07*')
            for ff in files:
                of = ff.replace('2006','2013').replace(orbit,orbit+'-M2R12K').replace('-g5nr','-m2r12k')
                od = os.path.dirname(of)
                if not os.path.exists(od):
                    os.makedirs(od)
                if os.path.islink(of):
                    os.unlink(of)
                os.symlink(rootDir+ff,of)
        if (os.path.isfile(dd) and ('polar07' in dd)):
            of = dd.replace('2006','2013').replace(orbit,orbit+'-M2R12K').replace('-g5nr','-m2r12k')
            od = os.path.dirname(of)
            if not os.path.exists(od):
                os.makedirs(od)
            if os.path.islink(of):
                os.unlink(of)
            os.symlink(rootDir+ff,of)
