#!/usr/bin/env python
"""
Plot curtain from actual aircraft tracks given in ICARTT files.
"""

import sys
import os
from curtain import *
from grads import gacm

if __name__ == "__main__":

    sdir = '/Users/adasilva/iesa/kaq/sampled/flight'
    meteo =  sdir + '/KORUSAQ-GEOS5-METEO-AIRCRAFT_MODEL_DATE_R0.nc'
    chem  =  sdir + '/KORUSAQ-GEOS5-CHEM-AIRCRAFT_MODEL_DATE_R0.nc'
    ext   =  sdir + '/KORUSAQ-GEOS5-EXT532-AIRCRAFT_MODEL_DATE_R0.nc'
    ict = '../Tracks/KORUSAQ-AIRCRAFTTRK_MODEL_DATE_RA.ICT'

    if len(sys.argv) < 3:
        print "Usage: "
        print "       curtain_flight aircraft date"
        print "Example:" 
        print "       curtain_flight DC8 20160503"
        raise SystemExit, "Error: not enough arguments"
    else:
        aircraft = sys.argv[1]
        date = sys.argv[2]

    figTail = '.'+aircraft+'.'+ date + '.png'

    # Dataset location
    # ----------------
    meteo = meteo.replace('AIRCRAFT',aircraft).replace('DATE',date)
    chem = chem.replace('AIRCRAFT',aircraft).replace('DATE',date)
    ext = ext.replace('AIRCRAFT',aircraft).replace('DATE',date)
    ict = ict.replace('AIRCRAFT',aircraft.lower()).replace('DATE',date)

    # Load flight path and sampled data
    # ---------------------------------
    f = Curtain(meteo,chem,ext,ict,aircraft=aircraft,prs=False)

    # Plot them
    # ---------
    figure(figsize=(11,5.5))

    # Pre-load meteorology
    # --------------------
    f.loadMet()
    f.contourf(f.T,r'GEOS-5 Temperature [K]',
               cmap=cm.hot,figFile='tmpu'+figTail, N=32, extend='max')
    f.contourf(f.RH,r'GEOS-5 Relative Humidity',
               cmap=cm.RdBu_r,figFile='rh'+figTail, N=32, extend='max', vmin=0,vmax=100)
    try:
        f.contourf(f.CLOUD,r'GEOS-5 Cloud Fraction [%]',
                   cmap=cm.RdBu,figFile='cld'+figTail, N=32, extend='max')
    except:
        pass

def later():
    
    # Concentration
    # -------------
    f.load(f.fh_chm,['du','ss','bc','oc','so4','so2','co2','co','conbas','cobbgl','cobbae'],
           factor='airdens')
    f.cobbot = f.cobbgl-f.cobbae

    # Plot them
    f.contourf(1e9*f.du,r'GEOS-5 Dust Aerosol Concentration [$\mu$g/m$^3$]    ',
               cmap=gacm.jet_l,figFile='du'+figTail, N=32, extend='max')

    f.contourf(1e9*f.bc,r'GEOS-5 Black Carbon Aerosol Concentration [$\mu$g/m$^3$]',
               cmap=gacm.jet_l,figFile='bc'+figTail,N=32, extend='max')
    f.contourf(1e9*f.oc,r'GEOS-5 Organic Aerosol Concentration [$\mu$g/m$^3$]',
               cmap=gacm.jet_l,figFile='oc'+figTail,N=32, end='max')
    f.contourf(1e9*f.so4,r'GEOS-5 Sulfate Aerosol Concentration [$\mu$g/m$^3$]',
               cmap=gacm.jet_l,figFile='su'+figTail, N=32, extend='max')
    f.contourf(1e9*f.so2,r'GEOS-5 Sulfer Dioxide Concentration [$\mu$g/m$^3$]',
               cmap=gacm.jet_l,figFile='so2'+figTail, N=32, extend='max')

    f.contourf(1e9*f.co2,r'GEOS-5 CO$_2$ Concentration [$\mu$g/m$^3$]',
               cmap=gacm.jet_l,figFile='co2'+figTail, N=32, extend='max')

    f.contourf(1e9*f.co,r'GEOS-5 CO Concentration [$\mu$g/m$^3$]',
               cmap=gacm.jet_l,figFile='co'+figTail, N=32, extend='max')
    f.contourf(1e9*f.conbas,r'GEOS-5 CO FF Asia Concentration [$\mu$g/m$^3$]',
               cmap=gacm.jet_l,figFile='coffas'+figTail, N=32, extend='max')
    f.contourf(1e9*f.cobbae,r'GEOS-5 CO BB N Asia+Europe Concentration [$\mu$g/m$^3$]',
               cmap=gacm.jet_l,figFile='cobbae'+figTail, N=32, extend='max')
    f.contourf(1e9*f.cobbot,r'GEOS-5 CO BB Other Concentration [$\mu$g/m$^3$]',
               cmap=gacm.jet_l,figFile='cobbot'+figTail, N=32, extend='max')

 
    # Extinction
    # ----------
    f.load(f.fh_ext,['duext','ssext','bcext','ocext','suext'])
    f.ext  = f.duext+f.ssext+f.bcext+f.ocext+f.suext

    f.contourf(1000*f.ext,  'GEOS-5 Total Aerosol Extinction [1/km]',
               cmap=gacm.jet_l,figFile='ext'+figTail,N=32, extend='max')
    f.contourf(1000*f.duext,'GEOS-5 Dust Aerosol Extinction [1/km]',
               cmap=gacm.jet_l,figFile='duext'+figTail, N=32, extend='max')
    f.contourf(1000*f.bcext,'GEOS-5 Black Carbon Aerosol Extinction [1/km]',
               cmap=gacm.jet_l,figFile='bcext'+figTail,N=32, extend='max')
    f.contourf(1000*f.ocext,'GEOS-5 Organic Aerosol Extinction [1/km]',
               cmap=gacm.jet_l,figFile='ocext'+figTail,N=32, extend='max')
    f.contourf(1000*f.suext,'GEOS-5 Sulfate Aerosol Extinction [1/km]',
               cmap=gacm.jet_l,figFile='suext'+figTail, N=32, extend='max')

