#!/usr/bin/env python
"""
Plot curtain from FlightPlan text files.
"""

import sys
import os
from curtain import *
from grads import gacm

if __name__ == "__main__":

    sdir = '/nobackup/3/pcastell/iesa/kaq/sampled/plan'
    meteo =  sdir + '/KORUSAQ-GEOS5-METEO-AIRCRAFT_PLAN_DATE_R0.nc'
    chem  =  sdir + '/KORUSAQ-GEOS5-CHEM-AIRCRAFT_PLAN_DATE_R0.nc'
    ict = '../Plans/fltplan_aircraft_DATE.ict'
    csv = '../Plans/xsect_SECTION_isoT0.csv'

    if len(sys.argv) < 3:
        print "Usage: "
        print "       curtain_plan aircraft date"
        print "       curtain_plan section  iso_time"
        print "Example:" 
        print "       curtain_plan DC8 20160503"
        print "       curtain_plan jeju2seoul 2016-05-24T00:00:00"
        raise SystemExit, "Error: not enough arguments"
    else:
        aircraft = sys.argv[1]
        date = sys.argv[2]

    figTail = '.plan.'+aircraft+'.'+ date + '.png'

    # Dataset location
    # ----------------
    meteo = meteo.replace('AIRCRAFT',aircraft).replace('DATE',date)
    chem = chem.replace('AIRCRAFT',aircraft).replace('DATE',date)
    if 'T' in date:
        coords = csv.replace('SECTION',aircraft).replace('isoT0',date)
    else:
        coords = ict.replace('aircraft',aircraft.lower()).replace('DATE',date)

    # Load flight path and sampled data
    # ---------------------------------
    f = Curtain(meteo,chem,chem,coords,aircraft=aircraft)

    # Plot them
    # ---------
    figure(figsize=(11,5.5))
    clevs=linspace(0,0.2,32)
    clevs2 = linspace(0,0.7,32)

    # Pre-load meteorology
    # --------------------
    f.loadMet()
    f.T.mask[abs(f.T)>400.] = True # detect undef contaminated interp
    f.contourf(f.T,r'GEOS-5 Temperature [K]',
               cmap=cm.hot,figFile='tmpu'+figTail, N=32, extend='max')
    f.RH.mask[abs(f.RH)>400.] = True # detect undef contaminated interp
    f.contourf(f.RH,r'GEOS-5 Relative Humidity',
               cmap=cm.RdBu_r,figFile='rh'+figTail, N=32, extend='max', vmin=0,vmax=100)
    try:
        f.CLOUD.mask[abs(f.CLOUD)>400.] = True # detect undef contaminated interp
        f.contourf(f.CLOUD,r'GEOS-5 Cloud Fraction [%]',
                   cmap=cm.RdBu,figFile='cld'+figTail, N=32, extend='max')
    except:
        pass
        
    # Concentration
    # -------------
    f.load(f.fh_chm,['du','ss','bc','oc','so4','so2','co2','co','conbas','cobbgl','cobbae'],
           factor='airdens')
    f.cobbot = f.cobbgl-f.cobbae

    for spec in ['du','ss','bc','oc','so4','so2','co2','co','conbas','cobbgl','cobbae','cobbot']:
        f.__dict__[spec].mask[abs(f.__dict__[spec])>1e10] = True    # detect undef contaminated interp

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
    for spec in ['duext','ssext','bcext','ocext','suext']:
        f.__dict__[spec].mask[abs(f.__dict__[spec])>1e10] = True    # detect undef contaminated interp

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

