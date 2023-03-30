#!/usr/bin/env python3
"""
Plot curtain from actual aircraft tracks given in ICARTT files.
"""

import sys
import os
from curtain import *
from grads import gacm

if __name__ == "__main__":

    sdir = '/home/adasilva/iesa/aerosol/data/WE-CAN/sampled'
    nav = '/home/adasilva/iesa/aerosol/data/WE-CAN/Nav'
    meteo =  sdir + '/WECANrfNN-METEO_R0.nc'
    chem  =  sdir + '/WECANrfNN-CHEM_R0.nc'
    ext   =  sdir + '/WECANrfNN-EXT355_R0.nc'
    ict   =  nav + '/WECANrfNN.ict.gz'

    if len(sys.argv) < 2:
        print("Usage: ")
        print("       curtain_flight NN")
        print("Example:") 
        print("       curtain_flight 11")
        raise SystemExit("Error: not enough arguments")
    else:
        rf = sys.argv[1]

    figTail = '.rf' + rf + '.png'

    # Dataset location
    # ----------------
    meteo = meteo.replace('NN',rf)
    chem = chem.replace('NN',rf)
    ext = ext.replace('NN',rf)
    ict = ict.replace('NN',rf)

    # Load flight path and sampled data
    # ---------------------------------
    f = Curtain(meteo,chem,ext,ict,aircraft='C130',prs=False)

    # Plot them
    # ---------
    figure(figsize=(11,5.5))

    # Pre-load meteorology
    # --------------------
    f.loadMet()
    f.contourf(f.T,r'GEOS-5 Temperature [K]',
               cmap=cm.hot,figFile='WECAN.tmpu'+figTail, N=32, extend='max')
    f.contourf(f.RH,r'GEOS-5 Relative Humidity',
               cmap=cm.RdBu_r,figFile='WECAN.rh'+figTail, N=32, extend='max', vmin=0,vmax=100)
    #f.contourf(100*f.CLOUD,r'GEOS-5 Cloud Fraction [%]',
    #           cmap=gacm.Bone_l,figFile='WECAN.cld'+figTail,pblc='w', N=32, extend='max')

    # Concentration
    # -------------
    f.load(f.fh_chm,['du001','du002','du003','du004','du005',
                     'ocphilic','ocphobic','bcphilic','bcphobic',
                     'so4','so2','co2','co'], #,'conbas','cobbgl','cobbae'],
           factor='airdens')
    #f.cobbot = f.cobbgl-f.cobbae
    
    # Plot them
    f.contourf(1e9*(f.du001+f.du002+f.du003+f.du004+f.du005),
               r'GEOS-5 Dust Aerosol Concentration [$\mu$g/m$^3$]    ',
               cmap=gacm.jet_l,figFile='WECAN.du'+figTail, N=32, extend='max')
    f.contourf(1e9*(f.bcphobic+f.bcphilic),
               r'GEOS-5 Black Carbon Aerosol Concentration [$\mu$g/m$^3$]',
               cmap=gacm.jet_l,figFile='WECAN.bc'+figTail,N=32, extend='max')
    f.contourf(1e9*(f.ocphobic+f.ocphilic),
               r'GEOS-5 Organic Aerosol Concentration [$\mu$g/m$^3$]',
               cmap=gacm.jet_l,figFile='WECAN.oc'+figTail,N=32, end='max')
    f.contourf(1e9*f.so4,r'GEOS-5 Sulfate Aerosol Concentration [$\mu$g/m$^3$]',
               cmap=gacm.jet_l,figFile='WECAN.su'+figTail, N=32, extend='max')
    f.contourf(1e9*f.so2,r'GEOS-5 Sulfer Dioxide Concentration [$\mu$g/m$^3$]',
               cmap=gacm.jet_l,figFile='WECAN.so2'+figTail, N=32, extend='max')

    f.contourf(1e9*f.co2,r'GEOS-5 CO$_2$ Concentration [$\mu$g/m$^3$]',
               cmap=gacm.jet_l,figFile='WECAN.co2'+figTail, N=32, extend='max')

    f.contourf(1e9*f.co,r'GEOS-5 CO Concentration [$\mu$g/m$^3$]',
               cmap=gacm.jet_l,figFile='WECAN.co'+figTail, N=32, extend='max')

    """
    f.contourf(1e9*f.conbas,r'GEOS-5 CO FF Asia Concentration [$\mu$g/m$^3$]',
               cmap=gacm.jet_l,figFile='WECAN.coffas'+figTail, N=32, extend='max')
    f.contourf(1e9*f.cobbae,r'GEOS-5 CO BB N Asia+Europe Concentration [$\mu$g/m$^3$]',
               cmap=gacm.jet_l,figFile='WECAN.cobbae'+figTail, N=32, extend='max')
    f.contourf(1e9*f.cobbot,r'GEOS-5 CO BB Other Concentration [$\mu$g/m$^3$]',
               cmap=gacm.jet_l,figFile='WECAN.cobbot'+figTail, N=32, extend='max')
    """

    # Extinction
    # ----------
    f.load(f.fh_ext,['ext','backscat','depol','ext2back'])
    f.contourf(1000*f.ext,  'GEOS-5 Aerosol 355nm Extinction [1/km]',
               cmap=gacm.jet_l,figFile='WECAN.ext'+figTail,N=32, extend='max')
    f.contourf(1000*f.backscat,'GEOS-5 Aerosol 355nm Backscatter [1/km]',
               cmap=gacm.jet_l,figFile='WECAN.bacscat'+figTail, N=32, extend='max')
    f.contourf(f.depol,'GEOS-5 Aerosol Depolarization',
               cmap=gacm.jet_l,figFile='WECAN.depol'+figTail,N=32, extend='max')
    f.contourf(f.ext2back,'GEOS-5 Extinction:Backscatter Ratio',
               cmap=gacm.jet_l,figFile='WECAN.ext2back'+figTail,N=32, extend='max')

