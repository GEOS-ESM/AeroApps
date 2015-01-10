#!/usr/bin/env python
"""
Plot curtain from FlightPlan text files.
"""

import sys
import os
from curtain import *
from grads import gacm

if __name__ == "__main__":

    fname = 'SEAC4RS_130816_s3_er2_1_table_path.txt'
    fname = 'SEAC4RS_130816_s3_dc8_1c_table_path.txt'
    
    start = 'latest'

    if len(sys.argv) < 2:
        print "Usage: "
        print "       curtain_fp flight_path.txt [fcast_start]"
        print "Example:" 
        print "       curtain_fp SEAC4RS_130816_s3_dc8_1c_table_path.txt latest"
        print "       curtain_fp SEAC4RS_130816_s3_dc8_1c_table_path.txt 20130814_12"
        raise SystemExit, "Error: not enough arguments"
    else:
        fname = sys.argv[1]
        if len(sys.argv)>=3:
            start = sys.argv[2]

    figTail = '.'+os.path.basename(fname).replace('SEAC4RS_','').replace('_table_path.txt','')+'.png'
    # figTail = '.%s.%s.png'%(f.aircraft,str(f.Tyme[0]).split()[0])
     
    # Dataset location
    # ----------------
    root = 'http://opendap.nccs.nasa.gov/dods/GEOS-5/fp/0.25_deg/fcast'
    if start == 'latest':
        asm_Np = root+'/inst3_3d_asm_Np.latest'
        asm_Nv = root+'/inst3_3d_asm_Nv.latest'
        ext_Np = root+'/inst3_3d_ext_Np.latest'
        flx_Nx = root+'/tavg1_2d_flx_Nx.latest'
    else:
        asm_Np = root+'/inst3_3d_asm_Np/inst3_3d_asm_Np.'+start
        asm_Nv = root+'/inst3_3d_asm_Nv/inst3_3d_asm_Nv.'+start
        ext_Np = root+'/inst3_3d_ext_Np/inst3_3d_ext_Np.'+start
        flx_Nx = root+'/tavg1_2d_flx_Nx/tavg1_2d_flx_Nx.'+start

    # Load flight path
    # ----------------
    f = CurtainFP(fname)
      
    # Interpolate forecast
    # --------------------
    f.sampleExt(asm_Np,asm_Nv,ext_Np,flx_Nx)
    
    # Plot them
    # ---------
    figure(figsize=(11,5.5))
    clevs=linspace(0,0.2,32)
    clevs2 = linspace(0,0.7,32)

    # Either specify N = 32 or clevs = linspace(min,max,32)
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
