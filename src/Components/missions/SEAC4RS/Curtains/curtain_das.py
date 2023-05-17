#!/usr/bin/env python3
"""
Plot curtain from FlightPlan text files.
"""

import sys
from curtain import *

if __name__ == "__main__":

    fname = 'SEAC4RS_130816_s3_er2_1_table_path.txt'
    fname = 'SEAC4RS_130816_s3_dc8_1c_table_path.txt'

    start = 'latest'

    if len(sys.argv) < 2:
        print("Usage: ")
        print("       curtain_fp flight_path.txt ")
        print("Example:") 
        print("       curtain_fp SEAC4RS_130816_s3_dc8_1c_table_path.txt")
        print("       curtain_fp SEAC4RS_130816_s3_dc8_1c_table_path.txt ")
        raise SystemExit("Error: not enough arguments")
    else:
        fname = sys.argv[1]

    # Dataset location
    # ----------------
    root = 'http://opendap.nccs.nasa.gov/dods/GEOS-5/fp/0.25_deg/assim'
    asm_Nv = root+'/inst3_3d_asm_Nv'
    aer_Nv = root+'/inst3_3d_aer_Nv'
    flx_Nx = root+'/tavg1_2d_flx_Nx'

    # Load flight path
    # ----------------
    f = CurtainDC8(fname)

    # Interpolate forecast
    # --------------------
    f.sampleAer(asm_Nv,aer_Nv,flx_Nx)

def hold():
    
    # Plot them
    # ---------
    figure(figsize=(11,5.5))
    figTail = '.%s.%s.png'%(f.aircraft,str(f.Tyme[0]).split()[0])
    f.contourf(1000*f.ext,  'GEOS-5 Total Aerosol Extinction [1/km]',
               N=32,cmap=cm.jet,figFile='ext'+figTail)
    f.contourf(1000*f.duext,'GEOS-5 Dust Aerosol Extinction [1/km]',
               N=32,cmap=cm.jet,figFile='duext'+figTail)
    f.contourf(1000*f.bcext,'GEOS-5 Black Carbon Aerosol Extinction [1/km]',
               N=32,cmap=cm.jet,figFile='bcext'+figTail)
    f.contourf(1000*f.ocext,'GEOS-5 Organic Aerosol Extinction [1/km]',
               N=32,cmap=cm.jet,figFile='ocext'+figTail)
    f.contourf(1000*f.suext,'GEOS-5 Sulfate Aerosol Extinction [1/km]'
               ,N=32,cmap=cm.jet,figFile='suext'+figTail)
