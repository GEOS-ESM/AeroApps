"""
Reduces AVHRR according to avering kernel and speciation.
"""

import os
import sys
from pyobs.avhrr import AVHRR_L2B

from glob import glob
from datetime  import date, datetime, timedelta

if __name__ == "__main__":

    RootDir = '/nobackup/AVHRR/Level2/PATMOSX'
    MAerDir = '/nobackup/MERRAero'
    MDir = '/nobackup/MERRA'

    gas_x = MAerDir + '/inst2d_gaas_x.ddf'
    aer_x = MAerDir + '/tavg2d_aer_x.ddf'
    avk_x = MAerDir + '/inst2d_avk_x.ddf'
    ext_Nc = MAerDir + '/ext_Nc.ddf'
    int_x = MDir + '/int_Nx'
    slv_x = MDir + '/slv_Nx'

    RootDir = '/nobackup/AVHRR/Level2/PATMOSX'

    oneday = timedelta(seconds=24*60*60)

    for year in range(2007,2008):
        for doy in range(118,366):

            fname = 'AVHRR_NPZ/%d/avhrr.%d_%03d.npz'%(year,year,doy)

            if os.path.exists(fname):
                print "<> Skipping ", fname
                continue

            Files = sorted(glob(RootDir+'/%d/%03d/*_asc_*.hdf'%(year,doy)))

            a = AVHRR_L2B(Files,Verb=False)
            a.sampleG5(avk_x=avk_x)

            nobs_ = a.nobs

            # Retain only those obs that correspond to MODIS retrievals
            # --------------------------------------------------------
            I = (a.avk>0.90)

            a.reduce(I)
            a.speciate(aer_x)

            J = (a.fdu>0.60) |\
                (a.fss>0.80) |\
                (a.fcc>0.60) |\
                (a.fsu>0.60) 

            a.reduce(J)

            if a.nobs>10:
                a.sampleG5(gas_x=gas_x,int_x=int_x,slv_x=slv_x)
            else:
                print '%d %3d | %7d %5d | %5d %5d %5d %5d | %5d'%\
                    (year,doy,nobs_,a.nobs,0,0,0,0,0)
                continue

            ndu = len(a.lon[a.fdu>0.60])
            nss = len(a.lon[a.fss>0.80])
            ncc = len(a.lon[a.fcc>0.60])
            nsu = len(a.lon[a.fsu>0.60])

            print '%d %3d | %7d %5d | %5d %5d %5d %5d | %5d'%\
                (year,doy,nobs_,a.nobs,ndu,nss,ncc,nsu,ndu+nss+ncc+nsu)

            # Save NPZ file
            # -------------
            a.writeNPZ(fname)
