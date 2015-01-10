#!/usr/bin/env python

"""
  A Python script to create VLIDORT/OMI Level 3a files.
"""

import warnings
warnings.simplefilter('ignore',DeprecationWarning)
warnings.simplefilter('always',UserWarning)

from omi import OMI
from MAPL       import strTemplate

if __name__ == "__main__":

    ods_template = '/discover/nobackup/projects/gmao/iesa/aerosol/data/OMI/Level2/ODS/Y%y4/M%m2/omi.aero_tc8.obs.%y4%m2%d2.ods'

    for nymd in range(20080701,20080732):
        for nhms in ( 0, 30000, 60000, 90000, 120000, 150000, 180000, 210000 ):
            omi_file = strTemplate(ods_template,nymd=nymd,nhms=nhms)
            try:
                omi = OMI(omi_file,nymd,nhms)    
                print "%d %6d [ok] --- nobs = %d"%(nymd,nhms,omi.nobs)
            except:
                print "%d %6d [not ok]"%(nymd,nhms)

