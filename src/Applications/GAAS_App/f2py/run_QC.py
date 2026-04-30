#!/usr/bin/env python

import sys
import os

sys.path.append("/usr/share/lmod/lmod/init/")
sys.path.append('/home/vbuchard/workspace/JEDI/TEST/latest_code_tested/BUDDY_CHECK/AeroApps/install/lib/Python')

# Ensure LD_LIBRARY_PATH is set if needed (note: os.environ.get just reads it, 
# if you need to set it from within python you'd use os.environ['LD_LIBRARY_PATH'] = ...)
# os.environ['OMP_NUM_THREADS'] = '20'

from QC_generic import QC 

# input file
iodafile = '/home/vbuchard/workspace/JEDI/OBS_IODA/BUDDY_CHECK/IODA_nnr_postods_VIIRS_exp.20190710_1200z.nc4'


print(f"Starting QC process for {iodafile}...")

qc = QC(
    iodaFiles=iodafile,
    var_name='aerosolOpticalDepth',
    z_coord_name='obs_wavelength',
    log=True,
)

print(f"Done!")
