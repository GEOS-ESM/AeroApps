#!/usr/bin/env python3
import sys
from pyobs.sampler import TRAJECTORY
from pyobs.icartt import ICARTT
from pyobs.aop import G2GAOP

import numpy as np
from optparse   import OptionParser   # Command-line args  

config = './g2g_pm25.yaml'

def sample(ictFile,model='c180R_arcsix',collection='inst3d_aer_v'):

#   Get the ICARTT file describing the trajectory
    if ictFile.find('P3B') > 0:
        aircraft = 'P3B'
        i0 = ictFile.find('P3B')+4
    if ictFile.find('G3')  > 0:
        aircraft = 'G3'
        i0 = ictFile.find('G3')+3
    if ictFile.find('Learjet') > 0:
        aircraft = 'Learjet'
        i0 = ictFile.find('Learjet')+8
    m = ICARTT(ictFile)
    yyyymmdd = ict[i0:i0+8] 
    print(ictFile, aircraft, yyyymmdd)
    alt, lon, lat, tyme = m.Nav['Altitude'], m.Nav['Longitude'], m.Nav['Latitude'], m.Nav['Time']

# create a trajectory object based on the model field provided
    do_optics = 1
    fpdata = ['./'+model+'.'+collection+'.ddf','./'+model+'.prog.eta.ddf',
              './'+model+'.inst2d_hwl_x.ddf']
    if(model == 'fp'):
        fpdata = ['./'+model+'.'+collection+'.ddf','./'+model+'.inst3_3d_asm_Nv.ddf',
                  './'+model+'.inst1_2d_hwl_Nx.ddf']
        config = './m2_pm25.yaml'
    if(model == 'MERRA2'):
        fpdata = ['./'+model+'.'+collection+'.ddf','./'+model+'.inst3_3d_asm_Nv.ddf']
        config = './m2_pm25.yaml'
        if(collection == 'tavg1_2d_aer_Nx'):
            fpdata = ['./'+model+'.'+collection+'.ddf','./'+model+'.tavg1_2d_adg_Nx.ddf']
            do_optics = 0

    traj = TRAJECTORY(tyme,lon,lat,fpdata)
# sample the dataset along the trajectory, and return an xarray dataset
    traj_ds = traj.sample()

# create a placeholder to put the aircraft altitude data
    if(do_optics):
        traj_ds["GPS_ALTITUDE"] = traj_ds["PS"]  # use PS as a template
        attrs = {'long_name':'aircraft_gps_altitude', 'standard_name':'gps_altitude',
                 'short_name':'gps_altitude', 'units':'m'}
        traj_ds["GPS_ALTITUDE"].attrs.update(attrs)
        traj_ds["GPS_ALTITUDE"].values = alt

# create a nominal optics profile
    if(do_optics):
        Species = None
        if(collection == 'inst3d_aerdms_v'):
            Species = ['SU']
        optics = G2GAOP(traj_ds,config=config)
        ext = optics.getAOPext(wavelength=532,Species=Species)
        traj_ds["EXT532nm"] = ext["EXT"]
        traj_ds["SCA532nm"] = ext["SCA"]
        traj_ds["BSC532nm"] = ext["BSC"]
        traj_ds["DEPOL532nm"] = ext["DEPOL"]
    
# write data to a netcdf file
    outFile = './trj_samples/%s.%s.%s.'%(model,collection,aircraft)+yyyymmdd+'.nc'
    print(outFile)
    traj_ds.to_netcdf(outFile)

if __name__ == "__main__":
#-------------------------------------------------------
#   Parse the commandline
#-------------------------------------------------------

#   CHECK INPUT ARGUMENTS
    parser = OptionParser(usage="Usage: %prog [options] modelname date0",
                          version='xxx' )
    (options, args) = parser.parse_args()
 
#  GET OMI FILE FROM INPUT ARGUMENT LIST
    if len(args) == 3:
        model      = args[0]
        collection = args[1]
        ict        = args[2]
    else:
        parser.error("must have 3 argument: modelname date0")
        
    sample(ict,model=model,collection=collection)
    sys.exit()
