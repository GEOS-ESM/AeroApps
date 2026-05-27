#!/usr/bin/env python3
import sys
from pyobs.sampler import STATION
from pyobs.aop import G2GAOP
import numpy as np
import datetime
from optparse   import OptionParser   # Command-line args  

config = './g2g_pm25.yaml'

#   Define the stations to use
stations = ['Ragged_Point','Capo_Verde']
lons     = [-59.4,-22.9]
lats     = [13.2,16.7]


def sample(model='c180R_arcsix', collection='inst3d_aer_v', date0='20240501'):

#   Do for one day at a time, 3 hour chunks
    start = datetime.datetime.strptime(date0, "%Y%m%d")
    end   = start + datetime.timedelta(days=1)
    dates = [start + datetime.timedelta(hours=x*3) for x in range(0, (end-start).days*8)]
    tyme     = dates

# create a trajectory object based on the model field provided
    do_optics = 1
    fpdata = ['./'+model+'.'+collection+'.ddf','./'+model+'.prog.eta.ddf',
              './'+model+'.inst2d_hwl_x.ddf']
    config = './g2g_pm25.yaml'
    if(model == 'fp'):
        fpdata = ['./'+model+'.'+collection+'.ddf','./'+model+'.inst3_3d_asm_Nv.ddf',
                  './'+model+'.inst1_2d_hwl_Nx.ddf']
        config = './geos529_pm25.yaml'
    if(model == 'MERRA2'):
        fpdata = ['./'+model+'.'+collection+'.ddf','./'+model+'.inst3_3d_asm_Nv.ddf']
        config = './m2_pm25.yaml'
        if(collection == 'tavg1_2d_aer_Nx'):
            dates  = [start +datetime.timedelta(0,1800) + datetime.timedelta(hours=x*1) for x in range(0, (end-start).days*24)]
            tyme   = dates
            fpdata = ['./'+model+'.'+collection+'.ddf','./'+model+'.tavg1_2d_adg_Nx.ddf']
            do_optics = 0

    stn = STATION(stations,lons,lats,fpdata,time_range=[min(dates),max(dates)],verbose=True)

# sample the dataset along the trajectory, and return an xarray dataset
    stn_ds = stn.sample().compute()
    print(min(dates), max(dates))
    print(stn_ds["TOTEXTTAU"].size)
    
# create a nominal optics profile
    if(do_optics):
        Species = None
        if(collection == 'inst3d_aerdms_v'):
            Species = ['SU']
        optics = G2GAOP(stn_ds,config=config)
        ext = optics.getAOPext(wavelength=532,Species=Species)
        stn_ds["EXT532nm"] = ext["EXT"]
        stn_ds["SCA532nm"] = ext["SCA"]
        stn_ds["BSC532nm"] = ext["BSC"]
        stn_ds["DEPOL532nm"] = ext["DEPOL"]

# write data to a netcdf file
    outFile = './stn_samples/%s.%s'%(model,collection)+'.stations.%s.nc'%(date0)
    print(outFile)
    stn_ds.to_netcdf(path=outFile.format(stations),engine='h5netcdf')

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
        date0      = args[2]
    else:
        parser.error("must have 3 argument: modelname date0")
        
    sample(model=model,collection=collection,date0=date0)
    sys.exit()
