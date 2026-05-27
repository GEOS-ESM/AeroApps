#!/usr/bin/env python3
import sys
from pyobs.sampler import STATION
from pyobs.aop import G2GAOP
import numpy as np
import datetime
import os
from optparse   import OptionParser   # Command-line args  
from time import strftime, gmtime

config = './g2g_pm25.yaml'

#   Define the stations to use
stations = ['Resolute_Bay', 'Pituffik', 'Kaffeklubben_Island', 'Eureka',
            'Alert', 'Baffin_Bay', 'Belle', 'Fram_Strait', 'Lisee',
            'Northern_Arctic_Ocean', 'Station_Nord', 'Summit_Greenland',
            'Svalbard_Zeppelin']
lons     = [-94.97, -68.77, -31.2, -86.42, -62.51, -67.4, 87., 0., 28.0,
            -60., -16.66, -38.46, 11.88]
lats     = [74.71, 76.52, 83.62, 80.05, 82.45, 74.6, 87., 78.92, 87.,
            88., 81.6, 72.58, 78.9]


def sample(model='c180R_arcsix', collection='inst3d_aer_v', date0='20240501'):

#   Do for one day at a time, 3 hour chunks
    start = datetime.datetime.strptime(date0, "%Y%m%d")
    end   = start + datetime.timedelta(days=1)
    dates = [start + datetime.timedelta(hours=x*3) for x in range(0, (end-start).days*8)]
    tyme     = dates

# create a trajectory object based on the model field provided
    do_optics = 1
    do_level  = 1
    fpdata = ['./'+model+'.'+collection+'.ddf','./'+model+'.prog.eta.ddf',
              './'+model+'.inst2d_hwl_x.ddf']
    config = './g2g_pm25.yaml'
    modname = model
    if(model == 'fp'):
        fpdata = ['./'+model+'.'+collection+'.ddf','./'+model+'.inst3_3d_asm_Nv.ddf',
                  './'+model+'.inst1_2d_hwl_Nx.ddf']
        config = './geos529_pm25.yaml'
        if(collection == "inst3_3d_aer_Nv"):
            collname = "inst3-3d-AER-Nv"
        if(collection == "tavg1_2d_lfo_Nx"):
            fpdata = ['./'+model+'.'+collection+'.ddf']
            do_optics = 0
            do_level  = 0
            collname = "tavg1-2d-LFO-Nx"
        modname = "GEOS-FP"
    if(model == 'MERRA2'):
        fpdata = ['./'+model+'.'+collection+'.ddf','./'+model+'.inst3_3d_asm_Nv.ddf']
        config = './m2_pm25.yaml'
        modname = "MERRA2"
        if(collection == 'tavg1_2d_aer_Nx'):
            dates  = [start + datetime.timedelta(hours=x*1) for x in range(0, (end-start).days*24)]
            tyme   = dates
            fpdata = ['./'+model+'.'+collection+'.ddf','./'+model+'.tavg1_2d_adg_Nx.ddf']
            do_optics = 0

    stn = STATION(stations,lons,lats,fpdata,time_range=[min(dates),max(dates)],verbose=True)
    stn_ds = stn.sample()
#   Rename dimensions
    if(do_level):
        stn_ds = stn_ds.rename_dims(lev="level")
#   Rename variables
    if(do_level):
        stn_ds = stn_ds.rename_vars(lev="level", lon="longitude", lat="latitude")
    else:
        stn_ds = stn_ds.rename_vars(             lon="longitude", lat="latitude")
    
#   Add global attributes
    titlestr = f"{modname} model sampled at stations"
    url = "https://www-air.larc.nasa.gov/missions/etc/AtmosphericCompositionVariableStandardNames.pdf"
    urlvers = "September 26, 2025"
    keywords = "EARTH SCIENCE, ATMOSPHERE, AEROSOLS, EARTH SCIENCE SERVICES, MODELS, ATMOSPHERIC CHEMISTRY MODELS"
    stn_ds = stn_ds.assign_attrs(   ACVSNC_standard_name_URL    =url,
                                    ACVSNC_standard_name_version=urlvers,
                                    Conventions                 ="CF-1.13",
                                    format                      ="netCDF-4",
                                    history                     ="v1.0.0",
                                    institution                 ="Code 614 NASA GSFC",
                                    keywords                    =keywords,
                                    PI_contact                  ="Peter.R.Colarco@nasa.gov",
                                    PI_name                     ="Peter Colarco",
                                    ProcessingLevel             ="L4",
                                    project                     ="ARCSIX",
                                    source                      =modname,
                                    title                       =titlestr,
                                    VersionID                   ="R01",
                                    data_product_groups         ="",
                                    data_use_guideline          ="see: https://gmao.gsfc.nasa.gov/geos-systems/",
                                    file_originator             ="Peter Colarco",
                                    file_originator_contact     ="Peter.R.Colarco@nasa.gov",
                                    flight_start_date           =date0,
                                    last_modified_date          =strftime("%Y-%m-%d %H:%M:%S",gmtime()),
                                    measurement_platform        ="GEOS-FP",
                                    platform_identifier         ="GEOS-FP",
                                    time_coverage_end           =str(tyme[-1]),
                                    time_coverage_start         =str(tyme[0]))



# sample the dataset along the trajectory, and return an xarray dataset
    stn_ds = stn_ds.compute()

#   Make an output directory
    dirname = f"samples/ARCSIX/sampled/stations/{modname}"
    print(dirname)
    try:
        os.makedirs(dirname)
        print(f"Directory: {dirname} -- created")
    except FileExistsError:
        print(f"Directory: {dirname} -- already exists")
    except PermissionError:
        print(f"Permission denied: Unable to create '{dirname}'.")
    except Exception as e:
        print(f"An error occurred: {e}")
        sys.exit()
        
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
    outFile = f"./{dirname}/ARCSIX-{modname}-{collname}-stations_Model_{date0}.nc"
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
