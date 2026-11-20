#!/usr/bin/env python3
import sys
from pyobs.sampler import TRAJECTORY
from pyobs.icartt import ICARTT
from pyobs.aop import G2GAOP

import numpy as np
from optparse   import OptionParser   # Command-line args
import os
import sys
from time import strftime, gmtime

config = './g2g_pm25.yaml'

def sample(ictFile,model='fp',collection='inst3d_aer_v'):

# Get the model configuration
    do_optics = 1
    fpdata = ['./'+model+'.'+collection+'.ddf','./'+model+'.prog.eta.ddf',
              './'+model+'.inst2d_hwl_x.ddf']
    modname = model
    if(model == 'fp'):
        fpdata = ['./'+model+'.'+collection+'.ddf','./'+model+'.inst3_3d_asm_Nv.ddf',
                  './'+model+'.inst1_2d_hwl_Nx.ddf']
        if(collection == "inst3_3d_aer_Nv"):
            collname = "inst3-3d-AER-Nv"
        if(collection == "tavg1_2d_lfo_Nx"):
            do_optics = 0
            collname = "tavg1-2d-LFO-Nx"
        modname = "GEOS-FP"
        config = './g2g_pm25.yaml'
    if(model == 'MERRA2'):
        fpdata = ['./'+model+'.'+collection+'.ddf','./'+model+'.inst3_3d_asm_Nv.ddf']
        config = './m2_pm25.yaml'
        if(collection == 'tavg1_2d_aer_Nx'):
            fpdata = ['./'+model+'.'+collection+'.ddf','./'+model+'.tavg1_2d_adg_Nx.ddf']
            do_optics = 0

#   Get the ICARTT file describing the trajectory
    if ictFile.find('GV') > 0:
        aircraft = 'GV'
        i0 = ictFile.find('GV')+3
    if ictFile.find('ER2')  > 0:
        aircraft = 'ER2'
        i0 = ictFile.find('ER2')+4
    m = ICARTT(ictFile)
    yyyymmdd  = ict[i0:-4]
    dateout   = ict[i0:i0+4]+"-"+ict[i0+4:i0+6]+"-"+ict[i0+6:i0+8]
    print(ictFile, aircraft, yyyymmdd)

#   Make an output directory
    dirname = f"samples/INSPYRE/sampled/{aircraft}/{modname}/{dateout}"
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
        
# create a trajectory object based on the model field provided
    alt, lon, lat, tyme = m.Nav['Altitude'], m.Nav['Longitude'], m.Nav['Latitude'], m.Nav['Time']
    traj = TRAJECTORY(tyme,lon,lat,fpdata)
# sample the dataset along the trajectory, and return an xarray dataset
    traj_ds = traj.sample()

#   Rename dimensions
    traj_ds = traj_ds.rename_dims(lev="level")
#   Rename variables
    traj_ds = traj_ds.rename_vars(lev="level", lon="longitude", lat="latitude")
    
#   Add global attributes
    dd   = yyyymmdd[0:8]
    rnum = yyyymmdd[9:11]
    fltrck = os.path.basename(ictFile)
    titlestr = f"{modname} model sampled output along {fltrck}"
    url = "https://www-air.larc.nasa.gov/missions/etc/AtmosphericCompositionVariableStandardNames.pdf"
    urlvers = "September 26, 2025"
    keywords = "EARTH SCIENCE, ATMOSPHERE, AEROSOLS, EARTH SCIENCE SERVICES, MODELS, ATMOSPHERIC CHEMISTRY MODELS"
    traj_ds = traj_ds.assign_attrs( ACVSNC_standard_name_URL    =url,
                                    ACVSNC_standard_name_version=urlvers,
                                    Conventions                 ="CF-1.13",
                                    format                      ="netCDF-4",
                                    history                     ="v1.0.0",
                                    institution                 ="Code 614 NASA GSFC",
                                    keywords                    =keywords,
                                    PI_contact                  ="Peter.R.Colarco@nasa.gov",
                                    PI_name                     ="Peter Colarco",
                                    ProcessingLevel             ="L4",
                                    project                     ="INSPYRE",
                                    source                      =modname,
                                    title                       =titlestr,
                                    VersionID                   =rnum,
                                    data_product_groups         ="",
                                    data_use_guideline          ="see: https://gmao.gsfc.nasa.gov/geos-systems/",
                                    file_originator             ="Peter Colarco",
                                    file_originator_contact     ="Peter.R.Colarco@nasa.gov",
                                    flight_start_date           =dd,
                                    last_modified_date          =strftime("%Y-%m-%d %H:%M:%S",gmtime()),
                                    measurement_platform        ="GEOS-FP",
                                    platform_identifier         ="GEOS-FP",
                                    time_coverage_end           =str(tyme[-1]),
                                    time_coverage_start         =str(tyme[0]))
                                   

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
        optics = G2GAOP(traj_ds,config=config)
        ext = optics.getAOPext(wavelength=532,Species=Species)
        traj_ds["EXT532nm"] = ext["EXT"]
        traj_ds["SCA532nm"] = ext["SCA"]
        traj_ds["BSC532nm"] = ext["BSC"]
        traj_ds["DEPOL532nm"] = ext["DEPOL"]
    
# Preferred ICARTT name
    outFile = f"./{dirname}/INSPYRE-{modname}-{collname}-{aircraft}_Model_{yyyymmdd}.nc"
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
