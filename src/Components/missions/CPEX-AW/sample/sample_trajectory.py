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
from datetime import datetime
from traj_make_plot import get_model_configuration

config = './g2g_pm25.yaml'

def sample(ictFile,model='m21c',collection='aer_inst_3hr_glo_Nv'):

# Get the model configuration
    print(ictFile)
    yyyymmdd, dateout, modname, aircraft, campaign, cs, do_optics, fpdata, config, collname = get_model_configuration(ictFile,model=model,collection=collection)
    m = ICARTT(ictFile)

#   Make an output directory
    dirname = f"samples/{campaign}/sampled/{aircraft}/{modname}/{dateout}"
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
    print("Start trajectory: ",str(datetime.now()))
    alt, lon, lat, tyme = m.Nav['Altitude'], m.Nav['Longitude'], m.Nav['Latitude'], m.Nav['Time']
    traj = TRAJECTORY(tyme,lon,lat,fpdata,cs=cs)
# sample the dataset along the trajectory, and return an xarray dataset
    print("Start sample", str(datetime.now()))
    traj_ds = traj.sample()

#   Rename dimensions
    print("Rename Dimensions", str(datetime.now()))
#   Rename variables
    try:
        traj_ds = traj_ds.rename_vars(lev="level", lon="longitude", lat="latitude")
    except:
        traj_ds = traj_ds.rename_vars(lon="longitude", lat="latitude")
    
#   Add global attributes
    print("Add attributes", str(datetime.now()))
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
                                    project                     =f"{campaign}",
                                    source                      =modname,
                                    title                       =titlestr,
                                    VersionID                   =rnum,
                                    data_product_groups         ="",
                                    data_use_guideline          ="see: https://gmao.gsfc.nasa.gov/geos-systems/",
                                    file_originator             ="Peter Colarco",
                                    file_originator_contact     ="Peter.R.Colarco@nasa.gov",
                                    flight_start_date           =dd,
                                    last_modified_date          =strftime("%Y-%m-%d %H:%M:%S",gmtime()),
                                    measurement_platform        =modname,
                                    platform_identifier         =modname,
                                    time_coverage_end           =str(tyme[-1]),
                                    time_coverage_start         =str(tyme[0]))
                                   

# create a placeholder to put the aircraft altitude data
    if(do_optics):
        print("Start optics", str(datetime.now()))
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
    print("Write file", str(datetime.now()))
    outFile = f"./{dirname}/{campaign}-{modname}-{collname}-{aircraft}_Model_{yyyymmdd}.nc"
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
