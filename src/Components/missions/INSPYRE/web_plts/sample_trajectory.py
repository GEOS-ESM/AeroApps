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
import xarray as xr
import pyobs.xrctl  as xc
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.path as mpath
from matplotlib import cm #cm is colormap
import matplotlib.dates as mdates
import matplotlib.colors as colors
import matplotlib.ticker as mticker
import matplotlib
matplotlib.use('agg')


config = './g2g_pm25.yaml'

#Get the HALO RGB colors
import csv
from matplotlib.colors import LinearSegmentedColormap
with open('/home/pcolarco/lib/halo_colorbar.csv', newline='') as csvfile:
    reader = csv.DictReader(csvfile, delimiter=',',fieldnames=['r','g','b'])
    i = 0
    rgb = []
    for row in reader:
        rgb.append([int(row['r'])/255.,int(row['g'])/255.,int(row['b'])/255.])
cm = LinearSegmentedColormap.from_list('my_map', rgb)


def plot_curtain(fname,trj=-100.,model="fp",species=None):

    modeltitle='GEOS-FP'
    modname = model
    if(model == 'fp'):
        modname = "GEOS-FP"
        config = './geos543_pm25.yaml'
    elif(model == 'MERRA2'):
        config = './m2_pm25.yaml'
    else:
        print("STOP")
        sys.exit()            

#   From the filename get the forecast init and valid time
    f0dateout = fname[-33:-22]  # forecast init time
    fvdateout = fname[-21:-10]  # forecast valid time
    datedir   = fname[-33:-25]+"T"+fname[-24:-22]
    if trj > 0:
        lon  = np.arange(-120.,-79.25,0.25)
        lat  = np.full(len(lon),trj)
        alt  = np.full(len(lon),0.)
        tyme = np.full(len(lon),datetime.strptime(fvdateout,"%Y%m%d_%H"))
        trjstr = f"lat0{abs(int(trj)):0d}N"
        trjtstr= f"{abs(int(trj)):0d}N"
        xstr = "Longitude"
    else:
        lat  = np.arange(30.,70.25,0.25)
        lon  = np.full(len(lat),trj)
        alt  = np.full(len(lon),0.)
        tyme = np.full(len(lon),datetime.strptime(fvdateout,"%Y%m%d_%H"))
        trjstr = f"lon{abs(int(trj)):0d}W"
        trjtstr= f"{abs(int(trj)):0d}W"
        xstr = "Latitude"
    dd   = f"{f0dateout}+{fvdateout}"
    dirname = f"samples/INSPYRE/sampled/{modname}/{datedir}"
    samplefile = f"./{dirname}/INSPYRE-{modname}_Model_{dd}.traj_{trjstr}.nc"
    print(samplefile)

    Species = None
    speciestitle = 'Total'
    if(species == 'du'):
        Species = ['DU']
        speciestitle = 'Dust '
    if(species == 'ss'):
        Species = ['SS']
        speciestitle = 'Sea Salt '
    if(species == 'su'):
        Species = ['SU']
        speciestitle = 'Sulfate '
    if(species == 'cc'):
        Species = ['OC','BR','BC']
        if(model == 'res'):
            Species = ['OC','BR','BC']
        speciestitle = 'Carbonaceous '

    optics = G2GAOP(samplefile,config=config)
    asmfile = samplefile
    if isinstance(asmfile,xr.Dataset):
        asm = asmfile
    else:
        asm = xc.open_mfdataset(asmfile)
    z   = asm['H']
#   Get the extinction profile
    ext = optics.getAOPext(wavelength=532,Species=Species)
    if trj < 0:
        x = ext.latitude.values
    else:
        x = ext.longitude.values
    
    fig, ax = plt.subplots(figsize=(18, 8))
    plt.subplots_adjust(left=0.1,bottom=0.1,right=0.99,top=0.9)
    ax.set_facecolor("white")
    nx = len(x)
    nlev = ext.sizes['level']
    xx = np.repeat(x.reshape(nx,1),nlev,axis=1)
    plt.ylabel('Altitude [km]')
    plt.xlabel(xstr)
    clevs = np.arange(-2,0.,.02)
    im  = ax.contourf(xx,z/1000.,np.log10(ext.EXT),clevs,cmap=cm,extend='max')
    y = np.squeeze(z[:,71]/1000.)
    im2 = ax.plot(x,y,color="grey",lw=2)
    d = np.zeros(len(x))
    ax.fill_between(x,y,where=y>=d,color="beige")
    ax.set_ylim(0,20)
    ax.set_facecolor('black')

    for item in ([ax.title, ax.xaxis.label, ax.yaxis.label] +
             ax.get_xticklabels() + ax.get_yticklabels()):
        item.set_fontsize(20)

    cbar1 = plt.colorbar(im,ax=ax,ticks=[-2,-1,0],
                         format=mticker.FixedFormatter(['0.01', '0.1', '1']) )
    cbar1.ax.tick_params(labelsize=16)
    cbar1.set_label(label=f"{speciestitle} Aerosol Extinction [532 nm, km-1]",
                    size=16,rotation=270.,labelpad=25)
    plt.title(f"{modname} {speciestitle} Aerosol Extinction {f0dateout}+{fvdateout} at {trjtstr}", size=20)
    if species == None:
        species = "Total"
    ofname = f"{model}.{species}_extinction_curtain.{trjstr}.{f0dateout}+{fvdateout}00.png"
    print(ofname)
#    plt.show()
    plt.savefig(ofname)
    plt.close(fig)
    return

def sample(fname,trj=-100.,model='fp'):

    fname2 = fname.replace("_aer_","_asm_")
    fpdata = [fname,fname2]

#   From the filename get the forecast init and valid time
    f0dateout = fname[-33:-22]  # forecast init time
    fvdateout = fname[-21:-10]  # forecast valid time
    datedir   = fname[-33:-25]+"T"+fname[-24:-22]
    
#   Create the trajectory based on requested lat/lon
#   A hack for INSPYRE: look at "trj" value passed in
#    if trj < 0 then asking for a plot at constant longitude
#    if try > 0 then asking for a plot at constant latitude
    if trj > 0:
        lon  = np.arange(-120.,-79.25,0.25)
        lat  = np.full(len(lon),trj)
        alt  = np.full(len(lon),0.)
        tyme = np.full(len(lon),datetime.strptime(fvdateout,"%Y%m%d_%H"))
        trjstr = f"lat0{abs(int(trj)):0d}N"
    else:
        lat  = np.arange(30.,70.25,0.25)
        lon  = np.full(len(lat),trj)
        alt  = np.full(len(lon),0.)
        tyme = np.full(len(lon),datetime.strptime(fvdateout,"%Y%m%d_%H"))
        trjstr = f"lon{abs(int(trj)):0d}W"

    
    
# Get the model configuration
    do_optics = 1
    modname = model
    if(model == 'fp'):
        modname = "GEOS-FP"
        config = './geos543_pm25.yaml'
    elif(model == 'MERRA2'):
        config = './m2_pm25.yaml'
    else:
        print("STOP")
        sys.exit()            

#   Make an output directory
    dirname = f"samples/INSPYRE/sampled/{modname}/{datedir}"
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

    print(fpdata)
    kwargs = {'compat':'override'}
    traj = TRAJECTORY(tyme,lon,lat,fpdata,**kwargs)
# sample the dataset along the trajectory, and return an xarray dataset
    traj_ds = traj.sample()

#   Rename dimensions
    traj_ds = traj_ds.rename_dims(lev="level")
#   Rename variables
    traj_ds = traj_ds.rename_vars(lev="level", lon="longitude", lat="latitude")
    
#   Add global attributes
    dd   = f"{f0dateout}+{fvdateout}"
    fltrck = trjstr
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
                                    VersionID                   ="0",
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
    outFile = f"./{dirname}/INSPYRE-{modname}_Model_{dd}.traj_{trjstr}.nc"
    traj_ds.to_netcdf(outFile)
    print(f"Wrote: {outFile}")
    traj_ds.close()
    return

if __name__ == "__main__":
#-------------------------------------------------------
#   Parse the commandline
#-------------------------------------------------------

#   CHECK INPUT ARGUMENTS
    parser = OptionParser(usage="Usage: %prog [options] modelname date0",
                          version='xxx' )
    (options, args) = parser.parse_args()
 
#  GET OMI FILE FROM INPUT ARGUMENT LIST
    if len(args) == 1:
        fname      = args[0]
    else:
        parser.error("must have 1 argument: filename")
        
    sample(fname,trj=-100.)
    sample(fname,trj=-110.)
    sample(fname,trj=-120.)
    sample(fname,trj=-130.)
    sample(fname,trj=-80.)
    sample(fname,trj=40.)
    sample(fname,trj=50.)
    plot_curtain(fname,trj=-80.)
    plot_curtain(fname,trj=-100.)
    plot_curtain(fname,trj=-110.)
    plot_curtain(fname,trj=-120.)
    plot_curtain(fname,trj=-130.)
    plot_curtain(fname,trj=40.)
    plot_curtain(fname,trj=50.)
    plot_curtain(fname,trj=-80.,species="cc")
    plot_curtain(fname,trj=-100.,species="cc")
    plot_curtain(fname,trj=-110.,species="cc")
    plot_curtain(fname,trj=-120.,species="cc")
    plot_curtain(fname,trj=-130.,species="cc")
    plot_curtain(fname,trj=40.,species="cc")
    plot_curtain(fname,trj=50.,species="cc")

    sys.exit()
