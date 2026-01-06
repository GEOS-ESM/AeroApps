#!/usr/bin/env python3

'''
Sample model fields at MOD04 and MYD04 viewing
Assumes all model output available at 3 hourly
cadence same as MODIS files. Regrids MODIS NNR
to model grid assuming nearest neighbor. Model
is masked where MODIS has no data and a monthly
average is computed. This is slow if run for,
e.g., all the variables in the inst2d_hwl_x
collection. There is an optional argument
"varlist" that allows you to pass a list of
variables to sample.
ToDo:
- can this be made more efficient?
- output file is simply pasted over a copied
  dataset of one of the model files
- output file contains all variables, not
  just the ones sampled; should reduce
'''

import numpy as np
import os
import sys
import xarray as xr

def file_list(sat,yy,mm,atype="blend",obsdir='./'):
#   Return a list of file names (per month) in the obsdir directory
#   Only capture hour dated files (e.g., *_0000z.nc4)
    rc = 0
    obsfiles=None
    if not os.path.exists(obsdir):
        print(f"Observation directory {obsdir} does not exist!")
        rc = 1
        return obsfiles, rc
    obsfiles = sorted(os.listdir(obsdir))
    res1 = [x for x in obsfiles if any (y in x for y in [atype])]
    res2 = [x for x in res1 if any (y in x for y in ['z.nc4'])]
    obsfiles = res2
    if not obsfiles:
        print(f"No observation files found in directory {obsdir}")
        rc = 1
        return obsfiles, rc
    # Prepend directory name and open as a dataset
    i = 0
    for file in obsfiles:
        obsfiles[i] = f"{obsdir}{file}"
        i += 1
    return obsfiles, rc

def modelfile_list(yy,mm,expdir='./'):
    rc = 0
    expfiles=None
    if not os.path.exists(expdir):
        print(f"Observation directory {expdir} does not exist!")
        rc = 1
        return expfile, rc
    expfiles = sorted(os.listdir(expdir))
    # Only get every 3 hours
    hourlist = ['0000z','0300z','0600z','0900z','1200z','1500z','1800z','2100z']
    res2 = [x for x in expfiles if any (y in x for y in hourlist)]
    expfiles = res2
    if not expfiles:
        raise ValueError(f"No model files found in directory {expdir}")
    # Prepend directory name and open as a dataset
    i = 0
    for file in expfiles:
        expfiles[i] = f"{expdir}{file}"
        i += 1
    return expfiles, rc

def monthlymodel(yy, mm, sat, expid, varlist=None, atype="blend",
                 expdir="./", obsdir="./", weighted=False):
    outfile = None
    
#   Get the observations files and open as a dataset
    if not os.path.exists(obsdir):
        print(f"Observation directory {obsdir} does not exist!")
        return outfile
    obsfiles, rc = file_list(sat,yy,mm,atype=atype,obsdir=obsdir)
    if rc != 0:
        print(f"No observation files found in directory {obsdir}")
        return outfile
    obs_ds = xr.open_mfdataset(obsfiles,engine='netcdf4')
    tau_ = obs_ds["tau_"].squeeze()
    
#   Get the model files
    if not os.path.exists(expdir):
        print(f"Model directory {expdir} does not exist!")
        return outfile
    modfiles, rc = modelfile_list(yy,mm,expdir=expdir)
    if rc != 0:
        print(f"No model files found in directory {expdir}")
        return outfile

#   Get the model files to sample
    mod_ds = xr.open_mfdataset(modfiles,engine='netcdf4')
    taui_ = tau_.interp_like(mod_ds,method="nearest",assume_sorted=True)
    tauin_ = taui_.isnull(keep_attrs=True)
    out_ds = xr.open_dataset(modfiles[0],engine='netcdf4')

#   Do masking and time average - this is the part that takes long time
    if varlist is None:
        for varn in mod_ds.data_vars:
            print(varn)
            mod_ds[varn] = mod_ds[varn].where(~tauin_)
            out_ds[varn][0,:,:] = mod_ds[varn].mean(dim="time")
    else:
        for varn in varlist:
            mod_ds[varn] = mod_ds[varn].where(~tauin_)
            out_ds[varn][0,:,:] = mod_ds[varn].mean(dim="time")
        for varn in mod_ds.data_vars:
            if varn not in varlist:
                out_ds = out_ds.drop_vars([varn])

#   Store in output
    basename = os.path.basename(modfiles[0])[0:-18]
    outfile = expdir+basename+f"monthly.{sat}.{yy:04d}{mm:02d}.nc4"
    out_ds.to_netcdf(path=outfile,format="NETCDF4")

    return outfile


# main script
if __name__ == "__main__":
    varlist = ["TOTEXTTAU550","DUEXTTAU550","SSEXTTAU550","OCEXTTAU550",
               "BCEXTTAU550","SUEXTTAU550","NIEXTTAU550","TOTSCATAU550",
               "TOTANGSTR","PM","PM25"]
    yy = 2019
    expids  = ["c180R_v11.8.0_develop","c180R_v11.8.0_newbrcoptics",
               "c180R_v11.8.0_newdust","c180R_v11.8.0_newdust_nogvf"]
    for expid in expids:
        expdir = f"/home/pcolarco/geos_aerosols/pcolarco/{expid}/"
        obsdir = f"/home/pcolarco/geos_aerosols/pcolarco/"
        for mm in range(1,13):
            for sat in ["MOD04","MYD04"]:
                outfile=monthlymodel(yy,mm,sat,expid,
                        expdir=f"{expdir}/holding/inst2d_hwl_x/{yy:04d}{mm:02d}/",
                        obsdir=f"{obsdir}/{sat}/Level3/Y{yy:04d}/M{mm:02d}/",
                        varlist=varlist)
                print(outfile)

    yy = 2020
    expids  = ["c180R_v11.8.0_develop"]
    for expid in expids:
        expdir = f"/home/pcolarco/geos_aerosols/pcolarco/{expid}/"
        for mm in range(1,13):
            for sat in ["MOD04","MYD04"]:
                outfile=monthlymodel(yy,mm,sat,expid,
                        expdir=f"{expdir}/holding/inst2d_hwl_x/{yy:04d}{mm:02d}/",
                        obsdir=f"{obsdir}/{sat}/Level3/Y{yy:04d}/M{mm:02d}/",
                        varlist=varlist)
                print(outfile)
