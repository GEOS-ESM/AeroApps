#!/usr/bin/env python3

import numpy as np
import os
import sys
import xarray as xr

#script to compare model to monthly mean MODIS NNR Retrievals
# use blendlanddeep_MYD.py to generate blend files first

def file_list(sat,yy,mm,atype="blend",obsdir='./'):
    if not os.path.exists(obsdir):
        raise FileNotFoundError(f"Observation directory {obsdir} does not exist!")
    obsfiles = sorted(os.listdir(obsdir))
    res1 = [x for x in obsfiles if any (y in x for y in [atype])]
    res2 = [x for x in res1 if any (y in x for y in ['z.nc4'])]
    obsfiles = res2
    if not obsfiles:
        raise ValueError(f"No observation files found in directory {obsdir}")
    # Prepend directory name and open as a dataset
    i = 0
    for file in obsfiles:
        obsfiles[i] = f"{obsdir}{file}"
        i += 1
    return obsfiles

def monthlyobs(yy, mm, sat, obsdir="./", outdir="./",weighted=False):

#   Get the observations files and open as a dataset
    if not os.path.exists(obsdir):
        raise FileNotFoundError(f"Observation directory {obsdir} does not exist!")
    obsfiles = file_list(sat,yy,mm,obsdir=obsdir)
    if not obsfiles:
        raise ValueError(f"No observation files found in directory {obsdir}")
    # Extract first and last available days in the dataset
    dds = int(obsfiles[0][-12:-10])     # First day
    dde = int(obsfiles[-1][-12:-10])  # Last day
    num_days = dde - dds + 1
    print(f"First and Last day of the month with data: {dds}, {dde}")

#   Open dataset
    obs_ds = xr.open_mfdataset(obsfiles,engine='netcdf4')
#   open a single file as an output template
    out_ds = xr.open_dataset(obsfiles[0],engine='netcdf4')

#   Sum the counts
    out_ds["count_tau_"].values[0:,:,:] = obs_ds["count_tau_"].sum(dim="time").values.copy()
    if not weighted:
#       Simple average over time dimension
        out_ds["tau_"].values[0,:,:,:] = obs_ds['tau_'].mean(dim="time").values.copy()
    else:
#       Weighted average over time dimension
        tau_ = np.squeeze(obs_ds["tau_"].values.copy())
        cnt_ = np.squeeze(obs_ds["count_tau_"].values.copy())
        num_ = np.nansum(tau_*cnt_,axis=0)
        den_ = np.squeeze(out_ds["count_tau_"].values)
        tao_ = np.full_like(num_,np.nan)
        mask = den_ > 0
        if np.sum(mask) > 0:
            tao_[mask] = num_[mask]/den_[mask]
        out_ds["tau_"].values[0,0,:,:] = tao_
    
    basename = os.path.basename(obsfiles[0])[0:-18]
    outfile = outdir+basename+f"monthly.{yy:04d}{mm:02d}.nc4"
    if weighted:
        outfile = outdir+basename+f"monthly.weighted.{yy:04d}{mm:02d}.nc4"
    out_ds.to_netcdf(path=outfile,format="NETCDF4")
    return outfile


# main script
if __name__ == "__main__":
    yy = 2019
    sat = "MOD04"
    for mm in range(1,13):
        outfile=monthlyobs(yy,mm,sat,weighted=True,
                           obsdir=f"./{sat}/Level3/Y{yy:04d}/M{mm:02d}/",
                           outdir=f"./{sat}/Level3/Y{yy:04d}/M{mm:02d}/")
        print(outfile)
        outfile=monthlyobs(yy,mm,sat,
                           obsdir=f"./{sat}/Level3/Y{yy:04d}/M{mm:02d}/",
                           outdir=f"./{sat}/Level3/Y{yy:04d}/M{mm:02d}/")
        print(outfile)

    yy = 2019
    sat = "MYD04"
    for mm in range(1,13):
        outfile=monthlyobs(yy,mm,sat,weighted=True,
                           obsdir=f"./{sat}/Level3/Y{yy:04d}/M{mm:02d}/",
                           outdir=f"./{sat}/Level3/Y{yy:04d}/M{mm:02d}/")
        outfile=monthlyobs(yy,mm,sat,
                           obsdir=f"./{sat}/Level3/Y{yy:04d}/M{mm:02d}/",
                           outdir=f"./{sat}/Level3/Y{yy:04d}/M{mm:02d}/")
        print(outfile)
        print(outfile)











