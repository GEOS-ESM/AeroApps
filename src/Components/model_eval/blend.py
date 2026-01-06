#!/usr/bin/env python3
'''
Colarco, January 2026
Code to blend MODIS ocean/deep/land into a single file
using algorithm of Allie Collow:
- where deep blue and land exist use a count-weighted average
- where only one of deep blue or land exist use that value
- superimpose onto ocean grid so that where land/deep points
  exist they are used instead of ocean
Write to similar format file as inputs (but only the NNR field)
'''

import numpy as np
import os
import sys
import xarray as xr

#script to compare model to monthly mean MODIS NNR Retrievals
# use blendlanddeep_MYD.py to generate blend files first


def file_list(sat,yy,mm,atype="ocean",obsdir='/css/gmao/dp/gds/AeroObs/'):
    rc = 0
    obs_files=None
    obs_dir = f"{obsdir}nnr_003_{sat}_061/Level3/Y{yy:04d}/M{mm:02d}/"
    if not os.path.exists(obs_dir):
        print(f"Observation directory {obs_dir} does not exist!")
        rc = 1
        return obs_files, rc
    obs_files = sorted(os.listdir(obs_dir))
    res1 = [x for x in obs_files if any (y in x for y in [atype])]
    res2 = [x for x in res1 if any (y in x for y in ['z.nc4'])]
    obs_files = res2
    if not obs_files:
        print(f"No observation files found in directory {obs_dir}")
        rc = 1
        return obs_files, rc
    # Prepend directory name and open as a dataset
    i = 0
    for file in obs_files:
        obs_files[i] = f"{obs_dir}{file}"
        i += 1
    return obs_files, rc


def blend(ocean,land,deep,outdir="./"):
    if not os.path.exists(outdir):
        os.makedirs(outdir, exist_ok=True)  # Recursively create directories if they don't exist

#   Open the datasets -- need handle if not exist
    ocean_ds = xr.open_dataset(ocean)
    oceantau = ocean_ds["tau_"].values.copy()
    oceannum = ocean_ds["count_tau_"].values.copy()

    try:
        land_ds  = xr.open_dataset(land)
        landtau  = land_ds["tau_"].values.copy()
        landnum  = land_ds["count_tau_"].values.copy()
    except:
        print(f"Missing land file:{land}")
        landtau = np.full_like(oceantau,0.)
        landnum = np.full_like(oceannum,0)
        
    try:
        deep_ds  = xr.open_dataset(deep)
        deeptau  = deep_ds["tau_"].values.copy()
        deepnum  = deep_ds["count_tau_"].values.copy()
    except:
        print(f"Missing deep file:{deep}")
        deeptau = np.full_like(oceantau,0.)
        deepnum = np.full_like(oceannum,0)

#   Make an output dataset and remove extraneous fields
    out_ds = ocean_ds.drop_vars(["tau","count_tau","tau_fine","cloud"])

#   Blend algorithm
    total_land_deep_obs = landnum+deepnum
    land_deep_blend = np.full_like(deeptau, np.nan)

    # Screen for NaN
    deeptau[deepnum == 0] = 0
    landtau[landnum == 0] = 0
    
    # Only calculate where we have observations
    mask = total_land_deep_obs > 0
    if np.sum(mask) > 0:
        land_deep_blend[mask] = ((deeptau[mask] * deepnum[mask]) + 
                                 (landtau[mask] * landnum[mask])) / total_land_deep_obs[mask]

#    # Return NaN to still zero -- not sure I need this
#    land_deep_blend[land_deep_blend == 0] = np.nan
#    And she had some logic about using ocean only where blend was nan
    
        # And drop onto the ocean grid
        oceantau[mask] = land_deep_blend[mask]
        oceannum[mask] = total_land_deep_obs[mask]

#   Put on output dataset
    out_ds['tau_'].values = oceantau
    out_ds['count_tau_'].values = oceannum

#   Get output filename and write
    basename = os.path.basename(ocean)
    outfile  = outdir+basename.replace("ocean","blend")
    out_ds.to_netcdf(path=outfile,format="NETCDF4")
        
    return outfile    

if __name__ == "__main__":
    for sat in ["MYD04"]:
        for yy in np.arange(2000,2026):
            for mm in np.arange(1,13):
                files, rc = file_list(sat,yy,mm)
                if rc == 0:
                    for ocean in files:
                        land = ocean.replace("ocean","land")
                        deep = ocean.replace("ocean","deep")
                        outfile = blend(ocean,land,deep,outdir=f"./{sat}/Level3/Y{yy:04d}/M{mm:02d}/")
                        print(f"Wrote file {outfile}")
