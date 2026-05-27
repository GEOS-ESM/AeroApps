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
AC update (4/17/26): generates a mask for each
timestep to determine where grid boxes are available.
This mask is then used on all variables to speed up 
the processing.
'''


import numpy as np
import os
import sys
import xarray as xr
import re
from glob import glob
import dask
from dask.diagnostics import ProgressBar
import time

def file_list(sat, yy, mm, atype="blend", obsdir='./'):
    """
    Get list of observation files for a given satellite, year, and month
    """
    if not os.path.exists(obsdir):
        raise FileNotFoundError(f"Observation directory {obsdir} does not exist!")
    obsfiles = sorted(os.listdir(obsdir))
    res1 = [x for x in obsfiles if any(y in x for y in [atype])]
    res2 = [x for x in res1 if any(y in x for y in ['z.nc4'])]
    obsfiles = res2
    if not obsfiles:
        raise ValueError(f"No observation files found in directory {obsdir}")
    # Prepend directory name
    obsfiles = [f"{obsdir}{file}" for file in obsfiles]
    return obsfiles

def modelfile_list(yy, mm, expid, expdir='./'):
    """
    Get model files matching the pattern:
    {expid}.aer_inst_1hr_glo_L1152x721_slv.%y4-%m2-%d2T%h200Z.nc4
    Only files with 3-hourly timestamps (00, 03, 06, 09, 12, 15, 18, 21)
    """
    if not os.path.exists(expdir):
        raise FileNotFoundError(f"Model directory {expdir} does not exist!")
    
    # Build pattern for matching files
    pattern = f"{expid}.inst2d_hwl_x.{yy:04d}{mm:02d}*00z.nc4"
    
    # Use glob to find matching files
    search_pattern = os.path.join(expdir, pattern)
    all_files = sorted(glob(search_pattern))
    
    if not all_files:
        raise ValueError(f"No model files found matching pattern: {search_pattern}")
    
    # Filter for 3-hourly files only (00, 03, 06, 09, 12, 15, 18, 21)
    hourlist = ['0000z', '0300z', '0600z', '0900z', '1200z', '1500z', '1800z', '2100z']
    expfiles = [f for f in all_files if any(hour in f for hour in hourlist)]
    
    if not expfiles:
        raise ValueError(f"No 3-hourly model files found in directory {expdir}")
    
    print(f"Found {len(expfiles)} model files for {yy:04d}-{mm:02d}")
    
    return expfiles

def monthlymodel_optimized(yy, mm, sat, expid, varlist=None,
                           expdir="./", obsdir="./", outdir="./", 
                           weighted=False, verbose=True):
    """
    
    Parameters:
    -----------
    yy : int
        Year
    mm : int
        Month
    sat : str
        Satellite name (e.g., 'MOD04', 'MYD04')
    expid : str
        Experiment ID
    varlist : list, optional
        List of variable names to process. If None, processes all variables.
    expdir : str
        Directory containing model files
    obsdir : str
        Directory containing observation files
    outdir : str
        Directory for output files
    weighted : bool
        Not currently used
    verbose : bool
        Print detailed diagnostics
    
    Returns:
    --------
    outfile : str
        Path to output file
    """
    
    start_time = time.time()
    
    # Get files
    if not os.path.exists(obsdir):
        raise FileNotFoundError(f"Observation directory {obsdir} does not exist!")
    obsfiles = file_list(sat, yy, mm, obsdir=obsdir)
    if not obsfiles:
        raise ValueError(f"No observation files found in directory {obsdir}")
    
    if not os.path.exists(expdir):
        raise FileNotFoundError(f"Model directory {expdir} does not exist!")
    modfiles = modelfile_list(yy, mm, expid, expdir=expdir)
    
    if verbose:
        print(f"\nFound {len(obsfiles)} observation files for {sat} {yy:04d}-{mm:02d}")
        print(f"Found {len(modfiles)} model files")
    
    # Step 1: Load observation data
    if verbose:
        print(f"\n[1/5] Loading observation data...")
    t1 = time.time()
    obs_ds = xr.open_mfdataset(obsfiles, engine='netcdf4')
    tau_ = obs_ds["tau_"]
    
    # Squeeze out the 'lev' dimension if it exists and has size 1
    if 'lev' in tau_.dims and len(tau_.lev) == 1:
        tau_ = tau_.squeeze('lev', drop=True)
    
    if verbose:
        print(f"  Observation tau_ shape: {tau_.shape}")
        print(f"  Observation tau_ dims: {tau_.dims}")
        print(f"  Number of timesteps in obs: {len(tau_.time)}")
    
    # Step 2: Load ALL model data (not just first file) to get proper time coordinates
    if verbose:
        print(f"\n[2/5] Loading model data to get time coordinates...")
    t2 = time.time()
    
    # Load first file to get grid structure and variable list
    mod_sample = xr.open_dataset(modfiles[0], engine='netcdf4')
    
    # Determine variables to load
    if varlist is None:
        varlist = [v for v in mod_sample.data_vars]
    
    # Now load all model files with only the needed variables
    mod_ds = xr.open_mfdataset(modfiles, engine='netcdf4', 
                               data_vars=varlist,
                               parallel=True)
    
    if verbose:
        print(f"  Model data shape: time={len(mod_ds.time)}, lat={len(mod_ds.lat)}, lon={len(mod_ds.lon)}")
        print(f"  Loaded {len(varlist)} variables in {time.time()-t2:.1f}s")
    
    # Step 3: Interpolate observation mask to model grid
    if verbose:
        print(f"\n[3/5] Interpolating mask to model grid...")
    t3 = time.time()
    

    taui_ = tau_.interp_like(mod_ds, method="nearest", assume_sorted=True)
    
    if verbose:
        print(f"  Interpolated mask shape: {taui_.shape}")
        print(f"  Interpolated mask dims: {taui_.dims}")
    
    if verbose:
        print(f"  Loading 3D mask into memory...")
    mask = taui_.isnull().compute()
    
    if verbose:
        print(f"\n  === MASK VERIFICATION ===")
        print(f"  Mask shape: {mask.shape}")
        print(f"  Mask dimensions: {mask.dims}")
        print(f"  Mask memory size: {mask.nbytes / 1024**2:.1f} MB")
        
        # Check if mask varies by timestep
        if 'time' in mask.dims and len(mask.time) > 1:
            # Compare first and last timesteps
            mask_t0 = mask.isel(time=0)
            mask_t1 = mask.isel(time=1)
            mask_tlast = mask.isel(time=-1)
            
            diff_01 = (mask_t0 != mask_t1).sum().values
            diff_0last = (mask_t0 != mask_tlast).sum().values
            
            print(f"  ✓ Mask HAS time dimension with {len(mask.time)} timesteps")
            print(f"  Difference between timestep 0 and 1: {diff_01} pixels")
            print(f"  Difference between timestep 0 and last: {diff_0last} pixels")
            
            # Show percentage masked for a few timesteps
            print(f"\n  Percentage masked by timestep:")
            timestep_indices = [0, len(mask.time)//4, len(mask.time)//2, 3*len(mask.time)//4, -1]
            for t_idx in timestep_indices:
                pct = mask.isel(time=t_idx).sum().values / (mask.shape[-2] * mask.shape[-1]) * 100
                t_label = f"{t_idx}" if t_idx >= 0 else "last"
                print(f"    Timestep {t_label}: {pct:.1f}% masked")
        else:
            print(f"  ✗ WARNING: Mask does NOT have time dimension or only has 1 timestep!")
            print(f"     This means the same mask will be applied to all timesteps!")
        
        print(f"  Overall % masked: {mask.sum().values / mask.size * 100:.1f}%")
        print(f"  Mask computed in {time.time()-t3:.1f}s")
    
    # Close obs dataset, we don't need it anymore
    obs_ds.close()
    
    # Verify time alignment
    if verbose:
        if len(mask.time) != len(mod_ds.time):
            print(f"\n  WARNING: Time mismatch!")
            print(f"    Mask timesteps: {len(mask.time)}")
            print(f"    Model timesteps: {len(mod_ds.time)}")
        else:
            print(f"  ✓ Mask and model time dimensions aligned ({len(mask.time)} timesteps)")
    
    # Step 4: Apply mask and compute means for all variables at once
    if verbose:
        print(f"\n[4/5] Computing masked means for {len(varlist)} variables...")
    t4 = time.time()
    
    # Build a dictionary of masked means (computed together for efficiency)
    masked_means = {}
    for varn in varlist:
        if verbose:
            print(f"  Setting up {varn}...")
        # Apply mask and compute mean
        # Since mask has (time, lat, lon) and mod_ds[varn] has (time, lat, lon),
        # this is element-wise per-timestep masking
        masked_means[varn] = mod_ds[varn].where(~mask).mean(dim="time")
    
    # Compute all at once (this is much faster than one at a time)
    if verbose:
        print(f"\n  Computing all variables together...")
    with ProgressBar():
        computed_means = xr.Dataset(masked_means).compute()
    
    if verbose:
        print(f"  Computed all means in {time.time()-t4:.1f}s")
    
    # Step 5: Create output file
    if verbose:
        print(f"\n[5/5] Creating output file...")
    t5 = time.time()
    
    # Create output dataset from template
    out_ds = mod_sample.copy()
    
    # Update with computed means
    for varn in varlist:
        out_ds[varn][0, :, :] = computed_means[varn]
    
    # Drop variables not in varlist
    vars_to_drop = [v for v in out_ds.data_vars if v not in varlist]
    if vars_to_drop:
        out_ds = out_ds.drop_vars(vars_to_drop)
        if verbose:
            print(f"  Dropped {len(vars_to_drop)} unused variables")
    
    # Save output
    basename = os.path.basename(modfiles[0]).split('.inst2d_hwl_x')[0]
    outfile = os.path.join(outdir, f"{basename}.monthly.{sat}.{yy:04d}{mm:02d}.nc4")
    os.makedirs(outdir, exist_ok=True)
    
    if verbose:
        print(f"  Writing to: {outfile}")
    out_ds.to_netcdf(path=outfile, format="NETCDF4")
    if verbose:
        print(f"  Written in {time.time()-t5:.1f}s")
    
    # Cleanup
    mod_ds.close()
    mod_sample.close()
    out_ds.close()
    
    total_time = time.time() - start_time
    if verbose:
        print(f"\n✓ Total processing time: {total_time:.1f}s ({total_time/60:.1f} minutes)")
    
    return outfile


def monthlymodel(yy, mm, sat, expid, varlist=None,
                 expdir="./", obsdir="./", outdir="./", weighted=False):
    """
    Wrapper function to maintain backward compatibility
    Calls the optimized version with default settings
    """
    return monthlymodel_optimized(yy, mm, sat, expid, varlist=varlist,
                                  expdir=expdir, obsdir=obsdir, outdir=outdir,
                                  weighted=weighted, verbose=True)


# main script
if __name__ == "__main__":
    varlist = ["TOTEXTTAU", "DUEXTTAU", "SSEXTTAU", "OCEXTTAU",
               "BCEXTTAU", "SUEXTTAU", "NIEXTTAU", "BREXTTAU",
               "TOTANGSTR"]
    yy = 2019
    expids = ["c180R_v11.8.0_develop_newvolc"]
    
    for expid in expids:

        expdir = f"/home/pcolarco/geos_aerosols/pcolarco/{expid}/"
        obsdir = f"/home/pcolarco/geos_aerosols/pcolarco/"
        for mm in range(1, 13):
            for sat in ["MOD04", "MYD04"]:
                print(f"\n{'='*70}")
                print(f"Processing: {expid} | {yy:04d}-{mm:02d} | {sat}")
                print(f"{'='*70}")
                
                try:
                    outfile = monthlymodel_optimized(
                        yy, mm, sat, expid,
                        expdir=f"{expdir}holding/inst2d_hwl_x/{yy:04d}{mm:02d}/",
                        obsdir=f"{obsdir}/{sat}/Level3/Y{yy:04d}/M{mm:02d}/",
                        outdir=f"./{sat}sampled/{expid}/",
                        varlist=varlist,
                        verbose=True
                    )
                    print(f"\nSUCCESS: {outfile}\n")
                except Exception as e:
                    print(f"\nERROR processing {expid} {yy:04d}-{mm:02d} {sat}: {e}")
                    import traceback
                    traceback.print_exc()
                    continue
                
    print(f"\n{'='*70}")
    print("All processing complete!")
    print(f"{'='*70}")
