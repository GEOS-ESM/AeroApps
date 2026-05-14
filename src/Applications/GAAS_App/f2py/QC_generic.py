#!/usr/bin/env python

import os
import sys
import argparse
import yaml
import numpy as np
import xarray as xr

# Default YAML file for SQC parameters
SQC_Param = """
reset_allqc: no
reset_passive: no
reset_sqc: yes

Backg_Check:
  tau_bgh: 2
  tau_bgx: 6

Buddy_Check:
  tau_buddy: 0.1
  niter_max: 5
  search_rad: 1
  nbuddy_max: 100
  nstar: 0
  ls_h: 150000.0
  ls_v: 1000.0     
"""

# Example sigO, sigF config for variables that require external error definitions
DEFAULT_SIGS = """  
 354: {sigO: 0.18, sigF: 0.45}
 388: {sigO: 0.18, sigF: 0.45}
 440: {sigO: 0.18, sigF: 0.45}
 470: {sigO: 0.18, sigF: 0.45}
 550: {sigO: 0.18, sigF: 0.45}
 660: {sigO: 0.18, sigF: 0.45}
 870: {sigO: 0.18, sigF: 0.45}
 1200: {sigO: 0.18, sigF: 0.45}
 1600: {sigO: 0.18, sigF: 0.45}
 2100: {sigO: 0.18, sigF: 0.45}
"""

class QC(object):

    def __init__(self, iodaFiles, var_name, z_coord_name=None, err_config=None, config=None, log=False, verbose=True, output_file=None):
        if config is None:
            config = SQC_Param

        self.var_name = var_name
        self.z_coord_name = z_coord_name
        self.log = log
        self.verbose = verbose
        
        # Load external error config or fallback to defaults
        if err_config is None:
            if self.verbose: print(f'[x] No err_config provided. Falling back to DEFAULT_SIGS.')
            self.err_config = yaml.safe_load(DEFAULT_SIGS)
        elif isinstance(err_config, str):
            self.err_config = yaml.safe_load(err_config)
        else:
            self.err_config = err_config
         
        # Automatic output filename logic
        if output_file is None:
            if isinstance(iodaFiles, str):
                base_name = os.path.basename(iodaFiles)
                if "preQC" in base_name:
                    self.output_file = base_name.replace("preQC", "postQC")
                else:
                    name, ext = os.path.splitext(base_name)
                    self.output_file = f"{name}_postQC{ext}"
            else:
                self.output_file = "ioda_out_sqc.nc" # Fallback if passing a Dataset object
        else:
            self.output_file = output_file

        # Load IODA file
        if isinstance(iodaFiles, xr.Dataset):
            self.ioda = iodaFiles
        else:
            self.ioda = xr.open_datatree(iodaFiles, decode_times=False)
            
        if self.verbose: print(f'[x] Loaded IODA file. Target variable: {self.var_name}') 
        
        self.nobs = self.ioda.dims['Location']
        # Handle cases where Channel might not exist (e.g., 1D variables like PM2.5)
        self.nchannels = self.ioda.dims.get('Channel', 1)
        
        # Extract Original shapes for later reshaping back
        self.orig_shape = self.ioda.PreQc[self.var_name].shape
        print(f"[x] Observations: {self.nobs}, Channels/Levels: {self.nchannels}")

        # Compute OMF (Observation Minus Forecast) natively
        try:
            hofx_key = 'HofX' if 'HofX' in self.ioda else 'hofx'
            obs = self.ioda.ObsValue[self.var_name].data
            bkg = self.ioda[hofx_key][self.var_name].data
            
            if self.log:  
                # Example: log scaling (adjust constants if needed for other variables)
                self.omf_data = np.log(obs + 0.01) - np.log(bkg + 0.01)
            else:
                self.omf_data = obs - bkg
                
            # Ensure OMF is 2D: (nobs, nchannels) for processing
            self.omf_data = self.omf_data.reshape(self.nobs, self.nchannels)
        except (AttributeError, KeyError) as e:
            raise ValueError(f"Couldn't create omf field for {self.var_name}. Ensure ObsValue and HofX exist. Error: {e}")
       
        self.sqc = yaml.safe_load(config)
        
        # Extract QC flags and ensure 2D structure
        self.qcexcl = self.ioda.PreQc[self.var_name].data.copy().reshape(self.nobs, self.nchannels)
        self.qchist = self.ioda.HistQc[self.var_name].data.copy().reshape(self.nobs, self.nchannels)

        # Handle Vertical/Channel Coordinates
        if self.z_coord_name and self.z_coord_name in self.ioda.MetaData:
            self.z_coord = self.ioda.MetaData[self.z_coord_name].data
        else:
            # Fallback to an array of zeros if no vertical coordinate applies
            self.z_coord = np.zeros(self.nchannels)
            
        self.lev = np.tile(np.array([self.z_coord]), (self.nobs, 1))

        # Compute Cartesian coordinates
        self.lons = self.ioda.MetaData['longitude'].data
        self.lats = self.ioda.MetaData['latitude'].data
        xobs, yobs, zobs = self.ll2xyz(self.lons, self.lats)
        self.xobs = np.array(xobs, dtype=np.float64)
        self.yobs = np.array(yobs, dtype=np.float64)
        self.zobs = np.array(zobs, dtype=np.float64)
   
        # Core Workflow
        self._reset_qc()
        self._getErrVar()
        self._backgCheck()
        self._reorder_arrays(message="Reordering by exclusion status")
        self.regions, self.nmaxregions, self.iregbeg, self.ireglen = self._getregions()
        self._buddyCheck()

        # Update QC values in the IODA file (reshaping back to original shape)
        self.ioda.PreQc[self.var_name].data = self.qcexcl.reshape(self.orig_shape)
        self.ioda.HistQc[self.var_name].data = self.qchist.reshape(self.orig_shape)

        # Save output
        self._save_ioda(self.output_file)

    def _reset_qc(self):
        if self.sqc['reset_allqc']:
            self.qcexcl[:] = 0
            self.qchist[:] = 0
            if self.verbose: print('[x] Reset all exclusion and history marks')
        elif self.sqc['reset_sqc']:
            self.qchist[:] = 0
            if self.verbose: print('[x] Reset all history marks')

    def _reorder_arrays(self, sort_indices=None, sort_by_key=None, descend=False, message=None):
        if message is None:
           message = 'Reordering arrays'
        if self.verbose: print(f'[x] {message}')

        # Determine the sort indices
        if sort_indices is None:
           if sort_by_key is None:
              mask = np.all(self.qcexcl == 0, axis=1)
              sort_indices = np.argsort(~mask) 
              n_valid = np.sum(mask)             
           else:
              indx = np.arange(len(sort_by_key), dtype=np.int32)
              sort_indices = indx[np.argsort(sort_by_key[indx], kind='mergesort')[::-1 if descend else 1]]
              n_valid = None
        else:
           n_valid = None

        # --- FIX: Preserve ROOT coordinates, variables, and attributes ---
        # Get root dataset
        if hasattr(self.ioda, 'dataset'):
            root_ds = self.ioda.dataset.copy()
        else:
            root_ds = self.ioda.to_dataset().copy() if hasattr(self.ioda, 'to_dataset') else self.ioda.copy()

        # Reindex variables in the root dataset
        for var in root_ds.variables:
            if 'Location' in root_ds[var].dims:
                root_ds[var] = root_ds[var].isel(Location=sort_indices)

        # Initialize the new DataTree using the reindexed root dataset
        new_tree = xr.DataTree(dataset=root_ds)
        new_tree.attrs = self.ioda.attrs

        # Get child groups (compatability across DataTree versions)
        children = self.ioda.children if hasattr(self.ioda, 'children') else {k: v for k, v in self.ioda.items() if k != '/'}

        # Process each child group
        for group_name, group in children.items():
            if hasattr(group, 'dataset'):
                group_ds = group.dataset.copy()
            else:
                group_ds = group.to_dataset().copy() if hasattr(group, 'to_dataset') else group.copy()
                
            # Reindex variables in the child group
            for var in group_ds.variables:
                if 'Location' in group_ds[var].dims:
                    group_ds[var] = group_ds[var].isel(Location=sort_indices)
            
            # Add child dataset to tree and preserve group attributes
            new_tree[group_name] = xr.DataTree(dataset=group_ds)
            new_tree[group_name].attrs = group.attrs

        self.ioda = new_tree
        # ----------------------------------------------------------------
      
        # Reorder numpy arrays
        self.qcexcl = self.qcexcl[sort_indices]
        self.qchist = self.qchist[sort_indices]
        self.omf_data = self.omf_data[sort_indices]
        
        if hasattr(self, 'VarO'): self.VarO = self.VarO[sort_indices]
        if hasattr(self, 'VarF'): self.VarF = self.VarF[sort_indices]

        for attr in ['xobs', 'yobs', 'zobs']:
            if hasattr(self, attr):
               setattr(self, attr, getattr(self, attr)[sort_indices])
    
        if n_valid is not None and self.verbose:
           print(f'[x] {message} complete: {n_valid} valid observations moved to the front')
        return n_valid

    def _getErrVar(self):
        print(f'[x] Getting Error Variances for {self.var_name}')
        self.VarO = np.full((self.nobs, self.nchannels), -999.)
        self.VarF = np.full((self.nobs, self.nchannels), -999.)
        
        # Ensure we have the config, as it is strictly required for the Background Error
        if self.err_config is None:
            raise ValueError(f"No err_config provided. Background error (sigF) requires YAML config.")

        print(f'[x] Setting Background Error (VarF) from YAML configuration.')
        
        # 1. Determine where ObsError (VarO) will come from
        get_varo_from_yaml = True
        if 'ObsError' in self.ioda and self.var_name in self.ioda.ObsError:
            print(f'[x] Found ObsError natively in IODA file. Using it for VarO.')
            err_data = self.ioda.ObsError[self.var_name].data
            self.VarO = (err_data.reshape(self.nobs, self.nchannels)) ** 2
            get_varo_from_yaml = False
        else:
            print(f'[x] ObsError not in IODA file. Will fallback to YAML configuration for VarO.')

        # 2. Extract Variances from YAML
        for ch, channel_val in enumerate(self.z_coord):
            found = False
            for key in self.err_config:
                if abs(channel_val - float(key)) < 0.01:
                    # Always set VarF from YAML
                    self.VarF[:, ch] = self.err_config[key]['sigF']**2
                    
                    # Set VarO from YAML only if it wasn't found in the IODA file
                    if get_varo_from_yaml:
                        self.VarO[:, ch] = self.err_config[key]['sigO']**2
                        
                    found = True
                    break
                    
            if not found:
                print(f'[x] Warning: No error variances found in config for coordinate {channel_val}')

    def _backgCheck(self):
        var = self.VarF + self.VarO
        if np.any(var <= 0):
           raise ValueError("Prescribed O-F Variance is zero or negative.")
    
        tau = np.sqrt(self.omf_data**2 / var)

        extreme_outliers = tau > self.sqc['Backg_Check']['tau_bgx']
        self.qcexcl[extreme_outliers] = 21
    
        suspect_outliers = tau > self.sqc['Backg_Check']['tau_bgh']
        self.qchist[suspect_outliers] = 17
    
        nexcl_per_channel = np.sum(extreme_outliers, axis=0)
        nsuspect_per_channel = np.sum(suspect_outliers, axis=0)
    
        for channel in range(nexcl_per_channel.size):
            z_val = self.z_coord[channel] if self.nchannels > 1 else 'Single-Level'
            print(f'[x] Exclusions after background_check for channel/level {z_val}: {nexcl_per_channel[channel].item()}')
            print(f'[x] Suspects after background_check for channel/level {z_val}: {nsuspect_per_channel[channel].item()}')

    @staticmethod
    def ll2xyz(lon, lat): 
        xobs = np.cos(np.radians(lat)) * np.cos(np.radians(lon))
        yobs = np.cos(np.radians(lat)) * np.sin(np.radians(lon))
        zobs = np.sin(np.radians(lat))
        return xobs, yobs, zobs

    def _getregions(self): 
        import Icosahedron_
        regions = np.zeros(self.nobs, dtype=np.int32)
        (regions, nmaxregions) = Icosahedron_.py_icosahedron_regions(self.xobs, self.yobs, self.zobs)
        
        original_regions = regions.copy()
        self._reorder_arrays(sort_by_key=regions, message="Sorting observations by region")
        
        indx = np.arange(len(original_regions), dtype=np.int32)
        sort_indices = indx[np.argsort(original_regions[indx], kind='mergesort')]
        regions = original_regions[sort_indices]
        if self.verbose: print(f'[x] Total max regions: {nmaxregions}, regions with obs: {len(np.unique(regions))}') 
        
        iregbeg = np.zeros(nmaxregions, dtype=np.int32) 
        ireglen = np.zeros(nmaxregions, dtype=np.int32) 
    
        np.add.at(ireglen, regions - 1, 1)
        iregbeg[1:] = np.cumsum(ireglen[:-1]) 
        
        return regions, nmaxregions, iregbeg, ireglen

    def _buddyCheck(self):
        import Icosahedron_
        print(f'[x] Starting the buddy check')
        
        bc_conf = self.sqc['Buddy_Check']
        search_rad = bc_conf['search_rad'] 
        ls_h = bc_conf['ls_h']  
        ls_v = bc_conf.get('ls_v', 1000.0) 
        single_level = bc_conf.get('single_level', True)
        tau_buddy = bc_conf.get('tau_buddy', 0.1) 
        nbuddy_max = bc_conf.get('nbuddy_max', 100)
        seplim = bc_conf.get('seplim', 26.5) 

        reaccept = np.zeros((self.nobs, self.nchannels), dtype=bool)
    
        for channel in range(self.nchannels):
            z_val = self.z_coord[channel] if self.nchannels > 1 else 'Single-Level'
            if self.verbose: print(f'\n[x] Processing channel/level {z_val}')
            
            suspect = (self.qchist[:, channel] > 0) & (self.qcexcl[:, channel] == 0)
            if self.verbose:
                print(f"Suspect indices for channel/level {z_val}:", np.where(suspect)[0])

            for iteration in range(1, bc_conf['niter_max'] + 1): 
                if self.verbose: print(f'[x] Iteration {iteration} for channel/level {z_val}')
                
                if not np.any(suspect):
                   if self.verbose: print(f'[x] No suspect observations for channel/level {z_val} at iteration {iteration}')
                   break     
                    
                suspect_indices = np.where(suspect)[0]
                n_susp = len(suspect_indices)
                if self.verbose: print(f'[x] Number of suspects for channel/level {z_val} at iteration {iteration} is {n_susp}')

                obs_regions = np.zeros(self.nobs, dtype=int)
                for ireg in range(1, self.nmaxregions + 1):
                    if self.ireglen[ireg-1] > 0:
                       ibeg = self.iregbeg[ireg-1]             
                       iend = ibeg + self.ireglen[ireg-1] - 1
                       obs_regions[ibeg:iend+1] = ireg 
                       
                ki_susp = suspect_indices + 1 
                kr_susp = obs_regions[suspect_indices]  

                iregbeg_fortran = self.iregbeg + 1 
                
                reaccept[:, channel] = Icosahedron_.py_find_buddies(ki_susp, kr_susp, self.xobs, 
                        self.yobs, self.zobs, self.lats, self.lons, self.lev[:,channel], self.omf_data[:,channel],
                        self.VarF[:,channel], self.VarO[:,channel], self.qcexcl[:,channel], ls_h, ls_v, search_rad, tau_buddy,
                        single_level, nbuddy_max, iregbeg_fortran, self.ireglen, seplim)
             
                reaccepted_obs = np.where(reaccept[:,channel])[0]
                if len(reaccepted_obs) > 0:
                    if self.verbose: print(f"[x] Reaccepted {len(reaccepted_obs)} observations for channel/level {z_val}")
                    suspect[reaccept[:,channel]] = False
                else:
                    if self.verbose: print(f"[x] No observations reaccepted for channel/level {z_val} at iteration {iteration}")
                    break 

            # Update QC flags: un-reaccepted suspects get QCX = 17 
            self.qcexcl[suspect, channel] = 17
            
        print(f'[x] Buddy check completed')
        
    def _save_ioda(self, output_filename):
        print(f'[x] Updating Fill Values to JEDI standards...')
        
        # Standard JEDI / NetCDF Fill Values
        JEDI_FILL_FLOAT = 9.969209968386869e+36
        JEDI_FILL_INT32 = -2147483647
        JEDI_FILL_INT64 = -9223372036854775806

        # Gather all nodes (root + children) in the Datatree
        nodes = [self.ioda]
        if hasattr(self.ioda, 'children'):
            nodes.extend(self.ioda.children.values())
        else:
            nodes.extend([v for k, v in self.ioda.items() if k != '/'])

        # Apply JEDI compliant fill values to all variables
        for node in nodes:
            ds = node.dataset if hasattr(node, 'dataset') else node
            for var_name, da in ds.variables.items():
                
                # 1. Remove _FillValue from attributes to prevent xarray conflicts
                if '_FillValue' in da.attrs:
                    del da.attrs['_FillValue']
    
                # 2. Determine the TARGET dtype (what it will actually be saved as on disk), 
                # because xarray might have promoted ints to floats in memory due to NaNs.
                target_dtype = np.dtype(da.encoding.get('dtype', da.dtype))


                # 3. Set the correct _FillValue in the encoding dictionary based on data type
                if target_dtype.kind == 'f':  # float32, float64
                    da.encoding['_FillValue'] = target_dtype.type(JEDI_FILL_FLOAT)
                elif target_dtype == np.int64:
                    da.encoding['_FillValue'] = np.int64(JEDI_FILL_INT64)
                elif target_dtype.kind in ['i', 'u']: # int32, int16, etc.
                    da.encoding['_FillValue'] = np.int32(JEDI_FILL_INT32)

        print(f'[x] Saving updated IODA file to {output_filename}...')
        try:
            self.ioda.to_netcdf(output_filename)
            print(f'[x] Successfully saved: {output_filename}')
        except Exception as e:
            print(f'[x] Error saving IODA file: {e}')

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Perform Background and Buddy Check on IODA observations.")
    parser.add_argument("input_file", type=str, help="Path to the input IODA NetCDF file.")
    parser.add_argument("--var", type=str, required=True, help="Target variable name (e.g., pm2_5, aerosolOpticalDepth).")
    parser.add_argument("--z-coord", type=str, default=None, help="Vertical/Channel coordinate name in MetaData (e.g., obs_wavelength).")
    parser.add_argument("--err-config", type=str, default=None, help="Path to YAML file defining sigO/sigF if ObsError isn't in the IODA file.")
    parser.add_argument("-o", "--output", type=str, default=None, help="Output file name (optional, will auto-generate if omitted).")
    parser.add_argument("--log", action="store_true", help="Enable log transformation for OMF calculation.")
    
    args = parser.parse_args()
    
    # Load external Error Config if a file path is provided
    err_config_dict = None
    if args.err_config:
        with open(args.err_config, 'r') as f:
            err_config_dict = yaml.safe_load(f)
            
    print(f"Running QC on {args.input_file} for variable {args.var}...")
    
    qc = QC(
        iodaFiles=args.input_file, 
        var_name=args.var,
        z_coord_name=args.z_coord,
        err_config=err_config_dict,
        log=args.log, 
        output_file=args.output
    )
