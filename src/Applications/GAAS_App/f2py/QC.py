#!/usr/bin/env python
# coding: utf-8



import xarray as xr
import numpy as np
from numba import jit, prange
import yaml
import sys

sys.path.append('/home/vbuchard/workspace/JEDI/OBS_IODA/GMAOpyobs/install/lib/Python/')




# Default YAML file for SQC parameters
SQC_Param = """
reset_allqc: yes
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

# sigO, sigF for LAOD:
LAOD_sigs = """  
 354:   
   sigO: 0.18  
   sigF: 0.45
 388:   
   sigO: 0.18  
   sigF: 0.45
 470:   
   sigO: 0.18  
   sigF: 0.45
 550:   
   sigO: 0.18  
   sigF: 0.45
 660:   
   sigO: 0.18  
   sigF: 0.45
 870:   
   sigO: 0.18  
   sigF: 0.45
 1200:   
   sigO: 0.18  
   sigF: 0.45
 1600:   
   sigO: 0.18  
   sigF: 0.45
 2100:   
   sigO: 0.18  
   sigF: 0.45
"""


# In[3]:


class QC(object):
    def __init__(self, iodaFiles, config=None, log=True, verbose=True):
        if config is None:
            config = SQC_Param

        self.log = log
        self.verbose = verbose
    
         
        if isinstance(iodaFiles, xr.Dataset):
            self.ioda = iodaFiles
        else:
            self.ioda = xr.open_datatree(iodaFiles)

        # add condition for omf
        
        self.ioda.omf = self.ioda.ObsValue['aerosolOpticalDepth'] - self.ioda.HofX['aerosolOpticalDepth']
        print('type omf', self.ioda.omf)
        if (self.ioda.dims['Location'] == 0):
            raise ValueError("No observation %s, nothing to do" % self.ioda.dims['Location']) 
        self.nobs, self.nchannels  = self.ioda.dims['Location'], self.ioda.dims['Channel']
        print('nobs', self.nobs)
        self.sqc = yaml.safe_load(config)
        
        self.qcexcl = self.ioda.PreQc['aerosolOpticalDepth'].data
        self.qchist = self.ioda.HistQc['aerosolOpticalDepth'].data

        # use lev in buddy check (original code has the wavelength) but BC could be used for vertically resolved 
        # obs 
        obs_wavelength = self.ioda.MetaData['obs_wavelength'].data
        print('obs wavelength', obs_wavelength)
        # reshape lev as (nlocs, nch) for the buddy check, repeating the same value for each obs if 2D variable
        self.lev = np.tile(obs_wavelength, self.nobs)
        self.lons, self.lats = self.ioda.MetaData['longitude'], self.ioda.MetaData['latitude']
        xobs, yobs, zobs = self.ll2xyz(self.lons, self.lats)
        self.xobs = np.array(xobs, dtype=np.float64)
        self.yobs = np.array(yobs, dtype=np.float64)
        self.zobs  = np.array(zobs, dtype=np.float64)
   
        self.I_old = ((self.qcexcl==0) & (self.qchist==17))
        print('I old reaccept', self.I_old, self.I_old.shape)
        self.reac_old = np.where((self.qcexcl==0) & (self.qchist==17))[0]
        print('find index reaccepted values in old code:', np.where((self.qcexcl==0) & (self.qchist==17))[0])
   #     print('parameters for comparisons:', 'lat =', self.ioda.MetaData['latitude'][self.reac_old].data, ' lon = ', self.ioda.MetaData['longitude'][self.reac_old].data,
   #           'omf =', self.ioda.omf[self.reac_old].data, 'obs value=', self.ioda.ObsValue['aerosolOpticalDepth'][self.reac_old].data, 'xobs = ', self.xobs[self.reac_old], 'qchi= ', self.qchist[self.reac_old])
    
        self._reset_qc()
        self._getErrVar()
        self._backgCheck()
     
        
   #     print('find index reaccepted values in old code after background:', np.where((self.qcexcl==0) & (self.qchist==17))[0])
   #     print('parameters for comparisons after background check:', 'lat =', self.ioda.MetaData['latitude'][self.reac_old].data, ' lon = ', self.ioda.MetaData['longitude'][self.reac_old].data,
   #          'omf =', self.ioda.omf[self.reac_old].data, 'obs value=', self.ioda.ObsValue['aerosolOpticalDepth'][self.reac_old].data, 'xobs = ', self.xobs[self.reac_old], 'qchi= ', self.qchist[self.reac_old])
             
        self._reorder_arrays(message="Reordering by exclusion status")
   
        self.regions, self.nmaxregions, self.iregbeg, self.ireglen = self._getregions()
#        print('regions', self.regions)
        self.reac_old = np.where(self.I_old == True)[0]
        print('find index reaccepted values in old code after region reorder:', np.where((self.qcexcl==0) & (self.qchist==17))[0])
 #       print('parameters for comparisons after region binning:', 'lat =', self.ioda.MetaData['latitude'][self.reac_old].data, ' lon = ', self.ioda.MetaData['longitude'][self.reac_old].data,
 #            'omf =', self.ioda.omf[self.reac_old].data, self.regions[self.reac_old], 'obs value=', 
 #             self.ioda.ObsValue['aerosolOpticalDepth'][self.reac_old].data, 'xobs = ', self.xobs[self.reac_old], 'qchi= ', self.qchist[self.reac_old])
        
        self._buddyCheck()

    def _reset_qc(self):
        if self.sqc['reset_allqc']:
            self.qcexcl[:] = 0
            self.qchist[:] = 0
            print('[x] reset all exclusion and history marks')
        elif self.sqc['reset_sqc']:
            self.qchist[:] = 0
            print('[x] reset all history marks')
                
    def _reorder_arrays(self,sort_indices=None, sort_by_key=None, descend=False, message=None):
        """
        Reorder all arrays based on provided indices or a key to sort by.
    
        Parameters:
        -----------
        sort_indices : numpy.ndarray, optional
            Pre-computed indices to use for reordering. If None, will be computed from sort_by_key.
        sort_by_key : numpy.ndarray, optional
            Array to use as the key for sorting. If None and sort_indices is None, will sort by qcexcl.
        descend : bool, optional
            Whether to sort in descending order. Default is False.
        message : str, optional
            Message to print before and after reordering. If None, a default message is used.  
        Returns:
        --------
        Number of valid observations if sorting by exclusion, otherwise None.
        """
        if message is None:
           if sort_by_key is None and sort_indices is None:
              message = f'Reordering by exclusion status'
           elif sort_by_key is not None:
              message = f'Reordering by provided key'
           else:
              message = f'Reordering by provided indices'
    
        print(f'[x] {message}')

        test_index = self.reac_old[0]  # Take first reaccepted observation as test case
        print("BEFORE REORDERING:")
        print(f"Test index {test_index}:")
        print(f"  Lat: {self.ioda.MetaData['latitude'][test_index].data}")
        print(f"  Lon: {self.ioda.MetaData['longitude'][test_index].data}")
        print(f"  OMF: {self.ioda.omf[test_index].data}")
        print(f"  XObs: {self.xobs[test_index]}")

        # Store omf array before reordering
        if hasattr(self.ioda, 'omf'):
           omf_array = self.ioda.omf.copy()
           print(f"  OMF: {omf_array[test_index].data}")
        else:
           # If omf is stored elsewhere, adjust this accordingly
           omf_array = None
           print("  OMF: Not found as direct attribute")
        
        # Determine the sort indices
        if sort_indices is None:
           if sort_by_key is None:
              # Default: sort by exclusion status (non-excluded obs first across all channels)
              mask = np.all(self.qcexcl == 0, axis=1)
              sort_indices = np.argsort(~mask)   #(True values first, then False)
              n_valid = np.sum(mask)             # Count number of non-excluded obs
           else:
              # Sort by the provided key
              indx = np.arange(len(sort_by_key), dtype=np.int32)
              sort_indices = indx[np.argsort(sort_by_key[indx], kind='mergesort')[::-1 if descend else 1]]
              n_valid = None
        else:
           n_valid = None
            
        # Track where our test index moves to
        new_position = np.where(sort_indices == test_index)[0][0]
        print(f"After reordering, index {test_index} will move to position {new_position}")

        # Create new DataTree
        new_tree = xr.DataTree()

        # Process each top-level group
        for group_name, group in self.ioda.items():
        # Create a corresponding group in the new tree
           new_tree[group_name] = xr.DataTree()    
           # Process variables in this group
           for var_name, var in group.items():
               if isinstance(var, xr.DataArray) and 'Location' in var.dims:
                  # Reorder arrays with Location dimension
                  new_tree[group_name][var_name] = var.isel(Location=sort_indices)
               else:
                  # Copy other arrays directly
                  new_tree[group_name][var_name] = var

        # Replace the original tree
        self.ioda = new_tree
      
        # Reorder qcexcl and qchist
        self.qcexcl = self.qcexcl[sort_indices]
        self.qchist = self.qchist[sort_indices]
        self.I_old = self.I_old[sort_indices]
        
        
        # Reorder VarO and VarF if they exist
        if hasattr(self, 'VarO'):
           self.VarO = self.VarO[sort_indices]
        if hasattr(self, 'VarF'):
           self.VarF = self.VarF[sort_indices]

        # Reorder coordinate arrays if they exist
        for attr in ['xobs', 'yobs', 'zobs']:
            if hasattr(self, attr):
               setattr(self, attr, getattr(self, attr)[sort_indices])
        # Reattach omf after reordering if it existed
        if omf_array is not None:
        # Reorder omf array
           self.ioda.omf = omf_array.isel(Location=sort_indices)

        # After reordering
        print("AFTER REORDERING:")
        print(f"Original test index {test_index} is now at {new_position}:")
        print(f"  Lat: {self.ioda.MetaData['latitude'][new_position].data}")
        print(f"  Lon: {self.ioda.MetaData['longitude'][new_position].data}")
        print(f"  OMF: {self.ioda.omf[new_position].data}")
        print(f"  XObs: {self.xobs[new_position]}")
       
    
        if n_valid is not None:
           print(f'[x] {message} complete: {n_valid} valid observations moved to the front')
        else:
           print(f'[x] {message} complete')

         # Check if our boolean mask tracking works
  #      reac_indices = np.where(self.I_old)[0]
 #       if len(self.I_old.shape) > 1:
  #      # Flatten or use the first column
  #        reac_indices = np.where(self.I_old[:, 0])[0]
  #      else:
  #        reac_indices = np.where(self.I_old)[0]
  #      print(f"Reaccepted indices after reordering: {reac_indices[:5]}...")
  #      print("First reaccepted observation values:")
  #      if len(reac_indices) > 0:
  #         idx = reac_indices[0]
  #         print(f"  Lat: {self.ioda.MetaData['latitude'][idx].data}")
  #         print(f"  Lon: {self.ioda.MetaData['longitude'][idx].data}")
  #         print(f"  OMF: {self.ioda.omf[idx].data}")
  #         print(f"  XObs: {self.xobs[idx]}")
      
        # Return the number of valid observations for future reference
        return n_valid
  

    def _getErrVar(self):
        print(f'[x] Getting Error Variances')
        self.sigs = yaml.safe_load(LAOD_sigs if self.log else AOD_sigs)
        
        self.VarO = np.full((self.ioda.dims['Location'], self.ioda.dims['Channel']), -999.)
        self.VarF = np.full((self.ioda.dims['Location'], self.ioda.dims['Channel']), -999.)
        
        for ch, channel in enumerate(self.ioda.MetaData.obs_wavelength):
            found = False
            for wav in self.sigs:
                if abs(channel.data - float(wav)) < 0.01:
                    self.VarO[:, ch] = self.sigs[wav]['sigO']**2
                    self.VarF[:, ch] = self.sigs[wav]['sigF']**2
                    found = True
                    break
            if not found:
                print(f'[x] Warning: no sigO sigF provided for the observed wavelength {channel}')

    def _backgCheck(self):
        var = self.VarF + self.VarO
        if np.any(var <= 0):
           raise ValueError("Prescribed O-F Variance is zero")
    
        # Calculate tau for all channels
        tau = np.sqrt(self.ioda.omf**2 / var)

        # Extreme outliers
        extreme_outliers = tau > self.sqc['Backg_Check']['tau_bgx']
        self.qcexcl[extreme_outliers] = 21
    
        # Suspect outliers
        suspect_outliers = tau > self.sqc['Backg_Check']['tau_bgh']
        self.qchist[suspect_outliers] = 17
    
        # Count exclusions and suspects for each channel
        nexcl_per_channel = np.sum(extreme_outliers, axis=0)
        nsuspect_per_channel = np.sum(suspect_outliers, axis=0)
    
        for channel in range(nexcl_per_channel.size):
            nexcl = nexcl_per_channel[channel].item()  # Convert to Python scalar
            nsuspect = nsuspect_per_channel[channel].item()  # Convert to Python scalar
            print(f'[x] number of exclusions after background_check for channel {channel} is {nexcl}')
            print(f'[x] number of suspects after background_check for channel {channel} is {nsuspect}')

        print('shape qc after background check:', self.qchist.shape)

        # Return suspect_outliers for use in buddy check
        return suspect_outliers
        
    @staticmethod
    def ll2xyz(lon, lat): # for one obs lat-lon, the original code loop over nobs to create vectors
        """Convert longitude and latitude to Cartesian coordinates (x, y, z) on the unit sphere."""
        xobs = np.cos(np.radians(lat)) * np.cos(np.radians(lon))
        yobs = np.cos(np.radians(lat)) * np.sin(np.radians(lon))
        zobs = np.sin(np.radians(lat))
        print('xobs', xobs.shape)
        return xobs, yobs, zobs

    def _getregions(self): 
        import BuddyCheck_
        regions = np.zeros(self.nobs, dtype=np.int32)
              
   #     print(BuddyCheck_.py_icosahedron_regions.__doc__)
        
        (regions, nmaxregions) = BuddyCheck_.py_icosahedron_regions(self.xobs, self.yobs, self.zobs)
        # Sort by region (the size of regions is nobs for the sorting and reordering compared to 
        # original fortran code)

        # Store a copy of regions before reordering
        original_regions = regions.copy()
      
        self._reorder_arrays(sort_by_key=regions, message="Sorting observations by region")
  #      print('regions',regions)
        # Reorder the regions array to match the other arrays
       
        indx = np.arange(len(original_regions), dtype=np.int32)
        sort_indices = indx[np.argsort(original_regions[indx], kind='mergesort')]
        regions = original_regions[sort_indices]
        print('nregions', nmaxregions) 
        
         # Set integer region pointers      
        iregbeg = np.zeros(nmaxregions, dtype=np.int32)  # beginning of a region
        ireglen = np.zeros(nmaxregions, dtype=np.int32)  # number of obs in a region
    
        # Count observations in each region
        np.add.at(ireglen, regions - 1, 1)
    
        # Set region beginning indices
        iregbeg[1:] = np.cumsum(ireglen[:-1]) # 0 based for python
        print('region beg and len', iregbeg, ireglen)
   #     print('number of regions with obs', len(np.unique(regions))) # number of regions with data
        
        return regions, nmaxregions, iregbeg, ireglen

#    def _issuspect(self):
    def _buddyCheck(self):
        import BuddyCheck_
       
        print(f'[x] Starting the buddy check')
        
        search_rad = self.sqc['Buddy_Check']['search_rad']  # unit length scale
        ls_h = self.sqc['Buddy_Check']['ls_h']  # 0.15e6 (m) in original code -> combined 150km search
        ls_v = self.sqc['Buddy_Check'].get('ls_v', 1000.0)  # Vertical length scale
        single_level = self.sqc['Buddy_Check'].get('single_level', True)
        nbuddy_max = self.sqc['Buddy_Check'].get('nbuddy_max', 100)
        seplim = self.sqc['Buddy_Check'].get('seplim', 26.5)  # Separation limit between regions

        # Initialize arrays for buddy check results
        reaccept = np.zeros((self.nobs, self.nchannels), dtype=bool)
    
        for channel in range(self.nchannels):
            print(f'[x] Processing wavelength {channel}')
            
            # first list of suspect after the background check
            suspect = (self.qchist[:, channel] > 0) & (self.qcexcl[:, channel] == 0)
            print(f"Suspect indices for channel {channel}:", np.where(suspect)[0])

            for iteration in range(1, self.sqc['Buddy_Check']['niter_max'] + 1): 
                print(f'[x] Iteration {iteration} for channel {channel}')
            
                # Skip if no suspect observations
                if not np.any(suspect):
                   print(f'[x] No suspect observations for channel {channel} at iteration {iteration}')
                   break     
                    
                # Vectorized preparation of suspect data
                ki_susp = []
                kr_susp = []

                # Create arrays to hold region indices for each observation
                obs_regions = np.zeros(self.nobs, dtype=int)

                # Populate the obs_regions array
                for ireg in range(1, self.nmaxregions + 1):
                    if self.ireglen[ireg-1] > 0:
                       ibeg = self.iregbeg[ireg-1]
                       iend = ibeg + self.ireglen[ireg-1] - 1
                       obs_regions[ibeg:iend+1] = ireg

                # Find indices where I_old is True (reaccepted observations)
                #reaccepted_indices = np.where(self.I_old)[0]
                #print(f"Reaccepted indices old code before BC: {reaccepted_indices[:5]}...")

                
                # Find which of the reaccepted observations are also suspect
                #suspect_mask = suspect[reaccepted_indices]
                #suspect_indices = reaccepted_indices[suspect_mask]
                suspect_indices = np.where(suspect)[0]
                print(f"Suspect indices before BC: {suspect_indices[:5]}...")
                #suspect_indices = np.where(suspect)[0]
                # print(f"Suspect indices: {suspect_indices[:5]}...")

                # Get the regions for these suspect observations
                ki_susp = suspect_indices  + 1 # These are the actual indices in the dataset, + 1 if index 0 in python -> 1 in fortran, f2py handles the return of reaccept indices with -1
                kr_susp = obs_regions[suspect_indices]  # These are the region numbers

                print(f"ki_susp shape: {ki_susp.shape}")
                print(f"First few ki_susp: {ki_susp[:10]}")
                print(f"First few kr_susp: {kr_susp[:10]}")

                n_susp = len(ki_susp)
                print(f'[x] number of suspects for channel {channel} at iteration {iteration} is {n_susp}')
                # note: lev is pressure level of obs--> AOD single level will put lev arrays to 1 for now, original
                # code has the wavelengths in it
                self.iregbeg = self.iregbeg + 1 # for fortran based index starting at 1
                print('beg', self.iregbeg)
                #print('before buddy check call', 'lat =', self.ioda.MetaData['latitude'][self.reac_old].data, ' lon = ', self.ioda.MetaData['longitude'][self.reac_old].data)
     #           'omf =', self.ioda.omf[self.reac_old].data, self.regions[self.reac_old], 'xobs = ', self.xobs[self.reac_old], 'yobs = ',self.yobs[self.reac_old], 
     #           'ki_susp', ki_susp, 'kr_susp', kr_susp)
                lats, lons = self.ioda.MetaData['latitude'].data, self.ioda.MetaData['longitude'].data    
                reaccept[:, channel] = BuddyCheck_.py_find_buddies(ki_susp, kr_susp, self.xobs, 
                        self.yobs, self.zobs, lats, lons, self.lev, self.ioda.omf,
                        self.VarF, self.VarO, self.qcexcl, ls_h, ls_v, search_rad, single_level, nbuddy_max, self.iregbeg, self.ireglen, 
                        seplim)
             
                print("reaccept:", reaccept[0:10,0])
                reaccepted_obs = np.where(reaccept[:,channel])[0]
                if len(reaccepted_obs) > 0:
                    print(f"[x] Reaccepted {len(reaccepted_obs)} observations for channel {channel}")
                    # Update suspect list for next iteration
                    suspect[reaccept[:,channel]] = False
                else:
                    print(f"[x] No observations reaccepted for channel {channel} at iteration {iteration}")
                    break  # No observations reaccepted, so stop iterating
                
         
    print(f'[x] Buddy check completed')

          
        












