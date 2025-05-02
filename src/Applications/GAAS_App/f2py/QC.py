#!/usr/bin/env python

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
  niter_max: 2
  search_rad: 1
  nbuddy_max: 100
  nstar: 0
  ls_h: 150000.0
  ls_v: 100.0     
  seplim: 26.5  
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
        self.lev = np.tile(obs_wavelength, (self.nobs,1))
        print('lev', self.lev.shape,self.lev[0:5,0:1])

        self._reset_qc()
        self._getErrVar()
        self._backgCheck()
        self._reorder_arrays(message="Reordering by exclusion status")
        # next regions divisions, call fortran code
        lons, lats = self.ioda.MetaData['longitude'], self.ioda.MetaData['latitude']
        xobs, yobs, zobs = self.ll2xyz(lons, lats)
        self.xobs = np.array(xobs, dtype=np.float64)
        self.yobs = np.array(yobs, dtype=np.float64)
        self.zobs  = np.array(zobs, dtype=np.float64)
        self.regions, self.nmaxregions, self.iregbeg, self.ireglen = self._getregions()
        print('regions', self.regions)
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
        

        # Function to reorder a single DataArray
        def reorder_dataarray(da):
           if 'Location' in da.dims:
               return da.isel(Location=sort_indices)
           return da

        # Reorder the entire DataTree
        reordered_tree = xr.DataTree()
        for path, data_array in self.ioda.items():
            if isinstance(data_array, xr.DataArray):
               reordered_tree[path] = reorder_dataarray(data_array)
            else:
               reordered_tree[path] = data_array

        # Handle omf separately if it exists
        if hasattr(self.ioda, 'omf'):
           omf = reorder_dataarray(self.ioda.omf)
        else:
           omf = None

        # Replace the original DataTree with the reordered one
        self.ioda = reordered_tree

        # Add omf back to the reordered DataTree if it existed
        if omf is not None:
           self.ioda.omf= omf
            
        # Reorder qcexcl and qchist
        self.qcexcl = self.qcexcl[sort_indices]
        self.qchist = self.qchist[sort_indices]
        
        # Reorder VarO and VarF if they exist
        if hasattr(self, 'VarO'):
           self.VarO = self.VarO[sort_indices]
        if hasattr(self, 'VarF'):
           self.VarF = self.VarF[sort_indices]

        # Reorder coordinate arrays if they exist
        for attr in ['xobs', 'yobs', 'zobs']:
            if hasattr(self, attr):
               setattr(self, attr, getattr(self, attr)[sort_indices])
    
        if n_valid is not None:
           print(f'[x] {message} complete: {n_valid} valid observations moved to the front')
        else:
           print(f'[x] {message} complete')
        
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
              
        print(BuddyCheck_.py_icosahedron_regions.__doc__)
        
        (regions, nmaxregions) = BuddyCheck_.py_icosahedron_regions(self.xobs, self.yobs, self.zobs)
        # Sort by region (the size of regions is nobs for the sorting and reordering compared to 
        # original fortran code)
        self._reorder_arrays(sort_by_key=regions, message="Sorting observations by region")
        print('regions',regions)
        print('nregions', nmaxregions) 
        
         # Set integer region pointers      
        iregbeg = np.zeros(nmaxregions, dtype=np.int32)  # beginning of a region
        ireglen = np.zeros(nmaxregions, dtype=np.int32)  # number of obs in a region
    
        # Count observations in each region
        np.add.at(ireglen, regions - 1, 1)
    
        # Set region beginning indices
        iregbeg[0] = 1
        iregbeg[1:] = np.cumsum(ireglen[:-1]) + 1 # (for fortran compatibility)
        print('region beg and len', iregbeg, ireglen)
        print('number of regions with obs', len(np.unique(regions))) # number of regions with data
        
        return regions, nmaxregions, iregbeg, ireglen

    def _buddyCheck(self):
        import BuddyCheck_
       
        print(BuddyCheck_.py_find_buddies.__doc__)
        print(f'[x] Starting the buddy check')
        
        search_rad = self.sqc['Buddy_Check']['search_rad']  # unit length scale
        ls_h = self.sqc['Buddy_Check']['ls_h']  # 0.15e6 (m) in original code -> combined 150km search
        ls_v = self.sqc['Buddy_Check'].get('ls_v', 0.0)  # Vertical length scale
        single_level = self.sqc['Buddy_Check'].get('single_level', False)
        nbuddy_max = self.sqc['Buddy_Check'].get('nbuddy_max', 100)
        seplim = self.sqc['Buddy_Check'].get('seplim', 26.5)  # Separation limit between regions

        # Initialize arrays for buddy check results
        reaccept = np.zeros((self.nobs, self.nchannels), dtype=bool)
    
        for channel in range(self.nchannels):
            print(f'[x] Processing wavelength {channel}')

            # note: after the reorderring(reorder_dataarray), the dimensions are now stored within each group rather 
            # than at the top level of the datatree
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

                # Find all suspect observations
                suspect_indices = np.where(suspect)[0]
                # will add + 1 to handle not having 0 in fortran

                ki_susp = suspect_indices + 1 
                print('ki_susp', ki_susp[0:10])
                kr_susp = obs_regions[suspect_indices] # keep using 0-based when accessing python array

                n_susp = len(ki_susp)
                print(f'[x] number of suspects for channel {channel} at iteration {iteration} is {n_susp}')
                # note: lev is pressure level of obs--> AOD single level will put lev arrays to 1 for now, original
                # code has the wavelengths in it
        
                print("ki_susp:", type(ki_susp), ki_susp.shape)
                print("kr_susp:", type(kr_susp), kr_susp.shape)
                print("xobs:", type(self.xobs), self.xobs.shape)
                print("yobs:", type(self.yobs), self.yobs.shape)
                print("zobs:", type(self.zobs), self.zobs.shape)
                print("lev:", type(self.lev), self.lev.shape)
                print("omf:", type(self.ioda.omf), self.ioda.omf.shape)
                print("VarF:", type(self.VarF), self.VarF.shape)
                print("VarO:", type(self.VarO), self.VarO.shape)
                print("ls_h:", type(ls_h))
                print("ls_v:", type(ls_v))
                print("search_rad:", type(search_rad), search_rad)
                print("single_level:", type(single_level))
                print("nbuddy_max:", type(nbuddy_max), nbuddy_max)
                print("iregbeg:", type(self.iregbeg), self.iregbeg.shape)
                print("ireglen:", type(self.ireglen), self.ireglen.shape)
                print("seplim:", type(seplim), seplim)
                reaccept = BuddyCheck_.py_find_buddies(ki_susp, kr_susp, self.xobs, self.yobs, self.zobs, self.lev, self.ioda.omf,
                        self.VarF, self.VarO, ls_h, ls_v, search_rad, single_level, nbuddy_max, self.iregbeg, self.ireglen, 
                        seplim)
            
                print('reaccept', np.where(reaccept==True))         
    print(f'[x] Buddy check completed')

          
