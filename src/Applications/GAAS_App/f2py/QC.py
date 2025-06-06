#!/usr/bin/env python

import xarray as xr
import numpy as np
import yaml
import sys

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
 440:   
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

    def __init__(self, iodaFiles, config=None, log=True, transform_log=True, verbose=True):
        if config is None:
            config = SQC_Param

        self.log = log
        self.transform = transform_log

        self.verbose = verbose
         
        if isinstance(iodaFiles, xr.Dataset):
            self.ioda = iodaFiles
        else:
            self.ioda = xr.open_datatree(iodaFiles)

        # Check if the 'omf' attribute exists in self.ioda
        if not hasattr(self.ioda, 'omf'):
           # If it doesn't exist, create it by computing the difference
            try:
                if self.transform:  
                   self.ioda.omf = np.log(self.ioda.ObsValue['aerosolOpticalDepth']+0.01) - np.log(1000*self.ioda.hofx['aerosolOpticalDepth']+0.01)
                   print('ok', self.ioda.ObsValue['aerosolOpticalDepth'],self.ioda.hofx['aerosolOpticalDepth']) 
                else:
                   self.ioda.omf = self.ioda.ObsValue['aerosolOpticalDepth'] - self.ioda.hofx['aerosolOpticalDepth']  #already in log in file 
            except (AttributeError, KeyError) as e:
                print(f"Couldn't create omf field: {e}")
       
        # Handle the case where the required attributes/keys don't exist
        if (self.ioda.dims['Location'] == 0):
            raise ValueError("No observation %s, nothing to do" % self.ioda.dims['Location']) 
        self.nobs, self.nchannels  = self.ioda.dims['Location'], self.ioda.dims['Channel']
        print(f" [x]: The observation file had {self.nobs} observations and {self.nchannels} wavelengths")
        
        # Get the sqc parameters from yaml file
        self.sqc = yaml.safe_load(config)
        
        # Get QCexcl and QChist
        self.qcexcl = self.ioda.PreQc['aerosolOpticalDepth'].data
        self.qchist = self.ioda.HistQc['aerosolOpticalDepth'].data

        # use lev in buddy check (original code has the wavelength) but BC could be used for vertically resolved obs 
        self.obs_wavelength = self.ioda.MetaData['obs_wavelength'].data
        print('obs wavelength', self.obs_wavelength)
        # reshape lev as (nlocs, nch) for the buddy check, repeating the same value for each obs if 2D variable
        self.lev = np.tile(np.array([self.obs_wavelength]), (self.nobs, 1))
        print('lev values', self.lev.shape)        

        # Compute the cartesians coordinates
        self.lons, self.lats = self.ioda.MetaData['longitude'], self.ioda.MetaData['latitude']
        xobs, yobs, zobs = self.ll2xyz(self.lons, self.lats)
        self.xobs = np.array(xobs, dtype=np.float64)
        self.yobs = np.array(yobs, dtype=np.float64)
        self.zobs  = np.array(zobs, dtype=np.float64)
   
        # Reset the QC flag if needed according to yaml
        self._reset_qc()

        # Get Obs and background variances
        self._getErrVar()

        # First, select the observation with a Background check - statistical outliers, qchist and qcexcl are updated
        self._backgCheck()
     
        # Reordering foloowing the background check with non-excluded data moved to the front      
        self._reorder_arrays(message="Reordering by exclusion status")
   
        # Bining the observation using the icosahedron geometry regions
        self.regions, self.nmaxregions, self.iregbeg, self.ireglen = self._getregions()

        # Finally doing the buddy check using suspect observations  
        self._buddyCheck()

        # Updating the QC values in the IODA file

#---------------  
    def _reset_qc(self):

        if self.sqc['reset_allqc']:
            self.qcexcl[:] = 0
            self.qchist[:] = 0
            print('[x] reset all exclusion and history marks')
        elif self.sqc['reset_sqc']:
            self.qchist[:] = 0
            print('[x] reset all history marks')
#---------------                
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

        test_index = 2  # Take a supsect observation as test case
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
           # If omf is stored elsewhere
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

        # Return the number of valid observations for future reference
        return n_valid
#-------------  
    def _getErrVar(self):
        print(f'[x] Getting Error Variances')

        # May add an option if VarO is available in ioda file
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
                print(f'[x] Warning: no sigO sigF provided for the observed wavelength {self.obs_wavelength[channel]}')
#-------------
    def _backgCheck(self):

        """ 
        Checks for statistical outliers:
        Parameters: - Obs and background variances ( self.VarO and self.VarF ) 
        Returns: suspect (qch = 17) and rejected obs (qcx = 21)   

        Marks observations for which the squared omf residual exceeds a pre-defined
        multiple of the presumed variance of the residual defined in yaml (tau_bgx and tau_bgh)
        """
        
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
            print(f'[x] number of exclusions after background_check for channel {self.obs_wavelength[channel]} is {nexcl}')
            print(f'[x] number of suspects after background_check for channel {self.obs_wavelength[channel]} is {nsuspect}')

        print('shape qc after background check:', self.qchist.shape)

#---------------        
    @staticmethod
    def ll2xyz(lon, lat): 
        """Convert longitude and latitude to Cartesian coordinates (x, y, z) on the unit sphere."""
        xobs = np.cos(np.radians(lat)) * np.cos(np.radians(lon))
        yobs = np.cos(np.radians(lat)) * np.sin(np.radians(lon))
        zobs = np.sin(np.radians(lat))
        return xobs, yobs, zobs

#---------------
    def _getregions(self): 
        """ Uses Icosahedron fortran subroutine to return the region number, max regions and  
        the index of the start of a region and the number of obs in each region
        """
        
        import BuddyCheck_
        regions = np.zeros(self.nobs, dtype=np.int32)
              
        (regions, nmaxregions) = BuddyCheck_.py_icosahedron_regions(self.xobs, self.yobs, self.zobs)
        
        # Sort the datatree and all parameters(qc, Var.. ) by region (the size of regions is nobs for the sorting and reordering)

        # Store a copy of regions before reordering
        original_regions = regions.copy()
      
        self._reorder_arrays(sort_by_key=regions, message="Sorting observations by region")
        
        # Reorder the regions array to match the other arrays
        indx = np.arange(len(original_regions), dtype=np.int32)
        sort_indices = indx[np.argsort(original_regions[indx], kind='mergesort')]
        regions = original_regions[sort_indices]
        print(f'[x] Total number of regions is {nmaxregions} and the number of regions with obs is {len(np.unique(regions))}') 
        
        # Set integer region pointers      
        iregbeg = np.zeros(nmaxregions, dtype=np.int32)  # beginning of a region
        ireglen = np.zeros(nmaxregions, dtype=np.int32)  # number of obs in a region
    
        # Count observations in each region
        np.add.at(ireglen, regions - 1, 1)
    
        # Set region beginning indices
        iregbeg[1:] = np.cumsum(ireglen[:-1]) # 0 based for python
        
        return regions, nmaxregions, iregbeg, ireglen
#----------------
    def _buddyCheck(self):
        import BuddyCheck_
        """
        Checks each suspect observation (with any nonzero history mark) against its "buddies" 
        (nearby observations of the same variable with a zero history mark). If a suspect observation passes the buddy check, it joins the 
        pool of buddies for the next round. This procedure is iterated until the pool stabilizes, or until a maximum number of iterations
        is reached.

        """
        print(f'[x] Starting the buddy check')
        
        search_rad = self.sqc['Buddy_Check']['search_rad']  # unit length scale
        ls_h = self.sqc['Buddy_Check']['ls_h']  # 0.15e6 (m) in original code -> combined 150km search
        ls_v = self.sqc['Buddy_Check'].get('ls_v', 1000.0)  # Vertical length scale
        single_level = self.sqc['Buddy_Check'].get('single_level', True)
        tau_buddy = self.sqc['Buddy_Check'].get('tau_buddy', 0.1) # tolerance parameter: 0.1->stringent the obs O-F needs to be close to the predicted O-F from the buddies
        nbuddy_max = self.sqc['Buddy_Check'].get('nbuddy_max', 100)
        seplim = self.sqc['Buddy_Check'].get('seplim', 26.5)  # Separation limit between regions

        # Initialize arrays for buddy check results
        reaccept = np.zeros((self.nobs, self.nchannels), dtype=bool)
    
        for channel in range(self.nchannels):
            print(f'[x] Processing wavelength {self.obs_wavelength[channel]}')
            
            # first list of suspect after the background check
            suspect = (self.qchist[:, channel] > 0) & (self.qcexcl[:, channel] == 0)
            print(f"Suspect indices for channel {self.obs_wavelength[channel]}:", np.where(suspect)[0])

            for iteration in range(1, self.sqc['Buddy_Check']['niter_max'] + 1): 
                print(f'[x] Iteration {iteration} for channel {self.obs_wavelength[channel]}')
            
                # Skip if no suspect observations
                if not np.any(suspect):
                   print(f'[x] No suspect observations for channel {self.obs_wavelength[channel]} at iteration {iteration}')
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
                
                self.iregbeg = self.iregbeg + 1 # for fortran based index starting at 1
                suspect_indices = np.where(suspect)[0]
                print(f"Suspect indices before BC: {suspect_indices[:5]}...")

                # Get the regions for these suspect observations
                ki_susp = suspect_indices  + 1 # These are the actual indices in the dataset, + 1 if index 0 in python -> 1 in fortran, f2py handles the return of reaccept indices with -1
                kr_susp = obs_regions[suspect_indices]  # These are the region numbers

                print(f"ki_susp shape: {ki_susp.shape}")
                print(f"First few ki_susp: {ki_susp[:10]}")
                print(f"First few kr_susp: {kr_susp[:10]}")

                n_susp = len(ki_susp)
                print(f'[x] number of suspects for channel {self.obs_wavelength[channel]} at iteration {iteration} is {n_susp}')
                # note: lev is pressure level of obs--> AOD single level will put lev arrays to 1 for now, original
                # code has the wavelengths in it
                
                lats, lons = self.ioda.MetaData['latitude'].data, self.ioda.MetaData['longitude'].data    

                reaccept[:, channel] = BuddyCheck_.py_find_buddies(ki_susp, kr_susp, self.xobs, 
                        self.yobs, self.zobs, lats, lons, self.lev[:,channel], self.ioda.omf[:,channel],
                        self.VarF[:,channel], self.VarO[:,channel], self.qcexcl[:,channel], ls_h, ls_v, search_rad, tau_buddy,
                        single_level, nbuddy_max, self.iregbeg, self.ireglen,seplim)
             
                print("reaccept:", reaccept[0:10,0])
                reaccepted_obs = np.where(reaccept[:,channel])[0]
                if len(reaccepted_obs) > 0:
                    print(f"[x] Reaccepted {len(reaccepted_obs)} observations for channel {self.obs_wavelength[channel]}")
                    # Update suspect list for next iteration
                    suspect[reaccept[:,channel]] = False
                else:
                    print(f"[x] No observations reaccepted for channel {self.obs_wavelength[channel]} at iteration {iteration}")
                    break  # No observations reaccepted, so stop iterating

            # Update QC flags after BC, all supected obs not reaccepted have QCX = 17 
            self.qcexcl[suspect, channel] = 17
    print(f'[x] Buddy check completed')

          
        












