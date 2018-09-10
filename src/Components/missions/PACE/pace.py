#!/usr/bin/env python

"""
    Utility to create GEOS-5 Collections on PACE L1B granule.
    Based on trj_sampler.py and geo_sampler.py

    Arlindo da Silva, November 2014.
    P. Castellanos, Aug 2016 - Modified from GEO to LEO
    P. Castellanos, 2018 - Modified from MODIS to PACE
"""

import os

from   datetime        import datetime, timedelta
from   dateutil.parser import parse  as isoparser

import numpy as np

from   netCDF4         import Dataset
from   glob            import glob



SDS = {'longitude'  : 'geolocation_data',
       'latitude'   : 'geolocation_data',
       'ev_mid_time': 'scan_line_attributes' }


ALIAS = {'longitude': 'lon',
         'latitude' : 'lat',
         'ev_mid_time': 'midTime'}


class PACE(object):
    """
    Generic container for PACE SDS
    """

    def __init__(self, Path,verb=False):

        # Read each granule, appending them to the list
        # ---------------------------------------------
        if type(Path) is list:
            if len(Path) == 0:
               self.nobs = 0
               print "WARNING: Empty PACE_L1B object created"
               return
            else:
                self.nobs = len(Path)
        else:
           Path = [Path, ]
           self.nobs = 1

        # Create empty lists for SDS to be read from file
        # -----------------------------------------------
        self.granules = Path
        self.verb = verb
        self.SDS  = SDS.keys()
        self.SDSg = SDS

        for name in self.SDS:
           self.__dict__[name] = []

        self.scanStart = []

        # Read each granule, appending them to the list
        # ---------------------------------------------
        self._readList(Path)

        # Alias
        for sds in self.SDS:
            self.__dict__[ALIAS[sds]] = self.__dict__[sds]

        # Create corresponding python time
        # --------------------------------
        nscan,npixel = self.lon[0].shape
        self.tyme = []        
        for scanStart,midTime,lon in zip(self.scanStart,self.midTime,self.lon):
            scanStart  = isoparser(scanStart.strftime('2006-%m-%dT00:00:00'))
            tyme       = np.array([scanStart + timedelta(seconds=t) for t in midTime])    
            tyme.shape = (nscan,1)
            tyme       = np.repeat(tyme,npixel,axis=1)
            tyme       = np.ma.array(tyme)
            tyme.mask  = lon.mask
            self.tyme.append(tyme)
                    

    def _readList(self,List):
        """
        Recursively, look for files in list; list items can
        be files or directories.
        """
        for item in List:
            if os.path.isdir(item):      self._readDir(item)
            elif os.path.isfile(item):   self._readGranule(item)
            else:
                print "%s is not a valid file or directory, ignoring it"%item

#---
    def _readDir(self,dir):
        """Recursively, look for files in directory."""
        for item in os.listdir(dir):
            path = dir + os.sep + item
            if os.path.isdir(path):      self._readDir(path)
            elif os.path.isfile(path):   self._readGranule(path)
            else:
                print "%s is not a valid file or directory, ignoring it"%item

#---
    def _readGranule(self,filename):
        """Reads one PACE granule."""

        # Don't fuss if the file cannot be opened
        # ---------------------------------------
        try:
            if self.verb:
                print "[] Working on "+filename
            nc = Dataset(filename)
        except:
            if self.verb > 2:
                print "- %s: not recognized as an netCDF file"%filename
            return 

        # Read select variables (do not reshape)
        # --------------------------------------
        self.scanStart.append(isoparser(nc.time_coverage_start))
        for sds in self.SDS:            
            group = self.SDSg[sds]
            v = nc.groups[group].variables[sds][:]
            if not hasattr(v,'mask'):
                v = np.ma.array(v)
                v.mask = np.zeros(v.shape).astype('bool')
            self.__dict__[sds].append(v) 


class LEVELBCS(PACE):
    """ Read in PACE LevelB data for the cloud simulator"""
    def __init__(self,Path, SDS,verb=False):
        # Read one granule
        # ---------------------------------------------
        if type(Path) is list:
            if len(Path) == 0:
               raise Exception("WARNING: Empty PACE_LEVELB path")
        else:
           Path = [Path, ]

        # Create empty lists for SDS to be read from file
        # -----------------------------------------------
        self.granules = Path
        self.verb = verb
        self.SDS  = SDS

        for name in self.SDS:
           self.__dict__[name] = []

        self.scanStart = []

        # Read each granule, appending them to the list
        # ---------------------------------------------
        self._readList(Path)

        # Alias
        for sds in self.SDS:
            if sds in ALIAS:
                self.__dict__[ALIAS[sds]] = self.__dict__[sds]


        # Create corresponding python time
        # --------------------------------
        if hasattr(self,'midTime'):
            nscan,npixel = self.lon[0].shape
            self.tyme = []        
            for scanStart,midTime,lon in zip(self.scanStart,self.midTime,self.lon):
                scanStart  = isoparser(scanStart.strftime('2006-%m-%dT00:00:00'))
                tyme       = np.array([scanStart + timedelta(seconds=int(t)) for t in midTime])    
                tyme.shape = (nscan,1)
                tyme       = np.repeat(tyme,npixel,axis=1)
                tyme       = np.ma.array(tyme)
                tyme.mask  = lon.mask
                self.tyme.append(tyme)
            self.tyme = np.ma.concatenate(self.tyme)

        # convert lists to arrays        
        for sds in self.SDS:
            if sds in ALIAS:
                sds = ALIAS[sds]
            print 'Concatenating ',sds
            self.__dict__[sds] = np.ma.concatenate(self.__dict__[sds])
            if len(self.__dict__[sds].mask.shape) == 0:
                self.__dict__[sds].mask = np.zeros(self.__dict__[sds].shape).astype(bool)


#---
    def _readGranule(self,filename):
        """Reads one PACE LevelB granule."""

        # Exit if the file cannot be opened
        # ---------------------------------------
        try:
            if self.verb:
                print "[] Working on "+filename
            nc = Dataset(filename)
        except:
            raise Exception("- %s: not recognized as an netCDF file"%filename)

        # Read select variables (do not reshape)
        # --------------------------------------
        if len(self.scanStart) == 0:
            self.scanStart.append(isoparser(nc.time_coverage_start))
        for sds in self.SDS:            
            if len(self.__dict__[sds]) == 0:
                # Don't fuss if you can't find it
                try:
                    v = np.squeeze(nc.variables[sds][:])
                    if not hasattr(v,'mask'):
                        v = np.ma.array(v)
                        v.mask = np.zeros(v.shape).astype('bool')
                    self.__dict__[sds].append(v) 
                except:
                    pass
            
# ---
def granules ( path, t1, t2):
    """
    Returns a list of PACE granules for a given product at given time.
    On input,

    path      ---  mounting point for the PACE files
    t1        ---  starting time (timedate format)
    t2        ---  ending time (timedate format)
    """

    # Find PACE granules in the time range
    # ------------------------------------------
    dt = timedelta(minutes=5)
    t = datetime(t1.year,t1.month,t1.day,t1.hour,0,0)
    Granules = []
    while t <= t2:
        if t >= t1:
            doy = t.timetuple()[7]
            basen = "%s/Y%04d/M%02d/D%02d/OCI%04d%03d%02d%02d00.L1B_PACE.nc"\
                     %(path,t.year,t.month,t.day,t.year,doy,t.hour,t.minute)
            
            try:
                filen = glob(basen)[0]
                Granules += [filen,]
#               print " [x] Found "+filen
            except:
                pass
        t += dt

    if len(Granules) == 0:
        print "WARNING: no PACE granules found for %s through %s"%(str(t1), str(t2))

    return Granules

# ---
def granulesLB ( path, t1, t2, coll):
    """
    Returns a list of PACE LevelB granules for a given product at given time.
    On input,

    path      ---  mounting point for the PACE files
    t1        ---  starting time (timedate format)
    t2        ---  ending time (timedate format)
    coll      ---  GEOS data collection
    """

    # Find PACE granules in the time range
    # ------------------------------------------
    if type(coll) is str:
        coll = [coll]
    dt = timedelta(minutes=5)
    Granules = []
    for collname in coll:
        t = datetime(t1.year,t1.month,t1.day,t1.hour,0,0)
        while t <= t2:
            if t >= t1:
                basen = "%s/Y%04d/M%02d/D%02d/pace-g5nr.lb.%s.%04d%02d%02d_%02d%02d00.nc4"\
                         %(path,t.year,t.month,t.day,collname,t.year,t.month,t.day,t.hour,t.minute)
                
                try:
                    filen = glob(basen)[0]
                    Granules += [filen,]
    #               print " [x] Found "+filen
                except:
                    pass
            t += dt

    if len(Granules) == 0:
        print "WARNING: no PACE LevelB granules found for %s through %s"%(str(t1), str(t2))

    return Granules

# ---
def granulesLBN ( path, t1, t2, coll):
    """
    Returns a list of PACE LevelB granules for a given product at given time.
    On input,

    path      ---  mounting point for the PACE files
    t1        ---  starting time (timedate format)
    t2        ---  ending time (timedate format)
    coll      ---  GEOS data collection
    """

    # Find PACE granules in the time range
    # ------------------------------------------
    if type(coll) is str:
        coll = [coll]
    dt = timedelta(minutes=5)    
    Granules = []
    for collname in coll:
        t = datetime(t1.year,t1.month,t1.day,t1.hour,0,0)
        while t <= t2:
            if t >= t1:
                basen = "%s/Y%04d/M%02d/D%02d/pace-g5nr.lb-nearest.%s.%04d%02d%02d_%02d%02d00.nc4"\
                         %(path,t.year,t.month,t.day,collname,t.year,t.month,t.day,t.hour,t.minute)
                
                try:
                    filen = glob(basen)[0]
                    Granules += [filen,]
    #               print " [x] Found "+filen
                except:
                    pass
            t += dt

    if len(Granules) == 0:
        print "WARNING: no PACE LevelB granules found for %s through %s"%(str(t1), str(t2))

    return Granules
    
#------------------------------------ M A I N ------------------------------------

if __name__ == "__main__":
    

    # PACE default
    # -------------------
    calculon = '/nobackup/3/pcastell/PACE/L1B'
    nccs = '/discover/nobackup/projects/gmao/osse2/pub/c1440_NR/OBS/PACE'
    if os.path.exists(nccs): 
        L1Root = nccs + '/L1B'
    elif os.path.exists(calculon): 
        L1Root = calculon + '/L1B'


    dt_minutes = 5 #  5 minutes for now
    iso_t1 = '2020-03-24T00:50:00'
    t1     = isoparser(iso_t1)
    t2     = t1 + timedelta(minutes=options.dt_minutes)


    # Get Granules to work on
    # -----------------------
    Granules = granules (L1Root, t1, t2)

    # Create (x,y,t) coordinates
    # --------------------------
    pace = PACE(Granules)

