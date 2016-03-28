"""
Reads GOCI Granules

This software is hereby placed in the public domain.
Arlindo.daSilva@nasa.gov

"""

import os
import sys
from types    import *

from numpy    import zeros, ones, concatenate, array, shape, arange, tile
from numpy    import flipud, isnan, append, empty

from datetime import date, datetime, timedelta
from glob     import glob

from pyhdf.SD import SD, HDF4Error
#---  

SDS = ['Longitude', 
         'Latitude', 
         'Observation_time_day',
         'Observation_time_hour_utc',
         'Observation_time_minute',
         'Observation_time_month',
         'Observation_time_year',
         'QA_AOD_550nm',
         'AOD_550nm']

ALIAS = dict ( Longitude = 'lon', 
               Latitude = 'lat', 
               AOD_550nm = 'aod550')

MISSING = 999.999

#...........................................................................
class GOCIHandle(object):
    """
    Generic container for GOCI.
    """
    def __init__(self,name):
        self.name = name

class GOCI(object):
    """
    This class implements the GOCI interface.
    """

    def __init__ (self,Path,Verb=0,only_good=True):
        """
        Reads individual 1-hour granules or a full day of GOCI files
        present on a given *Path* and returns a single object with
        all data concatenated for a given algorithm. On input, 

        Required parameters:
        Path -- can be a single file, a single directory, of a list
                  of files and directories.  Directories are
                  transversed recursively. If a non GOCI hdf 
                  file is encountered, it is simply ignored.
        Optional parameters:
        Verb      -- Verbose level:
                  0 - really quiet (default)
                  1 - Warns if invalid file is found
                  2 - Prints out non-zero number of aerosols in each file.

        """

        #      Initially are lists of numpy arrays for each granule
        #      ------------------------------------------------
        self.verb = Verb
        self.sat  = None    # Satellite name
        self.col  = None    # collection, e.g., 005
        self.sample = None  # will hold sampled variables later
        self.ica = None     # will hold ICA indices later
        self.SDS = list(SDS)      # name of data streams
        self.path = Path
        self.tyme = []      # Python time 

        # Create empty lists for SDS to be read from file
        # -----------------------------------------------
        for sds in self.SDS:
            sds_ = sds.replace(' ','_')
            self.__dict__[sds_] = []

        # Read each granule, appending them to the list
        # ---------------------------------------------
        if type(Path) is ListType:
            if len(Path) == 0:
                self.nobs = 0
                print "WARNING: Empty GOCI object created"
                return
        else:
            Path = [Path, ]
        self._readList(Path)

        # Make some variables for convenience
        # ----------------------------
        year   = self.Observation_time_year 
        month  = self.Observation_time_month
        day    = self.Observation_time_day 
        hour   = self.Observation_time_hour_utc
        minute = self.Observation_time_minute

        # Create python times
        # ----------------
        print '[] Create Python Times'
        for y,m,d,h,t in zip(year,month,day,hour,minute):
            DATE_START = datetime(y,m,d,h,0,0)
            Tyme = None
            for minute in t:
                tyme  = array((DATE_START,)*len(minute))
                tyme[isnan(minute)] = float('NaN')
                mins  = array([timedelta(seconds=int(60*mm)) for mm in minute[~isnan(minute)]])
                tyme[~isnan(minute)] = tyme[~isnan(minute)] + mins
                tyme.shape = (1,len(tyme))
                if Tyme is None:
                    Tyme = tyme
                else:
                    Tyme = append(Tyme,tyme,axis=0)

            self.tyme.append(Tyme)

        self.SDS.append('tyme')

        # Mean Granule Time
        # ---------------------
        NaN = isnan(concatenate(self.Observation_time_minute))
        dt = concatenate(self.tyme)[~NaN].max()-concatenate(self.tyme)[~NaN].min()
        self.gtime = dt/2 + concatenate(self.tyme)[~NaN].min()

        # Make each attribute a single numpy array
        # ----------------------------------------        
        for sds in self.SDS:
            sds_ = sds.replace(' ','_')
            try:
                self.__dict__[sds_] = concatenate(self.__dict__[sds_])
            except:
                print "Failed concatenating "+sds


        if only_good:
            # Strip NANs and Low QA
            # Returns a 1-D array with only good obs
            # ---------------------
            print '[] Filter GOCI for only good'
            self._noNANs()
            self._qaFilter()


        # Aliases for convenience
        # ----------------------------
        Alias = ALIAS.keys()
        for sds in self.SDS:
          sds_ = sds.replace(' ','_')
          if sds_ in Alias:
                self.__dict__[ALIAS[sds_]] = self.__dict__[sds_] 
#---
    def _qaFilter(self):
        data_sds = []
        for sds in self.SDS:

            if 'minute' in sds:
                data_sds.append(sds)
            elif 'Observation_time' in sds:
                pass
            elif 'QA_AOD_550nm' in sds:
                pass
            else:
                data_sds.append(sds)

        # Loop through datasets, remove QA<3
        for sds in data_sds:
            sds_ = sds.replace(' ','_')
            self.__dict__[sds_] = self.__dict__[sds_][self.QA_AOD_550nm == 3]

        self.QA_AOD_550nm = self.QA_AOD_550nm[self.QA_AOD_550nm == 3]

#---
    def _noNANs(self):
        data_sds = []
        for sds in self.SDS:
            if 'Observation_time' not in sds or 'minute' in sds:
                data_sds.append(sds)
                
        # Loop through datasets, look for NANs
        # Remove NANs from the other datasets and from itself
        for sds in data_sds:
            if sds != 'tyme':
                sds_ = sds.replace(' ','_')
                d    = self.__dict__[sds_]
                for ssds in data_sds:
                    if ssds != sds:
                        ssds_ = ssds.replace(' ','_')
                        self.__dict__[ssds_] = self.__dict__[ssds_][~isnan(d)]

                self.__dict__[sds_] = d[~isnan(d)]


#---
    def linearSampleFile(self, inFile, onlyVars=None):
        """
        Interpolates some or all variables of inFile.
        """
        from gfio import GFIOctl

        # Open file
        # ---------
        f = GFIOctl(inFile)
        if self.sample == None:
            self.sample = GOCIHandle(inFile)
        if onlyVars is None:
            onlyVars = f.vname

        # Flatten coordinates
        # -------------------
        nt, nr = self.lon.shape
        #tyme = tile(self.tyme,(nr,1)).T
        lons = self.lon.ravel()
        lats = self.lat.ravel()
        tymes = self.tyme.ravel()
        #tymes[:] = self.gtime # use mean granule time

        tymes[isnan(self.Observation_time_minute.ravel())] = self.gtime
                 
        # Loop over variables on file
        # ---------------------------
        for v in onlyVars:
            if self.verb:
                print " <> Linear sampling ", v
            var = f.sample(v,lons,lats,tymes,Verbose=self.verb)
            if len(var.shape) == 1:
                self.sample.__dict__[v] = var.reshape((nt,nr))
            elif len(var.shape) == 2:
                var = var.T # shape should be (nobs,nz)
                self.sample.__dict__[v] = var.reshape((nt,nr,-1))
            else:
                raise IndexError, 'variable <%s> has rannk = %d'%len(var.shape)

#---
    def nearestSampleFile(self, inFile, onlyVars=None):
        """
        Sample some or all variables of inFile. This simply returns the
        values of the gridbox the pixel fall into (nearest neighbor).
        """
        from gfio import GFIOctl

        # Instiate grads and open file
        # ----------------------------
        f = GFIOctl(inFile)
        if self.sample == None:
            self.sample = GOCIHandle(inFile)
        if onlyVars is None:
            onlyVars = f.vname

        # Flatten coordinates
        # -------------------
        nt, nr = self.lon.shape
        #tyme = tile(self.tyme,(nr,1)).T
        lons = self.lon.ravel()
        lats = self.lat.ravel()
        tymes = self.tyme.ravel()
        #tymes[:] = self.gtime # use mean granule time

        tymes[isnan(self.Observation_time_minute.ravel())] = self.gtime
                 
        # Loop over variables on file
        # ---------------------------
        for v in onlyVars:
            if self.verb:
                print " <> Nearest sampling ", v
            var = f.sample(v,lons,lats,tymes,
                                algorithm='nearest',Verbose=self.verb)
            self.var = var
            if len(var.shape) == 1:
                self.sample.__dict__[v] = var.reshape((nt,nr))
            elif len(var.shape) == 2:
                var = var.T # shape should be (nobs,nz)
                self.sample.__dict__[v] = var.reshape((nt,nr,-1))
            else:
                raise IndexError, 'variable <%s> has rannk = %d'%len(var.shape)

#---
    def getICAindx(self,inFile):
        """
        Get gridded related indices for ICA calculations.  It calculates
        grid-box indices useful to implement Independent COlumn
        Approximation (ICA) type of algorithms.
        """
        from gfio import GFIOctl

        # Instiate grads and open file
        # ----------------------------
        f = GFIOctl(inFile)
        if self.ica == None:
            self.ica = GOCIHandle(inFile)

        # Handle grid indices for ICA algorithms
        # --------------------------------------
        nt, nr = self.lon.shape
        if self.verb:
            print " <> Performing ICA index generation"
        iCoord, jCoord = f.coordNN(self.lon.ravel(),self.lat.ravel())
        iS, jS = iCoord.astype('S5'), jCoord.astype('S5'), 
        keys = [ ii+','+jj for ii,jj in zip(iS,jS) ]
        Indices = dict()
        for n in range(iCoord.size):
            Indices[keys[n]] = (iCoord[n],jCoord[n])

        # save his for later
        # ------------------
        self.ica.iCoord = iCoord.reshape((nt,nr))
        self.ica.jCoord = jCoord.reshape((nt,nr))
        self.ica.Indices = Indices
          
#---
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
        """Reads one GOCI granule."""

        # Don't fuss if the file cannot be opened
        # ---------------------------------------
        try:
            if self.verb:
                print "[] Working on "+filename
            hfile = SD(filename)
        except HDF4Error:
            if self.verb > 2:
                print "- %s: not recognized as an HDF file"%filename
            return 

        # Read select variables (reshape to allow concatenation later)
        # ------------------------------------------------------------
        for sds in self.SDS:

            sds_ = sds.replace(' ','_')
            v = hfile.select(sds).get()
            rank = len(v.shape)
            if rank ==2:
                self.__dict__[sds_].append(flipud(v))
            else:
                self.__dict__[sds_].append(v) 

        #       Satellite name
        #       --------------
        if self.sat is None:
            self.sat = 'GOCI'

        #       Collection
        #       ----------
        if self.col is None:
            col = 1
            self.col = "%03d"%col

#............................................................................

if __name__ == "__main__":

    if os.path.exists('/discover/nobackup/pcastell/GOCI/'):
        gocifile = ['/discover/nobackup/pcastell/GOCI/20160316/GOCI_YAER_AOP_20160316001643.hdf']
    else:
        gocifile  = ['/nobackup/3/pcastell/GOCI/20160316/']

    g = GOCI(gocifile, Verb=1,only_good=False)

    #inFile = '/home/adasilva/opendap/fp/opendap/assim/inst1_2d_hwl_Nx'
    inFile  = '/discover/nobackup/pcastell/workspace/vis/GOCI/inst1_2d_hwl_Nx'
    #g.linearSampleFile(inFile,onlyVars=('TOTEXTTAU',))
