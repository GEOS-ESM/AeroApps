"""
Reads climate modeling grid 30 arc second MCD43D BRDF files
Keeping nearest neighbor values for lat, lon

"""

import os
import sys
from numpy    import loadtxt, array, tile, where, concatenate, flipud
from numpy    import ones, arange, meshgrid
from datetime import date, datetime, timedelta
from glob     import glob
from pyhdf.SD import SD, HDF4Error
from glob     import glob

MISSING = 32.767


SDS = dict (
      LAND = ('BRDF_Albedo_Parameter1_Band1','BRDF_Albedo_Parameter1_Band2',
              'BRDF_Albedo_Parameter1_Band3','BRDF_Albedo_Parameter1_Band4',
              'BRDF_Albedo_Parameter1_Band5','BRDF_Albedo_Parameter1_Band6',
              'BRDF_Albedo_Parameter1_Band7',
              'BRDF_Albedo_Parameter2_Band1','BRDF_Albedo_Parameter2_Band2',
              'BRDF_Albedo_Parameter2_Band3','BRDF_Albedo_Parameter2_Band4',
              'BRDF_Albedo_Parameter2_Band5','BRDF_Albedo_Parameter2_Band6',
              'BRDF_Albedo_Parameter2_Band7',
              'BRDF_Albedo_Parameter3_Band1','BRDF_Albedo_Parameter3_Band2',
              'BRDF_Albedo_Parameter3_Band3','BRDF_Albedo_Parameter3_Band4',
              'BRDF_Albedo_Parameter3_Band5','BRDF_Albedo_Parameter3_Band6',
              'BRDF_Albedo_Parameter3_Band7'),

      QUAL = ('BRDF_Albedo_Quality',
              'Snow_BRDF_Albedo',
              'BRDF_Albedo_Ancillary', )
           )

ALIAS = dict ( BRDF_Albedo_Parameter1_Band1 = 'KISO_b1_645',
               BRDF_Albedo_Parameter1_Band2 = 'KISO_b2_856',
               BRDF_Albedo_Parameter1_Band3 = 'KISO_b3_465',
               BRDF_Albedo_Parameter1_Band4 = 'KISO_b4_553',
               BRDF_Albedo_Parameter1_Band5 = 'KISO_b5_1241',
               BRDF_Albedo_Parameter1_Band6 = 'KISO_b6_1629',
               BRDF_Albedo_Parameter1_Band7 = 'KISO_b7_2114',
               BRDF_Albedo_Parameter2_Band1 = 'KVOL_b1_645',
               BRDF_Albedo_Parameter2_Band2 = 'KVOL_b2_856',
               BRDF_Albedo_Parameter2_Band3 = 'KVOL_b3_465',
               BRDF_Albedo_Parameter2_Band4 = 'KVOL_b4_553',
               BRDF_Albedo_Parameter2_Band5 = 'KVOL_b5_1241',
               BRDF_Albedo_Parameter2_Band6 = 'KVOL_b6_1629',
               BRDF_Albedo_Parameter2_Band7 = 'KVOL_b7_2114',
               BRDF_Albedo_Parameter3_Band1 = 'KGEO_b1_645',
               BRDF_Albedo_Parameter3_Band2 = 'KGEO_b2_856',
               BRDF_Albedo_Parameter3_Band3 = 'KGEO_b3_465',
               BRDF_Albedo_Parameter3_Band4 = 'KGEO_b4_553',
               BRDF_Albedo_Parameter3_Band5 = 'KGEO_b5_1241',
               BRDF_Albedo_Parameter3_Band6 = 'KGEO_b6_1629',
               BRDF_Albedo_Parameter3_Band7 = 'KGEO_b7_2114',               
             )

#...........................................................................

class McD43D(object):
    """
    This class implements the MODIS LAND BRDF daily Level 3 products, MCD43D (30 arc second horz res), 
    """

    def __init__ (self, Path, qPath, lon, lat, Verb=1):
        """
        Reads files for one day of Level 3 MCD43D
        present on a given *Path* and returns an object with
        all 3 kernels coeff. On input, 

        Required parameters:
         Path -- for now a single file.  Eventually implement a single directory, or a list
                 of files and directories. 
                            
        """
        if type(lon) is list:
            lon = array(lon)
            lat = array(lat)

        # List of HDF files for a given date
        #-----------------------------------
        self.verb = Verb
        self.SDS = SDS['LAND'] 
        if type(Path) is str:
            self.Files = glob(Path)
        else:
            self.Files = Path

        if type(qPath) is str:
            self.qFile = glob(qPath)[0]
        else:
            self.qFile = qPath[0]   

        # From a list of lat and lon, return the 
        # dx, dy on the grid
        # -------------------------------------  
        self.nobs = len(lon)
        self._findNearest(Path,lon,lat)

        # Read BRDF kernel in a MODIS tile
        # --------------------------------- 
        self.read_BRDFQuality()
        self.read_BRDF()

        # Result
                  
       
#--- 
    def _findNearest(self,path,lon,lat):

       """Given a list of lat, lon, return numbers to find the 
          position of the nearest neighbor on the grid (dx,dy)
       """
       dLon = 30.*1./3600
       dLat = 30.*1./3600
       Lon0 = -180 - dLon
       Lat0 = -90 + dLat

       self.dx = (0.5+(lon-Lon0)/dLon).astype(int)
       self.dy = (0.5+(lat-Lat0)/dLat).astype(int)

       if self.verb:
          print 'dx','dy', self.dx,self.dy
      
#---
    def read_BRDF(self):
        """Reads MCD43D file with Level 3 BRDF kernels for each MODIS band."""

        # Create empty lists for SDS to be read from file
        # -----------------------------------------------
        for name in self.SDS:
           self.__dict__[name] = [] 

        BRDF = MISSING * ones((len(self.SDS),self.nobs)) 
        for fn in self.Files:        
            print '<>Reading ',fn
            try:
                if self.verb:
                    print "[] Working on "+fn
                hfile = SD(fn)
            except HDF4Error:
                if self.verb > 2:
                    print "- %s: not recognized as an HDF file"%filename
                return 

            # Read select variables (reshape to allow concatenation later)
            # ------------------------------------------------------------
            for sds in self.SDS:  
                if self.verb:
                  print 'sds',self.SDS.index(sds) 
                try:                
                  v = hfile.select(sds).get()           
                  a = hfile.select(sds).attributes()
                  if a['scale_factor']!=1.0 or a['add_offset']!=0.0:
                      v = a['scale_factor'] * v + a['add_offset']
                  if self.verb:
                    print array(self.dx), BRDF.shape, BRDF[self.SDS.index(sds),:], v.shape 

                  v = flipud(v)
                  if self.qa[array(self.dy), array(self.dx)] < 1:
                    BRDF[self.SDS.index(sds),:] = v[array(self.dy), array(self.dx)]
                  else:
                    BRDF[self.SDS.index(sds),:] = MISSING

                except:
                  pass

        for sds in self.SDS:  
           self.__dict__[sds] = BRDF[self.SDS.index(sds),:]  
           if sds in ALIAS.keys():
               self.__dict__[ALIAS[sds]] = self.__dict__[sds] 
         
#---
    def read_BRDFQuality(self):
        """Reads MCD43D file with QA Flag for Level 3 BRDF kernels for all MODIS bands."""

        fn = self.qFile        
        print '<>Reading ',fn
        try:
            if self.verb:
                print "[] Working on "+fn
            hfile = SD(fn)
        except HDF4Error:
            if self.verb > 2:
                print "- %s: not recognized as an HDF file"%filename
            return 

        # Read QA variable
        # ----------------
        sds = 'BRDF_Quality' 
        v = hfile.select(sds).get()           

        self.qa = flipud(v)

#---

#............................................................................

if __name__ == "__main__":


    path = '/nobackup/3/pcastell/MODIS/MCD43D/Y2006/M12/MCD43D01.A2006365*'
    qpath = '/nobackup/3/pcastell/MODIS/MCD43D/Y2006/M12/MCD43D31.A2006365*'
      
    lon = array([-2.,-120.,15.2,17.2,170.1])
    lat = array([88.,40.,-20.,-20.,-55.5])

    # lon = arange(-180,180,1)
    # lat = arange(-90,90,1)
    # lon,lat = meshgrid(lon,lat)
    ex = McD43D(path,qpath,lon.flatten(),lat.flatten())

