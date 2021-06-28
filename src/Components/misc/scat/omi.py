"""
  OMI class extends ODS.
"""

import OMI_  # f2py extensions


import os
from pyods import ODS
from numpy import save, zeros, ones, linspace, any, arange, savez, load
#from pylab import *

from scat.scat_binObs_    import binobs2d, binobs3d
from gfio       import GFIO

ktPS = 33
ktAOD = 45
ktALBEDO = 47
ktAAOD = 49
ktSOLAR_ZENITH = 54
ktRELAT_AZYMUTH = 55
ktSENSOR_ZENITH = 56
ktATYPE = 60
ktREFLECTIVITY = 78
ktRADIANCE = 80
ktAI = 81

MISSING = 999.999


class OMI(object):

    def __init__ (self, filename,nymd=None,nhms=None,only_good=True):
        """
        Reads an ODS or NPZ file and creates an OMI object with the following
        attributes:

             nobs           ---  number of profiles
             nch            ---  number of channels
             lon            ---  longitudes in degrees        (nobs)
             lat            ---  latitudes in degrees         (nobs)
             channels       ---  channels in nm               (nch)

             ps             ---  surface pressure in Pascal   (nobs)
             solar_zenith   ---  solar zenith angle           (nobs)
             relat_azymuth  ---  relative azymuth angle          (nobs)
             sensor_zenith  ---  sensor zenith angle          (nobs)

             albedo         ---  surface albedo               (nobs,nch)
             reflectivity   ---  reflectivity                 (nobs,nch)
             radiance       ---  radiance                     (nobs,nch)
           
        """

#       Read in the data
#       ----------------
        ext = os.path.splitext(filename)[1] 
        if ext == '.ods':
            self._ReadODS (filename,nymd,nhms,only_good)
        elif ext == '.npz':
            self._ReadNPZ (filename)
        else:
            raise IOError, "Unknown extension %s, it must be either .ods or .npz"%ext 

#       Create space for simulations (ODS does not have it)
#       ---------------------------------------------------
        if not self.__dict__.__contains__('radiance_'):
           self.radiance_= MISSING * ones(self.radiance.shape)
        if not self.__dict__.__contains__('ai_'):
           self.ai_ = MISSING * ones(self.ai.shape)

#---
    def _ReadODS (self, filename=None,nymd=None,nhms=None,only_good=True):
        """
        Reads an ODS file and creates an OMI object.
        """
        
#       Load ODS file
#       -------------
        ods = ODS(filename,nymd,nhms,only_good)

#       Fundamental data types (not retrieved)
#       --------------------------------------
        self.ps            = ods.select(kt=ktPS)
        self.solar_zenith  = ods.select(kt=ktSOLAR_ZENITH)
        self.relat_azymuth = ods.select(kt=ktRELAT_AZYMUTH)
        self.sensor_zenith = ods.select(kt=ktSENSOR_ZENITH)
        self.albedo        = ods.select(kt=ktALBEDO)
        self.reflectivity  = ods.select(kt=ktREFLECTIVITY)
        self.radiance      = ods.select(kt=ktRADIANCE)
        self.ai            = ods.select(kt=ktAI,lev=354) # Notice 1 channel only

#       Consistency checks: use radiance and ps as reference; these
#       must be defined for all soundings
#       -----------------------------------------------------------
        nobs = self.ps.nobs
        nch  = self.radiance.nobs / self.ps.nobs 
        
        if self.radiance.nobs%self.ps.nobs != 0:
            raise ValueError, "radiance.size=% is not a multiple of ps.size=%d"%\
             (self.radiance.nobs,self.ps.nobs)
 
#       Pull out relevant arrays out of ODS object
#       ------------------------------------------
        self.nymd         = nymd
        self.nhms         = nhms   
        self.nobs         = nobs
        self.nch          = nch
        self.ks           = self.ps.ks
        self.lon          = self.ps.lon
        self.lat          = self.ps.lat
        self.channels     = self.radiance.lev.reshape(nobs,nch)[0,:]

#       Construct hash table with index of each ks
        KS = {}
        for k in range(nobs):
            KS[str(self.ks[k])] = k 

        self.ps            = 100 * self.ps.obs  # convert to Pascal

#       NOTE: in order to pull out retrieved quantities one needs to
#             to match the sounding index ks because not all locations
#             have a sucessfull retrieval. To do.

#       Get obs, adding missing if sounding is missing
#       ----------------------------------------------
        self.solar_zenith  = _getObs('solar_zenith',self.solar_zenith,KS)
        self.relat_azymuth = _getObs('relative_azimuth',self.relat_azymuth,KS)
        self.sensor_zenith = _getObs('sensor_zenith',self.sensor_zenith,KS)
        self.ai            = _getObs('AI',self.ai,KS)

        self.albedo        = _getObsNch('albedo',self.albedo,KS,self.channels)
        self.reflectivity  = _getObsNch('reflectivity',self.reflectivity,KS,self.channels)
        self.radiance      = _getObsNch('radiance',self.radiance,KS,self.channels)

#---
    def _ReadNPZ (self, filename ):
        """
        Reads an NPZ file (written by method writeu) and creates an OMI object.        
        """
        npz = load(filename)

        self.nymd, self.nhms, self.nobs, self.nch, version = npz['meta']

        self.lon = npz['lon']
        self.lat = npz['lat']
        self.ps = npz['ps']
        self.ks = npz['ks']
        self.channels = npz['channels']
        self.solar_zenith = npz['solar_zenith']
        self.relat_azymuth = npz['relat_azymuth']
        self.sensor_zenith = npz['sensor_zenith']
        self.albedo = npz['albedo']
        self.reflectivity = npz['reflectivity']
        self.radiance = npz['radiance']
        self.radiance_ = npz['radiance_']
        self.ai = npz['ai']
        self.ai_ = npz['ai_']

#--
    def vlidort(self, nMom, nPol, 
                tau, ssa, g,pmom, pe, ze, te,
                scalar=False, verbose=0, I=None):
        """
        Uses VLIDORT to compute radiances. On Input,

           nMom  --- number of moments
           nPol  --- number of non-zero elements of scattering phase matrix
           tau   --- aerosol optical depth
           ssa   --- aerosol single scattering albedo
           g     --- aerosol asymmetry factor
           pmon  --- aerosol scattering phase matrix
           pe    --- pressure at edges [Pa]
           ze    --- height at edges [m]
           te    --- temperature at edges [K]

        Optionally,

           scalar --- If True, perform a scalar calcular (default is vector)
           verbose -- Verbosity level
           I     --- index of subset of observations to process

        """

        nch = self.nch
        if I is None:
            I = range(self.lon.shape[0])

#       Scalar calculation
#       ------------------
        if scalar:
           radiance_, ai_, reflectivity_, rc = OMI_.scalar(self.nymd,self.nhms,self.channels,\
                                             tau,ssa,g, pe,ze,te, \
                                             self.lon[I],self.lat[I],self.ps[I],self.albedo[I],\
                                             self.solar_zenith[I],self.relat_azymuth[I],\
                                             self.sensor_zenith[I],MISSING,self.radiance[I], \
                                             self.ai[I], verbose)
#       Vector calculation
#       ------------------
        else:
           radiance_, ai_, reflectivity_, rc = OMI_.vector(self.nymd,self.nhms,self.channels,\
                                           tau, ssa, g, pmom, pe,ze,te, \
                                           self.lon[I],self.lat[I],self.ps[I],self.albedo[I],\
                                           self.solar_zenith[I],self.relat_azymuth[I],\
                                           self.sensor_zenith[I], MISSING,self.radiance[I], \
                                           self.ai[I],verbose)

        if rc != 0:
            raise ValueError, "on return from OMI_.scalar/vector, rc = "+str(rc)

#       Save results in object
#       ----------------------
        self.radiance_[I] = radiance_
        self.ai_[I] = ai_
        
#---
    def writeu(self,filename=None,dir='.',expid='vlidort',Verb=1):
        """
        Writes the un-gridded OMI object to a numpy npz file. 
        """

        if filename is None:
            filename = '%s/%s.omi.%d_%02dz.npz'%(dir,expid,self.nymd,self.nhms/10000)

#       Create simulated radiance_/ai_ if not there
#       -------------------------------------------
        if not self.__dict__.__contains__('radiance_'):
           self.radiance_= MISSING * ones(self.radiance.shape)
        if not self.__dict__.__contains__('ai_'):
           self.ai_ = MISSING * ones(self.ai.shape)

        version = 1 # File format version
        meta = [self.nymd,self.nhms,self.nobs,self.nch,version]
        savez(filename,
                            meta = meta,
                             lon = self.lon,
                             lat = self.lat,
                             ps  = self.ps,
                              ks = self.ks,
                        channels = self.channels,
                    solar_zenith = self.solar_zenith,
                   relat_azymuth = self.relat_azymuth,
                   sensor_zenith = self.sensor_zenith,
                          albedo = self.albedo,
                    reflectivity = self.reflectivity,
                        radiance = self.radiance,
                       radiance_ = self.radiance_,
                              ai = self.ai,
                             ai_ = self.ai_)

        if Verb >=1:
            print "[w] Wrote file "+filename

#---
    def writeg(self,filename=None,dir='.',expid='vlidort',refine=4,res=None,Verb=1):
       """
        Writes gridded OMI measurements and GEOS-5 VLIDORT simulations
        to file.

         refine  -- refinement level for a base 4x5 GEOS-5 grid
                       refine=1  produces a   4  x  5    grid
                       refine=2  produces a   2  x2.50   grid
                       refine=4  produces a   1  x1,25   grid
                       refine=8  produces a  0.50x0.625  grid
                       refine=16 produces a  0.25x0.3125 grid
        Alternatively, one can specify the grid resolution with a
        single letter:

         res     -- single letter denoting GEOS-5 resolution,
                       res='a'  produces a   4  x  5    grid
                       res='b'  produces a   2  x2.50   grid
                       res='c'  produces a   1  x1,25   grid
                       res='d'  produces a  0.50x0.625  grid
                       res='e'  produces a  0.25x0.3125 grid

                   NOTE: *res*, if specified, supersedes *refine*.

         Verb -- Verbose level:
                 0 - really quiet (default)
                 1 - Warns if invalid file is found
                 2 - Prints out non-zero number of fires in each file.


       """

#      Output grid resolution
#      ----------------------
       if res is not None:
           if res=='a': refine = 1 
           if res=='b': refine = 2
           if res=='c': refine = 4
           if res=='d': refine = 8
           if res=='e': refine = 16

#      Lat lon grid
#      ------------
       dx = 5. / refine
       dy = 4. / refine
       im = int(360. / dx)
       jm = int(180. / dy + 1)

       glon = linspace(-180.,180.,im,endpoint=False)
       glat = linspace(-90.,90.,jm)

       nch = self.nch
       nymd = self.nymd
       nhms = self.nhms

       vtitle = [ 'Albedo',
                  'Aerosol Index (Retrieved by OMI)',
                  'Aerosol Index (Simulated by GEOS-5)',
                  'Radiance (Measured by OMI)',
                  'Radiance (Simulated by GEOS-5)' ]

       vname  = ['albedo','ai', 'ai_', 'rad', 'rad_' ]
       vunits = [ '%',    '1',  '1',   '1',   '1' ]
       kmvar  = [nch,      0,    0,    nch,   nch ]

       title = 'Gridded OMI Measurements and VLIDORT Simulations from GEOS-5'
       source = 'NASA/GSFC/GMAO GEOS-5 Aerosol Group'
       contact = 'arlindo.dasilva@nasa.gov'

       if filename is None:
           filename = '%s/%s.omi.%d_%02dz.nc4'%(dir,expid,self.nymd,self.nhms/10000)

#      Create the file
#      ---------------
       f = GFIO()
       f.create(filename, vname, nymd, nhms,
                lon=glon, lat=glat, levs=self.channels, levunits='nm',
                vtitle=vtitle, vunits=vunits,kmvar=kmvar,amiss=MISSING,
                title=title, source=source, contact=contact)

#      Grid variable and write to file
#      -------------------------------
       f.write('albedo', nymd, nhms, 
               binobs3d(self.lon,self.lat,self.albedo,im,jm,MISSING) )
       f.write('ai', nymd, nhms, 
               binobs2d(self.lon,self.lat,self.ai,im,jm,MISSING) )
       f.write('ai_', nymd, nhms, 
               binobs2d(self.lon,self.lat,self.ai_,im,jm,MISSING) )
       f.write('rad', nymd, nhms, 
               binobs3d(self.lon,self.lat,self.radiance,im,jm,MISSING) )
       f.write('rad_', nymd, nhms, 
               binobs3d(self.lon,self.lat,self.radiance_,im,jm,MISSING) )

       try:
           f.close()
       except:
           pass

       if Verb >=1:
           print "[w] Wrote file "+filename

#--

def _getObs(name,ods,KS):
    """
    Given an ODS structure, returns a 1D array with the "obs"
    attribute, with MISSING values used for the missing soundings.
    On input, *KS* is a hash with the index of each sounding index.
    """
    obs = MISSING * ones(len(KS))
    i = [ KS[str(s)] for s in ods.ks ] # index of each available sounding
    obs[i] = ods.obs
    return obs

def _getObsNch(name,ods,KS,channels):
    """
    Given an ODS structure, returns a 1D array with the "obs"
    attribute, with MISSING values used for the missing soundings.
    On input, *ns* is the number of sounds (number of pixels in the
    orbit), and *channels* are the channels.
    """
    ns  = len(KS)
    nch = channels.shape[0] 
    obs = MISSING * ones((ns,nch))
    j = 0
    for ch in channels:
        ods_ = ods.select(lev=ch)
        obs[:,j] = _getObs(name,ods_,KS)
        j = j + 1
    return obs
        
#.............................................................................................

if __name__ == "__main__":

    """Simple unit testing..."""

    from mie import getMieScal, getMieVect, getTopo
   
#   Read OMI data
#   -------------
    nymd = 20080629
    nhms = 0000
    omi = OMI('ods/orig/omi.aero_tc8.obs.20080629.ods',nymd,nhms)
#    omi = OMI('/nobackup/OMI/Level2/ODS/Y2008/M06/omi.aero_tc8.obs.20080611.ods',nymd,nhms)
  

#   Compute Mie parameters at the OMI locations
#   -------------------------------------------
    aerosols = 'data/a5_arctas_02.inst3d_aer_v.20080629_0000z.nc'
#    aerosols = '/nobackup/1/ARCTAS/Y2008/M06/d5_arctas_02.inst3d_aer_v.20080611_1500z.nc'

#   Vector with loop :
#   -------
    mobs = 100
    for i in range(0, omi.nobs, mobs):

             I = range(i,min(i+mobs,omi.nobs))
             nMom, nPol, tau, ssa, pe, ze, te, g, pmom = getMieVect(aerosols,omi.channels,\
             omi.lon[I],omi.lat[I],nymd,nhms)
            
             zs = getTopo('data/topography.1152x721.nc',omi.lon[I],omi.lat[I])
             he = ze + zs.reshape((1,)+zs.shape)
         
             omi.vlidort(nMom,nPol,tau,ssa,g,pmom,pe,he,te,verbose=2,scalar=False, I=I)


#   Write gridded output file
#   -------------------------
    omi.writeg(filename='omi_test_simulation_11_06_2008_0000_bis.nc')


