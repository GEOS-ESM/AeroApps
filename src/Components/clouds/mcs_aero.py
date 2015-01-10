"""

MODIS Cloud Simulator ---- Aerosol Preparation.

Get species of aerosols on pixel / DISORT-level grid.
Optionally evaluate simplified extinction, SSA, and assymmetry parameter on this grid.
This produces a companion aerosol file to that produced my mcs.py.

"""

from numpy import savez, ones, zeros, array, size, sort, sum, \
                  linspace, any, log, mean, tile, amin, amax, exp, \
                  float32

from   pyobs.mxd03  import MxD03, MxD03Handle, DATE_START
from   pyobs        import NPZ


from mcsRegrid_     import getedgecoords, interppressure, \
                           regridtracer, regridmet, edget, edgeq

from MAPL.constants import *
from MAPL           import config

import MieObs_
import os

# variables needed for calculations
LINEAR_VARS  = { 'asm_Nv': ('PHIS','T', 'DELP', 'QV') }
NEAREST_VARS = { 'aer_Nv': ('DU001', 'DU002', 'DU003', 'DU004', 'DU005',
                            'SS001', 'SS002', 'SS003', 'SS004', 'SS005',
                            'SO4', 'BCPHOBIC', 'BCPHILIC',
                                   'OCPHOBIC', 'OCPHILIC') }

# DISORT level edges [km]
DISORT_LEVELS = array([80, 60, 50, 45, 40, 35, 30, 25, 20,
                       18, 16, 15, 14, 13, 12, 11, 10, 9,
                       8, 7, 6, 5, 4, 3, 2, 1, 0])

# wavelengths at which optical parameters are required [nm]
CHANNELS = array([
  0.41, 0.44, 0.47, 0.55, 0.65, 0.86, 0.91, 0.94, 1.24, 1.38, 1.63,
  2.13, 3.7, 3.9, 6.2, 7.3, 8.5, 11.0, 12.0, 13.2, 13.4, 13.8, 14.2
  ]) * 1.e3

class MCSError(Exception):
    def __init__(self, value):
        self.value = value
    def __str__(self):
        return repr(self.value)
  
class MCS(MxD03):
    """
    Implements the MODIS cloud simulator.
    """

    def getGEOS(self, rootdir, npzLinear=None, npzNearest=None):
        """
        Sample GEOS-5 variables on MODIS 1km swath.        
        """

        asm_Nv = rootdir + '/inst3_3d_asm_Nv'
        aer_Nv = rootdir + '/inst3_3d_aer_Nv'

        self.getICAindx(rootdir + '/inst3_2d_asm_Nx')

        # Linear interpolation
        # --------------------
        self.linearSampleFile(asm_Nv,onlyVars=LINEAR_VARS['asm_Nv'])
        self.sample.ZS = self.sample.PHIS / MAPL_GRAV

        self.linear = self.sample
        self.sample = None

        # Nearest neighbors - 3D
        # ----------------------
        self.nearestSampleFile(aer_Nv,onlyVars=NEAREST_VARS['aer_Nv'])

        self.nearest = self.sample
        self.sample = None

        # Hack: to avoid extrapolation issues, pretend the GEOS-5 top is
        #       highter then it actually is
        # --------------------------------------------------------------
        self._fixPTOP()

        # Save NPZ files, if desired (useful for debugging)
        # -------------------------------------------------
        if npzLinear is not None:
            if self.verb: print " - Saving Linear NPZ file ",npzLinear
            self.linear.PTOP = self.PTOP
            savez(npzLinear,**self.linear.__dict__)

        if npzNearest is not None:
            if self.verb: print " - Saving Nearest NPZ file ",npzNearest
            self.nearest.PTOP = self.PTOP
            savez(npzNearest,**self.nearest.__dict__)

    def getGEOS_readFromNPZ(self, rootdir, npzLinear, npzNearest):
        """
        Read previously sampled GEOS-5 variables (on MODIS 1km swath).        
        """

        # index pixels to GEOS grid
        self.getICAindx(rootdir + '/inst3_2d_asm_Nx')

        # read npz files from earlier run
        if self.verb: print " - Reading NPZ files"
        self.linear  = NPZ(npzLinear)
        self.nearest = NPZ(npzNearest)

        # validate and store globals
        if not (self.linear.PTOP == self.nearest.PTOP):
          raise ValueError, 'inconsistent <PTOP> between NPZ files'
        else:
          self.PTOP = self.linear.PTOP

#---
    def regridCoords(self):
        """
        Get DISORT and GEOS-5 vertical coordinates needed for vertical
        regridding.
        """

        g5 = self.linear
        di = MxD03Handle('DISORT')
        ns, nf, km = g5.DELP.shape

        # Get GEOS-Edge coordinates
        # -------------------------
        if self.verb: print " - Getting GEOS-5 edge coordinates"
        g5.PE, g5.ZE = getedgecoords(g5.T,g5.QV,g5.DELP,g5.ZS,self.PTOP)

        g5.logPE = None # compute this on demand later

        # Get DISORT edge pressure
        # ------------------------
        if self.verb: print " - Getting DISORT edge pressure"
        di.ZE = 1000 * DISORT_LEVELS  # in meters
        di.PE,di.PS,di.KS,rc = interppressure(di.ZE,self.ZS,g5.PE,g5.ZE)
        if rc: raise MCSError, \
          'Error on return from interpPressure(), rc = <%d>'%rc
        # *ks* is the surface level, 1-offset as in Fortran
        # ps is the surface pressure which has been adjusted for terrain height


        di.DELP = di.PE[:,:,1:]-di.PE[:,:,0:-1]

        di.logPE = None # compute this on demand later

        self.di = di

    def regridTracer(self,q):
         """
         Regrid a tracer from GEOS-5 to disort vertical coordinates.
         This uses a mass conserving scheme.
         """
         
         di = self.di
         g5 = self.linear
         di_q = regridtracer(di.PE,di.KS,q,g5.PE)

         return di_q

    def regridP(self,q):
         """
         Regrid *q* from GEOS-5 to disort vertical coordinates, using pressure
         as vertical coordinate.
         """
         di = self.di
         g5 = self.linear
         
         di_q = regridmet(di.PE,di.KS,q,g5.PE)

         return di_q

    def regridLogP(self,q):
         """
         Regrid *q* from GEOS-5 to disort vertical coordinates, using log-P
         as vertical coordinate.
         """
         di = self.di
         g5 = self.linear
         if g5.logPE == None:
             g5.logPE = log(g5.PE)
         if di.logPE == None:
             di.logPE = log(di.PE)
         
         di_q = regridmet(di.logPE,di.KS,q,g5.logPE)

         return di_q

#---
    def putOnDISORTgrid(self,linearNames=[],nearestNames=[]):
      '''
      Put 3D names chosen on the vertical DISORT grid
      '''
      self._DISORTregrid(linearNames,self.linear)
      self._DISORTregrid(nearestNames,self.nearest)

#---
    def _DISORTregrid(self,Names,data):
      ''' Put 3D names on vertical DISORT grid '''
      for var in Names:
        # never regrid linear.PE, which is needed itself for the regridding
        if data is self.linear and var == 'PE': continue
        if self.verb: print "[] Regridding %s"%var
        if len(data.__dict__[var].shape) == 3:
          data.__dict__[var] = self.regridTracer(data.__dict__[var])

#---
    def getOpticalProperties(self,RHfile,RHname='RH'):
      '''calculate simple optical properties: ATAU,ASSA,AASM'''
      from netCDF4 import Dataset
      # retrieve RH
      if self.verb: print "[] Retrieve Relative Humidity"
      nc = Dataset(RHfile,'r',format='NETCDF4')
      he = nc.variables['height_edge'][:]
      RH = nc.variables[RHname][:,:,:]
      nc.close()
      if RH.shape != self.di.DELP.shape:
        raise ValueError, 'RH on file of different dimensions'
      if any(abs(he-DISORT_LEVELS) > 1e-4):
        raise ValueError, 'RH file has different edge heights'
      # do optical calculation
      if self.verb: print "[] Calculating optical properties"
      sn = self.nearest
      sn.ATAU,sn.ASSA,sn.AASM = getAOPscalar(self.di.DELP,RH,sn,CHANNELS)

#---
    def writeNC ( self, filename, format='NETCDF4', zlib=False,
                  linearNames=[], nearestNames=[],
                  title='GEOS-5 MODIS Cloud Simulator Aerosols',
                  rcVars='variables.rc', Verbose=True):
        """
        Write a NetCDF file with GEOS-5 aerosols on the swath grid.
        """
        from netCDF4 import Dataset

        # Registry of variable long names and units
        # -----------------------------------------
        cf = config.Config(rcVars)

        # Open NC file
        # ------------
        nc = Dataset(filename,'w',format=format)

        # Set global attributes
        # ---------------------
        nc.title = title
        nc.institution = 'NASA'
        nc.source = 'Global Model and Assimilation Office'
        nc.history = 'Created from GEOS-5 standard collections by mcs.py'
        nc.references = 'n/a'
        nc.comment = 'Contains GEOS-5 aerosols nearest sampled on a MODIS granule.'
        nc.contact = 'Arlindo da Silva <arlindo.dasilva@nasa.gov>'
        nc.Conventions = 'CF'
        nc.mod03 = self.path

        # Dimension sizes
        # ---------------
        Ns, Nf = self.lon.shape
        Nk = len(DISORT_LEVELS)-1
        Nw = len(CHANNELS)

        # Create dimensions
        # -----------------
        ns = nc.createDimension('nscans',Ns)
        nf = nc.createDimension('nframes',Nf)
        nk = nc.createDimension('height',Nk)
        ne = nc.createDimension('height_edge',Nk+1)
        nw = nc.createDimension('wavelength',Nw)

        # Coordinate variables
        # --------------------
        nscans = nc.createVariable('nscans','f4',('nscans',))
        nscans.long_name = 'Fake GrADS Longitude'
        nscans.comment = 'For use in GrADS, use the 2D "longitude" for coordinates'
        nscans.missing_value = MAPL_UNDEF
        nscans.units = 'degrees_north'
        lon1, lon2 = min(self.lon[0,0], self.lon[-1,-1]), max(self.lon[0,0],self.lon[-1,-1])
        if lon2-lon1>100: # dateline crossing
            if lon1<0:
                tmp = lon2
                lon2 = 360.+lon1
                lon1 = tmp

        nframes = nc.createVariable('nframes','f4',('nframes',))
        nframes.long_name = 'Fake GrADS Latitude'
        nframes.comment = 'For use in GrADS, use the 2D "latitude" for coordinates'
        nframes.units = 'degrees_east'
        nframes.missing_value = MAPL_UNDEF
        lat1, lat2 = min(self.lat[0,0],self.lat[-1,-1]), max(self.lat[0,0], self.lat[-1,-1])

        # Fake coordinate variables
        # -------------------------
        nframes[:] = linspace(lon1,lon2,Nf)
        nscans[:] = linspace(lat1,lat2,Ns)

        height = nc.createVariable('height','f4',('height',))
        height.standard_name = 'height'
        height.long_name = 'DISORT Layer Height Above Sea-level'
        height.missing_value = MAPL_UNDEF
        height.units = 'km'
        height.positive = 'up'
        height.axis = 'z'
        height[:] = 0.5 * (DISORT_LEVELS[1:]+DISORT_LEVELS[0:-1])

        heighte = nc.createVariable('height_edge','f4',('height_edge',))
        heighte.standard_name = 'height_edge'
        heighte.long_name = 'DISORT Edge Height Above Sea-level'
        heighte.missing_value = MAPL_UNDEF
        heighte.units = 'km'
        heighte.positive = 'up'
        heighte.axis = 'z'
        heighte[:] = DISORT_LEVELS[:]
        
        wavelength = nc.createVariable('wavelength','f4',('wavelength',))
        wavelength.standard_name = 'wavelength'
        wavelength.long_name = 'MODIS channel wavelength'
        wavelength.missing_value = MAPL_UNDEF
        wavelength.units = 'nm'
        wavelength[:] = CHANNELS[:]

        longitude = nc.createVariable('longitude','f4',('nscans','nframes'))
        longitude.standard_name = 'longitude'
        longitude.long_name = 'Pixel Longitude'
        longitude.units = 'degrees_east'
        longitude.missing_value = MAPL_UNDEF
        longitude[:,:] = self.lon[:,:]

        latitude = nc.createVariable('latitude','f4',('nscans','nframes'))
        latitude.standard_name = 'latitude'
        latitude.long_name = 'Pixel Latitude'
        latitude.units = 'degrees_north'
        latitude.missing_value = MAPL_UNDEF
        latitude[:,:] = self.lat[:,:]

        time = nc.createVariable('time','i8',('nscans',))
        time.standard_name = 'time'
        time.long_name = 'EV Start Time'
        time.units = 'seconds since 1993-01-01 00:00:00'
        time.missing_value = int(MAPL_UNDEF)
        time[:] = _secondsSince(self.time,DATE_START)

        # Linear variables
        # ----------------
        self._writeVars(nc,cf,linearNames,self.linear,
                        comment='Linearly sampled variable',
                        Verbose=Verbose,zlib=zlib)
        
        # Nearest variables
        # -----------------
        self._writeVars(nc,cf,nearestNames,self.nearest,
                        comment='Nearest neighbor sampled variable',
                        Verbose=Verbose,zlib=zlib)

        # Close the file
        # --------------
        nc.close()

    def _writeVars(self,nc,cf,Names,data,comment=None,Verbose=True,zlib=True):
        """
        Write variables for a given group (linear or nearest).
        """
        for var in sort(Names):
            if Verbose:
                print "[] Writing %s"%var
            rank = len(data.__dict__[var].shape)
            long_name, units = cf(var).split(';')
            if rank == 2:
                v = nc.createVariable(var,'f4',('nscans','nframes',),
                                      fill_value=MAPL_UNDEF,zlib=zlib)
            elif rank == 3:
                v = nc.createVariable(var,'f4',('nscans','nframes','height'),
                                      fill_value=MAPL_UNDEF,zlib=zlib)
            elif rank == 4:
                v = nc.createVariable(var,'f4',('wavelength','nscans','nframes','height'),
                                      fill_value=MAPL_UNDEF,zlib=zlib)
            else:
                raise MCSError, 'Invalid rank for variable <%s>'%var
            v.long_name = long_name
            v.units = units.strip()
            v.missing_value = MAPL_UNDEF
            if comment != None:
                v.comment = comment
            if rank == 2:
                v[:,:] = data.__dict__[var][:,:]
            if rank == 3:
                v[:,:,:] = data.__dict__[var][:,:,:]
            else:
                v[:,:,:,:] = data.__dict__[var][:,:,:,:] 

    def _fixPTOP(self):
        """
        To avoid extrapolation issues, pretend the GEOS-5 top is
        highter then it actually is
        """
        pTop  = 0.01                      # actual GEOS-5 top pressure, in hPa
        pTop_ = 0.001                     # extended top pressure in hPa
        dTOP  = (pTop - pTop_) * 100.     # in Pa
        self.PTOP = pTop_ * 100           # in Pa
        self.linear.DELP[:,:,0] += dTOP
        try:
            self.nearest.DELP[:,:,0] += dTOP
        except:
            pass # it may not be there, no fuss

#........................................................................

def _secondsSince(tyme, t_ref):
  """
  Return array with elapsed seconds since t_ref.
  """
  dt = tyme - t_ref
  for n in range(size(dt)):
    dt[n] = dt[n].total_seconds()
  return dt
    
#........................................................................

def getAOPscalar(delp,rh,species,channels,
      rcfile='Aod_MODIS_Clouds.rc',
      Verbose=False):
    """
    Compute (tau,ssa,g) given aer object.
    Assumes delp, rh, and species mixing ratios are (ns,nf,nk)
    species has uppercase vnames (see below)
      that are mixing ratios of aerosol species
    outputs are (nch,ns,nf,nk)
    """
 
    # Make sure channels is a numpy array
    # -------------------------------------
    channels = _toArray(channels)
    nch = len(channels)
 
    # required inputs
    # ---------------
    vnames = [ 'du001', 'du002', 'du003', 'du004', 'du005',
               'ss001', 'ss002', 'ss003', 'ss004', 'ss005',
               'BCphobic', 'BCphilic',
               'OCphobic', 'OCphilic',
               'SO4' ]
    nq = len(vnames)

    # Pack inputs
    # -----------
    # 3D inputs are (ns,nf,nk) but want (nk,nobs)
    ns,nf,nk = delp.shape; nobs = ns*nf
    dm = delp.reshape((nobs,nk)).T / 9.81
    qm = ones((nk,nq,nobs),dtype=float32)
    rhx = rh.reshape((nobs,nk)).T
    for n, v in zip(range(nq),vnames):
      qm[:,n,:] = dm * species.__dict__[v.upper()].reshape((nobs,nk)).T
 
    # Do the Mie calculation
    # ----------------------
    tau,ssa,g,rc = MieObs_.getaopscalar(rcfile,channels,pad(vnames),Verbose,qm,rhx)
    if rc != 0:
      print "<<<ERROR>>> on return from MieObs_.getaopscalar, rc = ", rc
      raise ValueError, 'cannot get Aerosol Optical Properties (scalar version)'

    # outputs currently (nk,nch,nobs)
    # put into (nch,ns,nf,nk) order
    # ----------------------------
    tau = tau.transpose((1,2,0)).reshape((nch,ns,nf,nk))
    ssa = ssa.transpose((1,2,0)).reshape((nch,ns,nf,nk))
    g   =   g.transpose((1,2,0)).reshape((nch,ns,nf,nk))
 
    return (tau,ssa,g)

#--
def _toArray(channels):
    """Return numpy array, if *channels* is a scalar."""
    if type(channels) is int:
      channels = (channels, )
    elif type(channels) is float:
      channels = (channels, )
    channels = array(channels)
    return channels

#--
def pad(names):
    """
    Make all strings in list *names* the same size for f2py's benefit.
    """
    return [ "%-16s"%v for v in names ]


#.......................... D E G U G ......................................
def calculon():
    rootdir = '/nobackup/fp/opendap'
    mxd03   = '/nobackup/MODIS/Level1/MYD03/sample/MYD03.A2012228.1200.005.2012229145830.hdf'    
    outdir  = os.sep.join((os.environ['NOBACKUP'],'4',os.environ['USER'],'MCS','output'))
    return (rootdir,mxd03,outdir)

def nccs():
    rootdir = '/home/pnorris/iesa/aerosol/data/MCS/opendap/e5110_nwp/20120908'
    mxd03   = '/discover/nobackup/gwind/WMO_cases/MYD03.A2012252.1730.006.2012253144917.hdf'
    tag     =  os.path.splitext(os.path.basename(mxd03))[0]
    outdir  = os.sep.join((os.environ['NOBACKUP'],'MCS','output_'+tag))
    return (rootdir,mxd03,outdir)

if __name__ == "__main__":

    # setup
    rootdir, mxd03, outdir = nccs()
    species = True  # output species?
    optical = False # output optical properties?

    # instantiate MxD03 granule
    m = MCS(mxd03,Verb=1)

#   # associate GEOS fields with granule
#   # optionally write npz files for future bypass of this step
#   m.getGEOS(rootdir, 
#       npzLinear =os.sep.join((outdir,'geos-aer-linear.npz')),
#       npzNearest=os.sep.join((outdir,'geos-aer-nearest.npz')))

    # OR ... read sampled GEOS fields from old npz files
    m.getGEOS_readFromNPZ(rootdir,
        os.sep.join((outdir,'geos-aer-linear.npz')),
        os.sep.join((outdir,'geos-aer-nearest.npz')))

    # check all OK
    if m.linear == None and m.nearest == None:
      raise RuntimeError, 'must sample first using getGEOS()'

    # Prepare for vertical regridding
    m.regridCoords()

    # Put on disort grid
    m.putOnDISORTgrid(linearNames=[],
      nearestNames=NEAREST_VARS['aer_Nv'])

    # optionally calculate simple optical properties
    if optical:
      RHfile = os.sep.join((outdir,'geos-icacl-TOTWPDF-GCOP-SKEWT.nc4'))
      m.getOpticalProperties(RHfile,RHname='RHMCS')

    # decide on nearest outputs
    outputNames = []
    if species: outputNames.extend(NEAREST_VARS['aer_Nv'])
    if optical: outputNames.extend(('ATAU','ASSA','AASM'))

    # Write the output file
    ofile = 'geos-mcs-aero'
    if species: ofile += '-species'
    m.writeNC(
       os.sep.join((outdir,ofile+'.nc4')),
       nearestNames=outputNames,
       title = 'GEOS-5 MODIS Cloud Simulator Aerosol file')
