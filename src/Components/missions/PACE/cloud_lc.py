"""
Geostationary Cloud Simulator for VLIDORT.
"""

import numpy as np
from numpy import savez, ones, zeros, array, size, sort, sum, \
                  linspace, any, log, mean, tile, amin, amax, exp, squeeze

from   pyobs.gcs03  import GCS03, GCS03Handle
from   pyobs        import NPZ

from datetime import datetime, timedelta
from dateutil.relativedelta import relativedelta

from   ica          import genICA, getRe, getTau, getQsat

from MAPL.constants import *
from MAPL           import config

import os
from collections    import OrderedDict
import multiprocessing
from MAPL.ShaveMantissa_ import shave32

# Variable names for output file
LINEAR_NAMES = ('PHIS','PS','U10M','V10M','TS','DELP','AIRDENS','O3')
NEAREST_NAMES = ('FRLAND',)
ICA_NAMES = ('TAUL','TAUI','REL','REI','RH')

# dictionary of possible modes and their description
MODES = \
  {'HOMOCLD-MAXRAN':     'homogeneous clouds, maximum-random overlap',
   'HOMOCLD-COSP':       'homogeneous clouds, COSP overlap',
   'TOTWPDF-GCOP-SKEWT': 'total water PDF, Gaussian copula overlap, Skewed-Triangle PDF'
  }


def unwrap_CallgenICAClump(arg, **kwarg):
    return genICAfastClump(*arg, **kwarg)

def genICAfastClump(ncols,DELP,T,QV,QL,QI,CLOUD,PTOP,mode,lon,lat):
    QV_, QL_, QI_, rc = genICA(ncols,DELP,T,QV,QL,QI,CLOUD,PTOP,mode,lon,lat)    

    return QV_, QL_, QI_, rc

def unwrap_CallgenICA(arg, **kwarg):
    return genICAfast(*arg, **kwarg)

def genICAfast(ncols,DELP,T,QV,QL,QI,CLOUD,PTOP,mode):
    QV_, QL_, QI_, rc = genICA(ncols,DELP,T,QV,QL,QI,CLOUD,PTOP,mode)    

    return QV_, QL_, QI_, rc    

#---
def shave(q,undef=MAPL_UNDEF,has_undef=1,nbits=12):
    """
    Shave variable. On input, nbits is the number of mantissa bits to keep
    out of maximum of 24.
    """
    # Determine shaving parameters
    # ----------------------------
    xbits = 24 - nbits
    shp = q.shape
    rank = len(shp)
    if rank == 2:  # yx
        chunksize = shp[0]*shp[1]
    elif rank == 3: # zyx
        chunksize = shp[1]*shp[2]
    else:
        raise ValueError, "invalid rank=%d"%rank

    # Shave it
    # --------
    qs, rc = shave32(q.ravel(),xbits,has_undef,undef,chunksize)
    if rc:
        raise ValueError, "error on return from shave32, rc=%d"%rc

    return qs.reshape(shp)

class GCSError(Exception):
    def __init__(self, value):
        self.value = value
    def __str__(self):
        return repr(self.value)
  
class GCS(GCS03):
    """
    Implements the Geostationary Cloud Simulator for VLIDORT.
    """

    def getGEOS(self, rootdir, tyme, npzLinear=None, npzNearest=None):
        """
        Sample GEOS-5 variables on geostationary image.        
        """

        const_x   = rootdir + '/const_2d_asm_x'
        met1_x    = rootdir + '/inst_2d_met1_x'
        DELP_v    = rootdir + '/inst_3d_DELP_v'
        AIRDENS_v = rootdir + '/inst_3d_AIRDENS_v'
        O3_v      = rootdir + '/inst_3d_O3_v'
        T_v       = rootdir + '/inst_3d_T_v'
        U_v       = rootdir + '/inst_3d_U_v'
        QV_v      = rootdir + '/inst_3d_QV_v'
        QL_v      = rootdir + '/inst_3d_QL_v'
        QI_v      = rootdir + '/inst_3d_QI_v'
        CLOUD_v   = rootdir + '/inst_3d_CLOUD_v'

#T      cld_v = rootdir + '/tavg3_3d_cld_v'

        self.getICAindx(const_x)

        # Linear interpolation
        # --------------------
        self.linearSampleFile(  const_x,tyme,onlyVars=(   'PHIS',))
        self.linearSampleFile(   met1_x,tyme,onlyVars=(     'PS', 'TS', 'U10M', 'V10M',))
        self.linearSampleFile(   DELP_v,tyme,onlyVars=(   'DELP',))
        self.linearSampleFile(AIRDENS_v,tyme,onlyVars=('AIRDENS',))
        self.linearSampleFile(     O3_v,tyme,onlyVars=(     'O3',))
        self.linearSampleFile(      T_v,tyme,onlyVars=(      'T',))
        self.linearSampleFile(      U_v,tyme,onlyVars=(      'U',))

        self.linear = self.sample
        self.sample = None

        # Nearest neighbors
        # -----------------
        self.nearestSampleFile(const_x,tyme,onlyVars=('FRLAND',))
        self.nearestSampleFile( DELP_v,tyme,onlyVars=(  'DELP',))
        self.nearestSampleFile(    T_v,tyme,onlyVars=(     'T',))
        self.nearestSampleFile(   QV_v,tyme,onlyVars=(    'QV',))
        self.nearestSampleFile(   QL_v,tyme,onlyVars=(    'QL',))
        self.nearestSampleFile(   QI_v,tyme,onlyVars=(    'QI',))
        self.nearestSampleFile(CLOUD_v,tyme,onlyVars=( 'CLOUD',))

#T r    self.nearestSampleFile(cld_v,tyme,onlyVars=('QILS','QIAN'))
        print 'warning: QILS, QIAN not available in NatureRun output'			#T
        print '         Disabling anvil specific effective radii weighting'		#T
        # see reff() in mod_reff: setting QILS=QIAN=0 disables anvil Re weighting	#T
        self.sample.QILS = np.zeros_like(self.sample.QI)				#T
        self.sample.QIAN = np.zeros_like(self.sample.QI)				#T

        self.nearest = self.sample
        self.sample = None

        # Hack: to avoid extrapolation issues, pretend the GEOS-5 top is
        #       highter then it actually is
        # --------------------------------------------------------------
        self._fixPTOP()

        # Save NPZ files, if desired
        # --------------------------
        if npzLinear is not None:
            if self.verb: print "Saving Linear NPZ file ",npzLinear
            self.linear.PTOP = self.PTOP
            savez(npzLinear,**self.linear.__dict__)

        if npzNearest is not None:
            if self.verb: print "Saving Nearest NPZ file ",npzNearest
            self.nearest.PTOP = self.PTOP
            savez(npzNearest,**self.nearest.__dict__)

    def getGEOS_readFromNPZ(self, rootdir, npzLinear, npzNearest):
        """
        Read previously sampled GEOS-5 variables (on GCS image).        
        """

        # index pixels to GEOS grid
        const_x = rootdir + '/const_2d_asm_x'
        self.getICAindx(const_x)

        # read npz files from earlier run
        if self.verb:
            print " <> Retrieving previously saved interpolation files"
        self.linear  = NPZ(npzLinear)
        self.nearest = NPZ(npzNearest)

        # just to remind user again
        print 'warning: QILS, QIAN not available in NatureRun output'		#T
        print '         Disabling anvil specific effective radii weighting'	#T

        # validate and store globals
        if not (self.linear.PTOP == self.nearest.PTOP):
          raise ValueError, 'inconsistent <PTOP> between NPZ files'
        else:
          self.PTOP = self.linear.PTOP

#---
    def genICA(self, mode, clump=True, **kwopts):
        """
        Calculate TAU, RE, and RH using ICA generated subcolumns.
        """
        n = self.nearest
        l = self.linear
        ica = self.ica
        ica.Indices = OrderedDict(ica.Indices)

        # Space for optical properties
        # ----------------------------
        shp = n.T.shape
        for v in ( 'QV', 'QL', 'QI', 'RH', 'TAUL', 'TAUI' ):
            ica.__dict__[v] = MAPL_UNDEF * ones(shp)
        # REL, REI currently generated for all pixels below
        # if this changes, add above and use [ok] for them like TAUs
        
        DELP  = []
        T     = []
        QV    = []
        QL    = []
        QI    = []
        CLOUD = []
        ncols = []
        PTOP  = []
        modeList = []
        lon      = []
        lat      = []

        # subcolumn generation defaults to failure
        # ----------------------------------------
        ok = np.zeros(shp[0],dtype=bool)
 
        if self.verb:
            print " <> Performing ICA calculations"

        # Loop over unique gridboxes
        # --------------------------
        for gckey in ica.Indices:            
            ipxs = ica.Indices[gckey]
            # each pixel gets a generated subcolumn
            ncols.append(len(ipxs))

            # Select profile
            # --------------
            # can use GEOS profile for any pixel for "nearest" quantities
            # here we just use the first pixel
            ip = ipxs[0]
            DELP.append(n.DELP[ip,:])
            T.append(n.T[ip,:])
            QV.append(n.QV[ip,:])
            QL.append(n.QL[ip,:])
            QI.append(n.QI[ip,:])
            CLOUD.append(n.CLOUD[ip,:])
            PTOP.append(self.PTOP)
            modeList.append(mode)
            lon.append(self.lon[ipxs])
            lat.append(self.lat[ipxs])


        #pool = multiprocessing.Pool(int(multiprocessing.cpu_count()*0.5)) 
        pool = multiprocessing.Pool(4) 

        if clump:
            args = zip(ncols,DELP,T,QV,QL,QI,CLOUD,PTOP,modeList,lon,lat)
            result = pool.map(unwrap_CallgenICAClump,args)
        else:
            args = zip(ncols,DELP,T,QV,QL,QI,CLOUD,PTOP,modeList)
            result = pool.map(unwrap_CallgenICA,args)

        pool.close()
        pool.join()


        # load subcolumns if success and generate relative humidities
        # -----------------------------------------------------------
        for ii,gckey in enumerate(ica.Indices): 
            ipxs = ica.Indices[gckey]

            QV_, QL_, QI_, rc = result[ii]
            if not rc:
              ok[ipxs] = True
              ica.QV[ipxs], ica.QL[ipxs], ica.QI[ipxs] = QV_, QL_, QI_
              ica.RH[ipxs] = QV_ / tile(getQsat(self.PTOP,DELP[ii],T[ii]),(ncols[ii],1))

        # warning
        # -------
        if any(~ok): print 'warning: ICA procedure failed to generate some subcolumns'

        # Effective radius -- use linear interp fields
        # --------------------------------------------
        # currently, these are valid for all columns even if some had generation failure
        ica.REL, ica.REI = getRe(l.DELP,l.T,l.U,n.QILS,n.QIAN,self.PTOP)

        # Tau -- use ICA sampled QL/QI with linear interp RE
        # --------------------------------------------------
        ica.TAUL[ok],ica.TAUI[ok] = getTau(
          l.DELP[ok],ica.REL[ok],ica.REI[ok],ica.QL[ok],ica.QI[ok])

#---
    def doICA(self,filename,mode,clump=True):
        """
        Read and write variables doing ICA sampling.
        """

        if self.linear == None:
          raise RuntimeError, 'must sample first using getGEOS()'

        # Calulate TAU and RE
        # -------------------
        self.genICA(mode,clump=clump)
        
        # Write the output file
        # ---------------------
        self.writeNC(filename, 
         title = 'GEOS-5 Geostationary Cloud Simulator for VLIDORT' 
           + ' with ICA sampling (%s)' % MODES[mode])

#---
    def writeNC ( self, filename, format='NETCDF4', zlib=True,
                  linearNames=LINEAR_NAMES, nearestNames=NEAREST_NAMES, icaNames=ICA_NAMES,
                  title='GEOS-5 Geostationary Cloud Simulator for VLIDORT',
                  rcVars='variablesGCS.rc', Verbose=True):
        """
        Write a NetCDF file with simulated GEOS-5 variables on the swath grid.
        """
        from netCDF4 import Dataset

        # Registry of variable long names and units
        # -----------------------------------------
        cf = config.Config(rcVars)

        # Open NC file
        # ------------
        hf = Dataset(filename,'w',format=format)

        # Set global attributes
        # ---------------------
        hf.title = title
        hf.institution = 'NASA'
        hf.source = 'Global Model and Assimilation Office'
        hf.history = 'Created from GEOS-5 standard collections by gcs_vlidort.py'
        hf.references = 'n/a'
        hf.comment = 'This file contains GEOS-5 cloud related parameters sampled on a Geostationary image.'
        hf.contact = 'Arlindo da Silva <arlindo.dasilva@nasa.gov>'
        hf.Conventions = 'CF'
        hf.fgeom = self.fgeom

        # Dimension sizes
        # ---------------
        NS, EW = self.orgshape
        NK = self.nearest.T.shape[1]

        # Create dimensions
        # -----------------
        ns = hf.createDimension( 'ns',NS)
        ew = hf.createDimension( 'ew',EW)
        nk = hf.createDimension('lev',NK)

        # Coordinate variables
        # --------------------

        # fake horizontal coords
        lons = hf.createVariable('ew','f4',('ew',),fill_value=MAPL_UNDEF)
        lons.long_name = 'Fake GrADS Longitude'
        lons.comment = 'For use in GrADS, use the 2D "longitude" for coordinates'
        lons.units = 'degrees_east'
        lons.missing_value = MAPL_UNDEF
        lon1, lon2 = min(self.lon[0], self.lon[-1]), max(self.lon[0],self.lon[-1])
        if lon2-lon1>100: # dateline crossing
            if lon1<0:
                tmp = lon2
                lon2 = 360.+lon1
                lon1 = tmp
        lons[:] = linspace(lon1,lon2,EW)

        lats = hf.createVariable('ns','f4',('ns',),fill_value=MAPL_UNDEF)
        lats.long_name = 'Fake GrADS Latitude'
        lats.comment = 'For use in GrADS, use the 2D "latitude" for coordinates'
        lats.missing_value = MAPL_UNDEF
        lats.units = 'degrees_north'
        lat1, lat2 = min(self.lat[0],self.lat[-1]), max(self.lat[0], self.lat[-1])
        lats[:] = linspace(lat1,lat2,NS)

        # actual horizontal coords
        longitude = hf.createVariable('longitude','f4',('ns','ew'),fill_value=MAPL_UNDEF)
        longitude.standard_name = 'longitude'
        longitude.long_name = 'Pixel Longitude'
        longitude.units = 'degrees_east'
        longitude.missing_value = MAPL_UNDEF
        tmp = np.zeros(self.orgshape)
        tmp[~self.offview] = self.lon
        tmp[ self.offview] = MAPL_UNDEF
        longitude[:,:] = tmp

        latitude = hf.createVariable('latitude','f4',('ns','ew'),fill_value=MAPL_UNDEF)
        latitude.standard_name = 'latitude'
        latitude.long_name = 'Pixel Latitude'
        latitude.units = 'degrees_north'
        latitude.missing_value = MAPL_UNDEF
        tmp = np.zeros(self.orgshape)
        tmp[~self.offview] = self.lat
        tmp[ self.offview] = MAPL_UNDEF
        latitude[:,:] = tmp

        # eta levels
        lev = hf.createVariable('lev','f4',('lev',),fill_value=MAPL_UNDEF)
        lev.standard_name = 'model_layers'
        lev.long_name = 'vertical level'
        lev.missing_value = MAPL_UNDEF
        lev.units = 'layer'
        lev.positive = 'down'
        lev.coordinate = 'eta'
        lev[:] = np.arange(1,NK+1)

        # actual time in seconds past hour
        time = hf.createVariable('scanTime','f4',('ns','ew'),fill_value=MAPL_UNDEF)
        time.standard_name = 'scanTime'
        time.long_name = 'scan time'
        time.units = 'seconds past start hour'
        time.missing_value = MAPL_UNDEF
        tmp = np.zeros(self.orgshape)
        tmp[~self.offview] = self.scanTime
        tmp[ self.offview] = MAPL_UNDEF
        time[:,:] = tmp

        # Linear variables
        # ----------------
        self._writeVars(hf,cf,linearNames,self.linear,
                        comment='Linearly sampled variable',
                        Verbose=Verbose,zlib=zlib)
        
        # Nearest variables
        # -----------------
        self._writeVars(hf,cf,nearestNames,self.nearest,
                        comment='Nearest neighbor sampled variable',
                        Verbose=Verbose,zlib=zlib)

        # ICA variables
        # -------------
        self._writeVars(hf,cf,icaNames,self.ica,
                        comment='ICA sampled variable',
                        Verbose=Verbose,zlib=zlib)

        # Close the file
        # --------------
        hf.close()

    def _writeVars(self,hf,cf,Names,data,comment=None,Verbose=True,zlib=True):
        """
        Write variables for a given group (linear, nearest or ica).
        """
        for var in sort(Names):
            if Verbose:
                print "[] Writing %s"%var
            rank = len(data.__dict__[var].shape)
            long_name, units = cf(var).split(';')
            if rank == 1:
                v = hf.createVariable(var,'f4',('ns','ew',),
                                      fill_value=MAPL_UNDEF,zlib=zlib)
            elif rank == 2:
                v = hf.createVariable(var,'f4',('lev','ns','ew'),
                                      fill_value=MAPL_UNDEF,zlib=zlib)
            else:
                raise GCSError, 'Invalid rank for variable <%s>'%var
            v.long_name = long_name
            v.units = units.strip()
            v.missing_value = MAPL_UNDEF
            if comment != None:
                v.comment = comment
            if rank == 1:
                tmp = np.zeros(self.orgshape)
                tmp[self.offview] = MAPL_UNDEF

                tmp[~self.offview] = data.__dict__[var]
                v[:] = shave(tmp)
            else:
                if var in ('REL','REI'):
                    t2d = data.__dict__[var] * 1e6 # -> [um]
                else:
                    t2d = data.__dict__[var]
        
                tmp = np.zeros(v.shape)                
                a   = np.where(self.offview)
                tmp[:,a[0],a[1]] = MAPL_UNDEF
                a   = np.where(~self.offview)
                for k in range(v.shape[0]):
                    tmp[k,a[0],a[1]] = t2d[:,k]

                v[:] = shave(tmp)

    def _fixPTOP(self):
        """
        To avoid extrapolation issues, pretend the GEOS-5 top is
        highter then it actually is
        """
        pTop  = 0.01                      # actual GEOS-5 top pressure, in hPa
        pTop_ = 0.001                     # extended top pressure in hPa
        dTOP  = (pTop - pTop_) * 100.     # in Pa
        self.PTOP = pTop_ * 100           # in Pa
        self.linear.DELP[:,0] += dTOP
        try:
            self.nearest.DELP[:,0] += dTOP
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
    
def albedo_dE_cons(tauc,mu0,g,f):
  """
  Delta-Eddington cloud albedo for conservative scattering
  albedo = albedo_dE_cons(tauc,mu0,g,f)
    tauc --- cloud optical depth
    mu0  --- cosine solar zenith angle, must be positive
    g    --- asymmetry parameter
    f    --- forward scattering fraction
  Suggest using g = 0.85 and f = g**2
  """
  if mu0 <= 0.: raise ValueError, 'require mu0 > 0'
  cmu0 = 1.5*mu0; efac = exp(-tauc*(1.-f)/mu0)
  return 1. - 2.*((1+cmu0)+(1-cmu0)*efac)/(4.+3.*(1.-g)*tauc)

#.......................... D E G U G ......................................

def showdict(d,name):
    print '<<<<<< dictionary', name, 'begin >>>>>>'
    for key in d:
      t = type(d[key])
      print '>>', key, t
      if t in (int, float, str, np.float32, np.float64):
        print '  scalar:', d[key]
      elif t in (tuple, list):
        print '  sequence:', d[key]
      elif t in (np.ndarray, np.ma.core.MaskedArray):
        v = d[key]
        print '  numpy.ndarray: shape', v.shape, v.dtype
        if v.ndim == 0:
          print '  scalar:', v
        elif v.ndim == 1:
          print '  vect: first 5 vals:', v[:min(5,v.size)]
          print '  min, max:', np.amin(v), np.amax(v)
        elif v.ndim == 2:
          print '  array: first column:', v[0]
          print '  min, max:', np.amin(v), np.amax(v)
      print
    print '<<<<<< dictionary', name, 'end >>>>>>'
    print
    print

#   # debug
#   showdict(g.__dict__, 'g')
#   showdict(g.linear .__dict__, 'g.linear')
#   showdict(g.nearest.__dict__, 'g.nearest')
#   showdict(g.ica    .__dict__, 'g.ica')


if __name__ == "__main__":

    # setup
    rootdir   = '/home/pcastell/NatureRun/data/GCS/opendap/native'
    starttime = datetime(2006,04,01,17,0,0)
    endtime   = datetime(2006,9,01,17,0,0)
    #dt        = timedelta(hours=1)
    dt        = relativedelta(months=1)
    layout    = '41'
    outroot   = '/discover/nobackup/projects/gmao/osse2/pub/c1440_NR/OBS/TEMPO/CLD_DATA/'
    #outroot   = '/discover/nobackup/pcastell/GCS/output-start-hour_20051231_18z_tempo.lg1.cld.invariant.410'

    groot     = '/discover/nobackup/projects/gmao/osse2/pub/c1440_NR/OBS/TEMPO/CLD_DATA/'

    #  END OF USER INPUT
    # --------------------
    ntiles = int(layout[0])*int(layout[1])
    #ntiles = 1
    tyme   = starttime
    while tyme <= endtime:
        outdir  = os.sep.join((outroot,'LevelB/Y'+tyme.strftime('%Y'),'M'+tyme.strftime('%m'),'D'+tyme.strftime('%d')))
        #outdir  = outroot
        if not os.path.exists(outdir): os.makedirs(outdir)

        #for ll in [2,3]:
        for ll in range(ntiles):
            laycode = '{}{}'.format(layout,ll)

            print '<> Working on ',tyme,laycode

            gcs03   = '{}/LevelG/invariant/tempo.lg1.invariant.{}.nc4'.format(outroot,laycode)

            # instantiate Geostationary granule
            g = GCS(gcs03, varTime=True, Verb=1)

            # associate GEOS fields with granule
            # optionally write npz files for future bypass of this step
            # OR ... read sampled GEOS fields from old npz files
            try:
                g.getGEOS_readFromNPZ(rootdir,
                  os.sep.join((outdir,'geos-linear.'+tyme.strftime('%Y%m%d_%H')+'z.'+laycode+'.npz')),
                  os.sep.join((outdir,'geos-nearest.'+tyme.strftime('%Y%m%d_%H')+'z.'+laycode+'.npz')))
            except:
                g.getGEOS(rootdir, tyme,
                    npzLinear =os.sep.join((outdir,'geos-linear.'+tyme.strftime('%Y%m%d_%H')+'z.'+laycode+'.npz')),
                    npzNearest=os.sep.join((outdir,'geos-nearest.'+tyme.strftime('%Y%m%d_%H')+'z.'+laycode+'.npz')))


            # mode of ICA? (from MODES)
        #   mode = 'HOMOCLD-MAXRAN' 
        #   mode = 'HOMOCLD-COSP' 
            mode = 'TOTWPDF-GCOP-SKEWT'

        #     # do the simulation and write out results
            ofile = os.sep.join((outdir,'tempo-g5nr-icacl-'+mode+'.'+tyme.strftime('%Y%m%d_%H')+'z.'+laycode+'.nc4'))
            print '<> creating ',ofile
            g.doICA(ofile,mode,clump=True)

            os.remove(os.sep.join((outdir,'geos-linear.'+tyme.strftime('%Y%m%d_%H')+'z.'+laycode+'.npz')))
            os.remove(os.sep.join((outdir,'geos-nearest.'+tyme.strftime('%Y%m%d_%H')+'z.'+laycode+'.npz')))

        tyme += dt
