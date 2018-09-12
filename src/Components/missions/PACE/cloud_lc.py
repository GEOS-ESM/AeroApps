"""
Cloud Simulator for PACE
"""
from optparse        import OptionParser
import numpy as np
from numpy import savez, ones, zeros, array, size, sort, sum, \
                  linspace, any, log, mean, tile, amin, amax, exp, squeeze
from   dateutil.parser import parse  as isoparser

from   pace  import granulesLB, granulesLBN, LEVELBCS
from   pyobs.gcs03 import GCS03
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

# Variables Names to read in
lSDS = ['PHIS','PS','TS','U10M','V10M','DELP','AIRDENS','O3','T','U']
nSDS = ['FRLAND','DELP','T','QV','QL','QI','CLOUD']
mSDS = ['longitude', 'latitude', 'ev_mid_time']

class HOLDER(LEVELBCS):
    """
    Generic container.
    """
    def __init__(self):
        pass


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
  
class PCS(LEVELBCS,GCS03):
    """
    Implements the PACE Cloud Simulator for VLIDORT.
    """
    def __init__(self,linearFiles,nearestFiles,const_x):

        self.ica = None
        # Read in sampled data
        LEVELBCS.__init__(self,nearestFiles,mSDS)
        self.clon = self.lon
        self.clat = self.lat
        self.offview = self.clon.mask
        self.nobs = np.sum(~self.offview)
        #store flattened lon/lat
        self.lon = self.clon[~self.offview]
        self.lat = self.clat[~self.offview]
        scanTime       = self.midTime
        scanTime.shape += (1,)
        nscan,npixel = self.clon.shape   
        scanTime       = np.repeat(scanTime,npixel,axis=1)        
        self.scanTime = scanTime[~self.offview]



        self.linear = HOLDER()
        self.nearest = HOLDER()
        self.linear.offview = self.offview
        self.nearest.offview = self.offview
        self.linear.nobs = self.nobs
        self.nearest.nobs = self.nobs
        print 'Reading Linear...'
        LEVELBCS.__init__(self.linear,linearFiles,lSDS,)
        print 'Reading Nearest...'
        LEVELBCS.__init__(self.nearest,nearestFiles,nSDS)

        # # Reshape arrays to (nobs,nz)
        # for sds in lSDS:
        #     print 'reshaping lnear ', sds
        #     v = self.linear.__dict__[sds]
        #     if len(v.shape) == 2:
        #         self.linear.__dict__[sds] = v[~self.offview]
        #     else:
        #         self.linear.__dict__[sds] = np.zeros([v.shape[0],self.nobs])
        #         for k in range(v.shape[0]):
        #             self.linear.__dict__[sds][k,:] = v[k,:,:][~self.offview]

        #         self.linear.__dict__[sds] = self.linear.__dict__[sds].T

        # for sds in nSDS:
        #     print 'reshaping nearest ', sds
        #     v = self.nearest.__dict__[sds]
        #     if len(v.shape) == 2:
        #         self.nearest.__dict__[sds] = v[~self.offview]
        #     else:
        #         self.nearest.__dict__[sds] = np.zeros([v.shape[0],self.nobs])
        #         for k in range(v.shape[0]):
        #             self.nearest.__dict__[sds][k,:] = v[k,:,:][~self.offview]

        #         self.nearest.__dict__[sds] = self.nearest.__dict__[sds].T


        self.getICAindx(const_x)

        print 'warning: QILS, QIAN not available in NatureRun output'           
        print '         Disabling anvil specific effective radii weighting'     
        # see reff() in mod_reff: setting QILS=QIAN=0 disables anvil Re weighting   
        self.nearest.QILS = np.zeros_like(self.nearest.QI)                
        self.nearest.QIAN = np.zeros_like(self.nearest.QI)                

        # Hack: to avoid extrapolation issues, pretend the GEOS-5 top is
        #       highter then it actually is
        # --------------------------------------------------------------
        self._fixPTOP()

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

        if self.linear is None:
          raise RuntimeError, 'must sample first using getGEOS()'

        # Calulate TAU and RE
        # -------------------
        self.genICA(mode,clump=clump)
        
        # # Write the output file
        # # ---------------------
        # self.writeNC(filename, 
        # title = 'GEOS-5 PACE Cloud Simulator for VLIDORT' 
        #   + ' with ICA sampling (%s)' % MODES[mode])

#---
    def writeNC ( self, filename, format='NETCDF4', zlib=True,
                  linearNames=LINEAR_NAMES, nearestNames=NEAREST_NAMES, icaNames=ICA_NAMES,
                  title='GEOS-5 PACE Cloud Simulator for VLIDORT',
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
        hf.history = 'Created from GEOS-5 standard collections by cloud_lc.py'
        hf.references = 'n/a'
        hf.comment = 'This file contains GEOS-5 cloud related parameters sampled on a PACE granule.'
        hf.contact = 'Patricia Castellanos <patricia.castellanos@nasa.gov>'
        hf.Conventions = 'CF'


        # Dimension sizes
        # ---------------
        self.orgshape = self.clon.shape
        NS, EW = self.orgshape
        NK = self.nearest.T.shape[1]

        # Create dimensions
        # -----------------        
        ew = hf.createDimension( 'ccd_pixels',EW)
        ns = hf.createDimension( 'number_of_scans',NS)
        nk = hf.createDimension('lev',NK)

        # Coordinate variables
        # --------------------

        # fake horizontal coords
        lons = hf.createVariable('ccd_pixels','f4',('ccd_pixels',),fill_value=MAPL_UNDEF)
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

        lats = hf.createVariable('number_of_scans','f4',('number_of_scans',),fill_value=MAPL_UNDEF)
        lats.long_name = 'Fake GrADS Latitude'
        lats.comment = 'For use in GrADS, use the 2D "latitude" for coordinates'
        lats.missing_value = MAPL_UNDEF
        lats.units = 'degrees_north'
        lat1, lat2 = min(self.lat[0],self.lat[-1]), max(self.lat[0], self.lat[-1])
        lats[:] = linspace(lat1,lat2,NS)

        # eta levels
        lev = hf.createVariable('lev','f4',('lev',),fill_value=MAPL_UNDEF)
        lev.standard_name = 'model_layers'
        lev.long_name = 'vertical level'
        lev.missing_value = MAPL_UNDEF
        lev.units = 'layer'
        lev.positive = 'down'
        lev.coordinate = 'eta'
        lev[:] = np.arange(1,NK+1)

        # actual horizontal coords
        longitude = hf.createVariable('longitude','f4',('number_of_scans','ccd_pixels'),fill_value=MAPL_UNDEF)
        longitude.standard_name = 'longitude'
        longitude.long_name = "Longitudes of pixel locations"
        longitude.units = 'degrees_east'
        longitude.missing_value = MAPL_UNDEF
        tmp = np.zeros(self.orgshape)
        tmp[~self.offview] = self.lon
        tmp[ self.offview] = MAPL_UNDEF
        longitude[:,:] = tmp

        latitude = hf.createVariable('latitude','f4',('number_of_scans','ccd_pixels'),fill_value=MAPL_UNDEF)
        latitude.standard_name = 'latitude'
        latitude.long_name = "Latitudes of pixel locations"
        latitude.units = 'degrees_north'
        latitude.missing_value = MAPL_UNDEF
        tmp = np.zeros(self.orgshape)
        tmp[~self.offview] = self.lat
        tmp[ self.offview] = MAPL_UNDEF
        latitude[:,:] = tmp


        # actual time in seconds past hour
        time = hf.createVariable('ev_mid_time','f4',('number_of_scans','ccd_pixels'),fill_value=MAPL_UNDEF)
        time.long_name = "Earth view mid time (seconds of day)"
        time.units = 'seconds'
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
                v = hf.createVariable(var,'f4',('number_of_scans','ccd_pixels',),
                                      fill_value=MAPL_UNDEF,zlib=zlib)
            elif rank == 2:
                v = hf.createVariable(var,'f4',('lev','number_of_scans','ccd_pixels'),
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



if __name__ == "__main__":

    # setup
    # starttime = '2006-03-24T00:50'
    # endtime   = '2006-03-24T00:50'
    dt        = relativedelta(minutes=5)
    # levelB    = '/discover/nobackup/projects/gmao/osse2/pub/c1440_NR/OBS/PACE/LevelB'
    # levelC    = '/discover/nobackup/projects/gmao/osse2/pub/c1440_NR/OBS/PACE/LevelC'
    # const_x = 'const_2d_asm_x'
    levelB = '/nobackup/3/pcastell/PACE/LevelB'
    levelC = '/nobackup/3/pcastell/PACE/LevelC'
    const_x = 'const_2d_asm_x_C'

    # mode of ICA? (from MODES)
    #   mode = 'HOMOCLD-MAXRAN' 
    #   mode = 'HOMOCLD-COSP' 
    mode = 'TOTWPDF-GCOP-SKEWT'


#   Parse command line options
#   --------------------------
    parser = OptionParser(usage="Usage: %prog [OPTIONS] iso_t1 [iso_t2]",
                          version='1.0.0' )

    parser.add_option("--levelB", dest="levelB", default=levelB,
              help="levelB directory (default=%s)"%levelB )

    parser.add_option("--levelC", dest="levelC", default=levelC,
              help='LevelC (output) directory (default=%s)'%levelC)    

    (options, args) = parser.parse_args()
    
    if len(args) == 2:
        iso_t1, iso_t2 = args
        starttime, endtime = (isoparser(iso_t1), isoparser(iso_t2))
    elif len(args) == 1:
        iso_t1 = args[0]
        iso_t2 = None
        starttime    = isoparser(iso_t1)
        endtime      = isoparser(iso_t1)
    else:
        parser.error("must have 1 or 2 arguments: iso_t1 [iso_t2]")


    tyme   = starttime
    while tyme <= endtime:
        outdir  = os.sep.join((options.levelC,'Y'+tyme.strftime('%Y'),'M'+tyme.strftime('%m'),'D'+tyme.strftime('%d')))
        if not os.path.exists(outdir): os.makedirs(outdir)

        print '<> Working on ',tyme

        # instantiate PACE granule
        coll = 'aer_Nv','asm_Nx','met_Nv','chm_Nv'
        granules = granulesLB(options.levelB,tyme,tyme,coll)
        coll = 'cld_nearest'
        granulesN = granulesLBN(options.levelB,tyme,tyme,coll)
        g = PCS(granules,granulesN,const_x)

        # do the simulation and write out results
        ofile = os.sep.join((outdir,'pace-g5nr.'+mode+'.'+tyme.strftime('%Y%m%d_%H%M00')+'.nc4'))
        print '<> creating ',ofile
        g.doICA(ofile,mode,clump=True)

        tyme += dt
