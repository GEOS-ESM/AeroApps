"""

MODIS Cloud Simulator.

"""

from numpy import savez, ones, zeros, array, size, sort, sum, \
                  linspace, any, log, mean, tile, amin, amax, exp
from scipy.stats.mstats import gmean

from   pyobs.mxd03  import MxD03, MxD03Handle, DATE_START
from   pyobs        import NPZ


from   ica          import genICA, getRe, getTau, getQsat
from mcsRegrid_     import getedgecoords, interppressure, \
                           regridtracer, regridmet, edget, edgeq

from MAPL.constants import *
from MAPL           import config

import os

LINEAR_VARS  = dict (
                 asm_Nx = ('PS', 'TS', 'U10M', 'V10M',),
                 asm_Nv = ('PHIS','T', 'U','DELP','O3','QV'),
               )

NEAREST_VARS = dict (
                 lnd_Nx = ('FRSNO',),
                 ocn_Nx = ('FRSEAICE',),
                 asm_Nv = ('T','DELP','QV','QI','QL',),
                 rad_Nv = ('CLOUD',),
                 ice_Nv = ('QILS', 'QIAN',),
#                cld_Nv = ('TAUCLI', 'TAUCLW', 'QILS', 'QIAN',),
               )

DISORT_LEVELS = array([80, 60, 50, 45, 40, 35, 30, 25, 20,
                       18, 16, 15, 14, 13, 12, 11, 10, 9,
                       8, 7, 6, 5, 4, 3, 2, 1, 0])

# Common Variable names for output file; TAU and RE are added later
LINEAR_NAMES = ( 'U10M', 'V10M', 'TS')
NEAREST_NAMES = ( 'FRSNO', 'FRSEAICE' )
ICA_NAMES = ()

# dictionary of possible modes and their description
MODES = \
  {'HOMOCLD-MAXRAN':     'homogeneous clouds, maximum-random overlap',
   'HOMOCLD-COSP':       'homogeneous clouds, COSP overlap',
   'TOTWPDF-GCOP-SKEWT': 'total water PDF, Gaussian copula overlap, Skewed-Triangle PDF'
  }

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

        asm_Nx = rootdir + '/inst3_2d_asm_Nx'
        lnd_Nx = rootdir + '/tavg1_2d_lnd_Nx'
        ocn_Nx = rootdir + '/tavg1_2d_ocn_Nx'

        asm_Nv = rootdir + '/inst3_3d_asm_Nv'
        rad_Nv = rootdir + '/tavg3_3d_rad_Nv'
        cld_Nv = rootdir + '/tavg3_3d_cld_Nv'

        self.getICAindx(asm_Nx)

        # Linear interpolation
        # --------------------
        self.linearSampleFile(asm_Nx,onlyVars=LINEAR_VARS['asm_Nx'])
        self.linearSampleFile(asm_Nv,onlyVars=LINEAR_VARS['asm_Nv'])
        self.sample.ZS = self.sample.PHIS / MAPL_GRAV

        self.linear = self.sample
        self.sample = None

        # Nearest neighbors - 2D
        # ----------------------
        self.nearestSampleFile(lnd_Nx,onlyVars=NEAREST_VARS['lnd_Nx'])
        self.nearestSampleFile(ocn_Nx,onlyVars=NEAREST_VARS['ocn_Nx'])

        # Nearest neighbors - 3D
        # ----------------------
        self.nearestSampleFile(rad_Nv,onlyVars=NEAREST_VARS['rad_Nv'])
        self.nearestSampleFile(asm_Nv,onlyVars=NEAREST_VARS['asm_Nv'])
        self.nearestSampleFile(cld_Nv,onlyVars=NEAREST_VARS['ice_Nv'])

        self.nearest = self.sample
        self.sample = None

        # Hack: to avoid extrapolation issues, pretend the GEOS-5 top is
        #       highter then it actually is
        # --------------------------------------------------------------
        self._fixPTOP()

        # Save NPZ files, if desired (useful for debugging)
        # -------------------------------------------------
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
        Read previously sampled GEOS-5 variables (on MODIS 1km swath).        
        """

        # index pixels to GEOS grid
        self.getICAindx(rootdir + '/inst3_2d_asm_Nx')

        # read npz files from earlier run
        self.linear  = NPZ(npzLinear)
        self.nearest = NPZ(npzNearest)

        # validate and store globals
        if not (self.linear.PTOP == self.nearest.PTOP):
          raise ValueError, 'inconsistent <PTOP> between NPZ files'
        else:
          self.PTOP = self.linear.PTOP

#---
    def genICA(self, mode, clump=True, debug=False, **kwopts):
        """
        Calculate TAU and RE using ICA generated subcolumns.
        """
        n = self.nearest
        l = self.linear
        ica = self.ica

        # Space for optical properties
        # ----------------------------
        shp = n.T.shape
        for v in ( 'QV', 'QL', 'QI', 'RHMCS' ):
            ica.__dict__[v] = MAPL_UNDEF * ones(shp)
 
        if debug:
          for v in ( 'MQV', 'MQL', 'MQI','TAULGM', 'TAUIGM' ):
            ica.__dict__[v] = MAPL_UNDEF * ones(shp)
          for v in ( 'ADECZ', 'ADECZM' ):
            ica.__dict__[v] = MAPL_UNDEF * ones(shp[:-1])

        if self.verb:
            print " <> Performing ICA calculations"

        # Loop over unique gridboxes
        # --------------------------
        m, ng = 0, len(ica.Indices.keys())
        p_ = -1
        for gi,gj in ica.Indices.values():

            I = (ica.iCoord==gi) & (ica.jCoord==gj)
            ii, jj = I.nonzero()        # pixels on this gridbox
            i, j = ii[0], jj[0]         # index of one of the profiles
            ncols = ii.size             # number of columns

            # Select profile
            # --------------
            DELP  = n.DELP[i,j,:]
            T     = n.T[i,j,:]
            QV    = n.QV[i,j,:]
            QL    = n.QL[i,j,:]
            QI    = n.QI[i,j,:]
            CLOUD = n.CLOUD[i,j,:]

            if self.verb:
              m += 1; p = int(100.*m/ng)
              if p > p_:
                print '    - working on (%4d,%4d) with ncols = %4d, %3d%%'%(gi,gj,ncols,p)
              p_ = p

            # Generate condensate subcolumns
            # ------------------------------
            if clump:
              ica.QV[I],ica.QL[I],ica.QI[I] = \
                genICA(ncols,DELP,T,QV,QL,QI,CLOUD,self.PTOP,mode,
                       self.Longitude[I], self.Latitude[I])
            else:
              ica.QV[I],ica.QL[I],ica.QI[I] = \
                genICA(ncols,DELP,T,QV,QL,QI,CLOUD,self.PTOP,mode)

            # generate relative humidities on subcolumns
            # ------------------------------------------
            ica.RHMCS[I] = ica.QV[I] / tile(getQsat(self.PTOP,DELP,T),(ncols,1))

        # Consistency check
        # -----------------
        if any(ica.QV==MAPL_UNDEF) or any(ica.QL==MAPL_UNDEF) or any(ica.QI==MAPL_UNDEF):
            raise MCSError, "Very strange: ICA procedure failed to generate some subcolumns"

        # Effective radius -- use linear interp fields
        # --------------------------------------------
        ica.REL, ica.REI = getRe(l.DELP,l.T,l.U,n.QILS,n.QIAN,self.PTOP)

        # Tau -- use ICA sampled QL/QI with linear interp RE
        # --------------------------------------------------
        ica.TAUL,ica.TAUI = getTau(l.DELP,ica.REL,ica.REI,ica.QL,ica.QI)

        # extra diagnostic variables for debugging
        if debug:

          if self.verb: print " <> extra degugging diagnostics"

          # average the ICA results back to the grid-box
          # --------------------------------------------
          m, p_ = 0, -1
          for gi,gj in ica.Indices.values():

            I = (ica.iCoord==gi) & (ica.jCoord==gj)
            ncols = I.nonzero()[0].size

            if self.verb:
              m += 1; p = int(100.*m/ng)
              if p > p_:
                print '    - working on (%4d,%4d) with ncols = %4d, %3d%%'%(gi,gj,ncols,p)
              p_ = p

            # water content means
            ica.MQV[I] = tile(mean(ica.QV[I],axis=0),(ncols,1))
            ica.MQL[I] = tile(mean(ica.QL[I],axis=0),(ncols,1))
            ica.MQI[I] = tile(mean(ica.QI[I],axis=0),(ncols,1))

            # TAU means
            ica.TAULGM[I] = tile(gmean(ica.TAUL[I],axis=0),(ncols,1))
            ica.TAUIGM[I] = tile(gmean(ica.TAUI[I],axis=0),(ncols,1))

            # conservative delta-Eddington albedo with zenith sun (ADECZ)
            # -----------------------------------------------------------
            mu0 = 1.; g = 0.85; f = g**2
            # [I] collapses horizontal dimensions to 0th-dim
            ica.ADECZ[I] = albedo_dE_cons(
              sum(ica.TAUL[I]+ica.TAUI[I],axis=1), # sum over z
              mu0,g,f)
            ica.ADECZM[I] = mean(ica.ADECZ[I])

          # layer water paths (DTWP, DCWP)
          # ------------------------------
          ica.DCWP  = (          ica.QL  + ica.QI ) * n.DELP / MAPL_GRAV
          ica.DTWP  = (ica.QV +  ica.QL  + ica.QI ) * n.DELP / MAPL_GRAV
          ica.DCWPM = (          ica.MQL + ica.MQI) * n.DELP / MAPL_GRAV
          ica.DTWPM = (ica.MQV + ica.MQL + ica.MQI) * n.DELP / MAPL_GRAV 

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
        if self.verb:
            print " - Getting GEOS-5 edge coordinates"
        g5.PE, g5.ZE = getedgecoords(g5.T,g5.QV,g5.DELP,g5.ZS,self.PTOP)

        g5.logPE = None # compute this on demand later

        # Get DISORT edge pressure
        # ------------------------
        if self.verb:
            print " - Getting DISORT edge pressure"
        di.ZE = 1000 * DISORT_LEVELS  # in meters
        di.PE,di.PS,di.KS,rc = interppressure(di.ZE,self.ZS,g5.PE,g5.ZE)
        if rc:
            raise MCSError, 'Error on return from interpPressure(), rc = <%d>'%rc

        di.DELP = di.PE[:,:,1:]-di.PE[:,:,0:-1]

        di.logPE = None # compute this on demand later

        # Adjust skin temperature for terrain
        # -----------------------------------
        g5.TS = g5.TS * (di.PS/g5.PS)**MAPL_KAPPA

        # Notice that *ks* is the surface level, 1-offset as in Fortran
        # ps is the surface pressure which has been adjusted for terrain height

        self.di = di

    def edgeVars(self):
        """
        Get T, QV and O3 on DISORT edge levels.
        """
        
        if self.verb:
            print " - Getting GEOS-5 edge variables"

        di = self.di
        g5 = self.linear

        # Edge temperature
        # ----------------
        di.TE, rc = edget(di.KS,di.PE,g5.PE,g5.T,g5.TS)
        if rc:
            raise MCSError, 'Error on return from edgeT(), rc = <%d>'%rc

        # Edge specific humidity
        # ----------------------
        di.QE, rc = edgeq(di.KS,di.PE,g5.PE,g5.QV)
        if rc:
            raise MCSError, 'Error on return from edgeQ() for QV, rc = <%d>'%rc
            
        # Convert O3 from mass mixing ratio to volume mixing ratio (ppmv)
        # ---------------------------------------------------------------
        M_air = 28.9        # kg / mol
        M_o3  = 3. * 15.999 # kg / mol
        mmr2vmr = (M_air / M_o3 ) * 1E6 # kg/kg -> ppmv 
        g5.O3 = mmr2vmr * g5.O3

        # Edge ozone
        # ----------
        di.O3E, rc = edgeq(di.KS,di.PE,g5.PE,g5.O3)
        if rc:
            raise MCSError, 'Error on return from edgeQ() for O3, rc = <%d>'%rc
            

    def regridTracer(self,q):
         """
         Regrid a tracer from GEOS-5 to disort vertical coordinates.
         This uses a mass conserving scheme.
         """
         
         di = self.di
         g5 = self.linear
         di_q = regridtracer(di.PE,di.KS,q,g5.PE)

         return di_q

    def regridTau(self,tau):
         """
         Regrid cloud optical depth from GEOS-5 to DISORT vertical coordinates.
         First we use a tracer mass conserving scheme to interpolate
         tau/DELP to the DISORT vertical grid, then multiply the result
         by DELP on the DISORT grid.
         """
         
         di = self.di
         g5 = self.linear
         di_q = regridtracer(di.PE,di.KS,tau/g5.DELP,g5.PE)

         di_tau = di_q * di.DELP
         di_tau[di_q==MAPL_UNDEF] = MAPL_UNDEF

         return di_tau

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

    def doICA(self,filename,mode,clump=True,debug=False):
        """
        Read and write variables doing ICA sampling.
        """

        if self.linear == None:
          raise RuntimeError, 'must sample first using getGEOS()'

        # Calulate TAU and RE
        # -------------------
        self.genICA(mode,clump=clump,debug=debug)
        
        # Prepare for vertical regridding
        # -------------------------------
        self.regridCoords()

        # Edge variables: T, QV and O3
        # ----------------------------
        self.edgeVars()

        # Write the output file
        # ---------------------
        icaNames = ('TAUL','TAUI','REL','REI','RHMCS')
        if debug:
          icaNames += ('DCWP','DTWP','DCWPM','DTWPM', 
            'TAULGM','TAUIGM','ADECZ', 'ADECZM')
        self.writeNC(filename, icaNames=icaNames,
          title = 'GEOS-5 MODIS Cloud Simulator' 
            + ' with ICA sampling (%s)' % MODES[mode])

#---
    def writeNC ( self, filename, format='NETCDF4', zlib=False,
                  linearNames=LINEAR_NAMES, nearestNames=NEAREST_NAMES,icaNames=ICA_NAMES,
                  title='GEOS-5 MODIS Cloud Simulator',
                  rcVars='variables.rc', Verbose=True):
        """
        Write a NetCDF file with simulated GEOS-5 variables on the swath grid.
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
        nc.comment = 'This file contains GEOS-5 cloud related parameters sampled on a MODIS granule swath.'
        nc.contact = 'Arlindo da Silva <arlindo.dasilva@nasa.gov>'
        nc.Conventions = 'CF'
        nc.mod03 = self.path

        # Dimension sizes
        # 
        Ns, Nf = self.lon.shape
        Nk = len(DISORT_LEVELS)-1

        # Create dimensions
        # -----------------
        ns = nc.createDimension('nscans',Ns)
        nf = nc.createDimension('nframes',Nf)
        nk = nc.createDimension('height',Nk)
        ne = nc.createDimension('height_edge',Nk+1)

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

        # Skip standard names for now....

        # Terrain Height
        # --------------
        ZS = nc.createVariable('ZS','f4',('nscans','nframes'))
        ZS.long_name = 'Terrain Height'
        ZS.units = 'm'
        ZS.missing_value = MAPL_UNDEF
        ZS[:,:] = self.ZS[:,:]

        # Terrain Pressure
        # ----------------
        PS = nc.createVariable('PS','f4',('nscans','nframes'))
        PS.long_name = 'Terrain Corrected Surface Pressure'
        PS.units = 'Pa'
        PS.missing_value = MAPL_UNDEF
        PS[:,:] = self.di.PS[:,:]

        # Surface Index
        # -------------
        KS = nc.createVariable('K_SURFACE','I4',('nscans','nframes'))
        KS.long_name = 'DISORT Surface Index (Fortran 1-offset)'
        KS.comments = 'Valid levels for each pixel are from 1 thru K_SURFACE'
        KS.units = '1'
        KS.missing_value = int(MAPL_UNDEF)
        KS[:,:] = self.di.KS[:,:]

        # Edge Pressure
        # -------------
        PE = nc.createVariable('PE','f4',('nscans','nframes','height_edge'))
        PE.long_name = 'DISORT Edge Pressures'
        PE.comments = 'Valid levels for each pixel are from 1 thru K_SURFACE'
        PE.units = 'Pa'
        PE.missing_value = MAPL_UNDEF
        PE[:,:,:] = self.di.PE[:,:,:]

        # Edge Temperature
        # ----------------
        TE = nc.createVariable('T','f4',('nscans','nframes','height_edge'))
        TE.long_name = 'Temperature at DISORT Edge Levels'
        TE.comments = 'Valid levels for each pixel are from 1 thru K_SURFACE'
        TE.units = 'K'
        TE.missing_value = MAPL_UNDEF
        TE[:,:,:] = self.di.TE[:,:,:]

        # Edge Specific Humidty
        # ----------------------
        QE = nc.createVariable('QV','f4',('nscans','nframes','height_edge'))
        QE.long_name = 'Specific Humidity at DISORT Edge Levels'
        QE.comments = 'Valid levels for each pixel are from 1 thru K_SURFACE'
        QE.units = 'Kg Kg-1'
        QE.missing_value = MAPL_UNDEF
        QE[:,:,:] = self.di.QE[:,:,:]

        # Edge Ozone
        # ----------
        O3E = nc.createVariable('O3','f4',('nscans','nframes','height_edge'))
        O3E.long_name = 'Ozone Mass Mixing Ratio at DISORT Edge Levels'
        O3E.comments = 'Valid levels for each pixel are from 1 thru K_SURFACE'
        O3E.units = 'Kg Kg-1'
        O3E.missing_value = MAPL_UNDEF
        O3E[:,:,:] = self.di.O3E[:,:,:]

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

        # ICA variables
        # -------------
        self._writeVars(nc,cf,icaNames,self.ica,
                        comment='ICA sampled variable',
                        Verbose=Verbose,zlib=zlib)

        # Close the file
        # --------------
        nc.close()

    def _writeVars(self,nc,cf,Names,data,comment=None,Verbose=True,zlib=True):
        """
        Write variables for a given group (linear, nearest or ica).
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
            else:
                raise MCSError, 'Invalid rank for variable <%s>'%var
            v.long_name = long_name
            v.units = units.strip()
            v.missing_value = MAPL_UNDEF
            if comment != None:
                v.comment = comment
            if rank == 2:
                v[:,:] = data.__dict__[var][:,:]
            else:
                if var in ('TAUL','TAUI',                       
                           'DCWP','DTWP','DCWPM','DTWPM','TAULGM','TAUIGM'):
                    v[:,:,:] = self.regridTau(data.__dict__[var])[:,:,:] 
                elif var in ('REL','REI'):
                    # convert to micron here
                    v[:,:,:] = self.regridLogP(1e6*data.__dict__[var])[:,:,:] 
                elif var in ('RHMCS'):
                    # for now, use logP interpolation for RH,
                    # typically a more smooth variable than qv.
                    v[:,:,:] = self.regridLogP(data.__dict__[var])[:,:,:] 
                else:
                    v[:,:,:] = self.regridTracer(data.__dict__[var])[:,:,:] 

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
def calculon():
    rootdir = '/nobackup/fp/opendap'
    mxd03   = '/nobackup/MODIS/Level1/MYD03/sample/MYD03.A2012228.1200.005.2012229145830.hdf'    
    outdir  = os.sep.join((os.environ['NOBACKUP'],'4',os.environ['USER'],'MCS','output'))
    return (rootdir,mxd03,outdir)

#def nccs():
#   rootdir = '/home/adasilva/iesa/aerosol/data/MCS/opendap'
#   mxd03   = '/home/adasilva/iesa/aerosol/data/MODIS/Level1/051/MOD03/2012/228/MYD03.A2012228.1200.005.2012229145830.hdf'    
#   outdir  = os.sep.join((os.environ['NOBACKUP'],'MCS','output'))

def nccs():
    rootdir = '/home/pnorris/iesa/aerosol/data/MCS/opendap/e5110_ugne_case2ia/20130112'
    mxd03   = '/discover/nobackup/gwind/WMO_cases/MYD03.A2013012.0450.006.2013012170123.hdf'
    tag =  os.path.splitext(os.path.basename(mxd03))[0]
    outdir  = os.sep.join((os.environ['NOBACKUP'],'MCS','output_'+tag))
    return (rootdir,mxd03,outdir)

if __name__ == "__main__":

    # setup
    rootdir, mxd03, outdir = nccs()
    os.mkdir(outdir)

    # instantiate MxD03 granule
    m = MCS(mxd03,Verb=1)

    # associate GEOS fields with granule
    # optionally write npz files for future bypass of this step
    m.getGEOS(rootdir, 
        npzLinear =os.sep.join((outdir,'geos-linear.npz')),
        npzNearest=os.sep.join((outdir,'geos-nearest.npz')))

#   # OR ... read sampled GEOS fields from old npz files
#   m.getGEOS_readFromNPZ(rootdir,
#       os.sep.join((outdir,'geos-linear.npz')),
#       os.sep.join((outdir,'geos-nearest.npz')))

    # mode of ICA? (from MODES)
#   mode = 'HOMOCLD-MAXRAN' 
#   mode = 'HOMOCLD-COSP' 
    mode = 'TOTWPDF-GCOP-SKEWT'

    # do the simulation and write out results
    ofile = os.sep.join((outdir,'geos-icacl-'+mode+'.nc4'))
    m.doICA(ofile,mode,clump=True)
