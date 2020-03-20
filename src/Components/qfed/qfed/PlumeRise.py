"""
A Python interface to CPTEC's Plume Rise Model.

TO DO: Create a base class name PLUMES and share code with the LockPlume.

"""

from numpy      import zeros, ones, meshgrid, linspace, any, \
                       pi, sin, cos, arccos, arange, array, \
                       savez, NaN, isnan

from scipy import optimize as opt

from datetime   import datetime, date, timedelta
from glob       import glob

from dozier     import DOZIER, granules
from PlumeRise_ import *       # f2py extension
from gfio       import GFIO

from pyobs           import NPZ, kde
from pyobs.binObs_   import binareas
from pyobs.minx      import MINXs
from MAPL.constants  import *

import eta

__VERSION__ = 2.1
__CVSTAG__  = '@CVSTAG'
__AMISS__   = 1.E+20

Bioma = [ 'Tropical Forest',
          'Extra-Tropical Forest',
          'Savanna',
          'Grassland' ]

DAY = timedelta(seconds=60*60*24)

#----------------------------------------------------------------------------------------------

#----
class MINXs_PR(MINXs):

    """
    Extension of the MINXs class adding the Freitas Plume Rise
    functionality. This class handles non-gridded, observation
    location fires.

    Parabolic Vertical Mass Distribution (VMD)
    ------------------------------------------
    If using the parabolic VMD as in getVMD() below, the (z_c,delta) parameters can
    be computed from the (z_i,z_f,z_d) parameters returned by plumevmd().

     a) Like Saulo:
                     z_c   = (z_f+z_i)/2
                     delta = (z_f-z_i)/2

     b) Preserve bottom half:
                     z_c   = z_d
                     delta = z_d - z_i

     c) Preserve upper half:
                     z_c   = z_d
                     delta = z_f - z_d

                       ---

    """
    
    def getPlume1(self,i,
                  hflux_kW=None,
                  frp_MW=None,
                  area_m2=None,afac=None,
                  Area=None,
                  rad2conv=5.,
                  Nominal=False,
                  which='z_a',
                  Verbose=False):
 
        """
        Runs the Plume Rise extension to compute the extent of the plume.
        On input,
        
           i         ---  index to operate on
           frp_MW    ---  fire radiative power in MW
           area_m2   ---  firea are im m2
           afac      ---  scale area by this amount
           Area      ---  if specified, PR model will see this value, but heat flux
                          will be computed the regular way, based on *area*
           rad2conv  ---  factor for converting radiative to convective heat fluxes.
           which     ---  Return one of z_i, z_f, z_d, z_a or all (Tupple), where
                          z_i -- height of maximum W (bottom of plume)
                          z_d -- height of maximum detrainment
                          z_a -- average height in (z_i,z_f), weighted by -dw/dz
                          z_f -- height where w<1 (top of plume)
           Nominal   ---  by default, areas and heat fluxes are from a
                          Dozier type algorithm, If Nominal=True,
                             area = 20e4 (20 ha)
                             hflux = rad2conv * FRP / area
           
	"""

        area = area_m2
        ptop = 1. # Pascal

        # Default areas, heat flux
        # ------------------------
        if Nominal:
            d_area = 20e4  #  typical fire size in m2 (20 ha)
        else: 
            d_area = 1e6 * self.mod14.farea[i]      # km2 --> m2

        # Area
        # ----
        if area == None:
            area = d_area # m2

        # FRP
        # ---
        if frp_MW == None:
            frp_MW = self.mod14.frp[i] # MW

        # Convective Heat Flux
        # --------------------
        if hflux_kW == None:
            hflux_kW = 1e3 * rad2conv * frp_MW / area # kW/m2

        # Scaled Area
        # -----------
        if afac!=None:
            if afac==0:
                afac = 1e-8
            area = afac * area
            hflux_kW = hflux_kW / afac
        
        u = self.sample.u[i]
        v = self.sample.v[i]
        T = self.sample.t[i]
        q = self.sample.qv[i]
        delp = self.sample.delp[i]

        if delp.min()<=0 or T.min()<=0:
            print "out of range: ", delp.min(), T.min() 
            return NaN

        # Ensure arrays are top-down as in GEOS-5
        # ---------------------------------------
        if delp[0] > delp[-1]:
            u = u[::-1]
            v = v[::-1]
            T = T[::-1]
            q = q[::-1]
            delp = delp[::-1]

        # Override firea area
        # -------------------
        if Area is not None:
            area = Area
                
        # Run plume rise model
        # --------------------
        # p,z,k,rc = plume(u,v,T,q,delp,ptop,hflux_kW,area)
        z_i,z_d,z_a,z_f,z,w,rc = plumevmd(u,v,T,q,delp,ptop,hflux_kW,area)
        if rc:
            raise ValueError, "error on return from <plume>, rc = %d"%rc

        if z_i==-1.: 
             z_i,z_d,z_a,z_f = (NaN,self.sample.pblh[i],NaN,NaN)
             
        if   which=='z_i': z = z_i
        elif which=='z_f': z = z_f
        elif which=='z_d': z = z_d
        elif which=='z_a': z = z_a
        else:              z = (z_i,z_d,z_a,z_f,z,w)
        
        if Verbose:
            print " ", self.tyme[i], "| %8.2f | %8.2f %8.2f | %8.2f %8.2f %8.2f %8.2f | %03d"%\
                (self.mod14.fdist[i], area/1e4, hflux_kW, z_i,z_d,z_a,z_f,i)

        return z

#---
    def getPlume(self,I=None,Verbose=True,**kwopts):
 
        """
        Runs the Plume Rise extension to compute the extent of each plume.
	    """

        z_plume = __AMISS__ * ones(self.N)

        R = arange(self.N)
        if I is not None:
            R = R[I]

        if Verbose:
            print ""
            print "                    Plume Height Estimates"
            print ""
            print "  --------------------|----------|-------------------|-------------------------------------------"
            print "                      | Distance |   Fire Properties |             Plume Height"
            print "    MINX Date/Time    |  to Fire |   Area   Heat Flx |    z_i      z_d      z_a      z_f"
            print "                      |    km    |    ha      kW     |    km       km       km       km"
            print "  --------------------|----------|-------------------|-------------------------------------------"
            
#       Loop over time
#       --------------
        for i in R:
            z_plume[i] = self.getPlume1(i,Verbose=Verbose,**kwopts)
        if Verbose:
            print "  --------------------|----------|-------------------|-------------------------------------------"

        return array(z_plume)

        if Verbose:
            print "  --------------------|----------|-------------------|----------"

#---
    def getFires(self,mod14_path='/home/adasilva/iesa/aerosol/data/MODIS/Level2/MOD14',
                      method='classic',
                      npzFile=None,Verbose=True):
        """
        Retrieves Level2 MOD14 fire data for each MINX fire.
        """
        from dozier import DOZIER

        self.mod14 = []
        d2r = pi / 180.
        a = MAPL_RADIUS/1000. # Earth radius in km
        dt = timedelta(seconds=5*60)

        self.mod14 = MOD14(self.N)
        
        if Verbose:
            print ""
            print "                   Fire Heat Flux Estimates"
            print ""
            print "  --------------------|----------|-------------------|----------"
            print "                      | Distance |    FRP Estimates  |   Fire"
            print "    MINX Date/Time    |  to Fire |   MINX      MODIS | Heat Flux"
            print "                      |    km    |    MW         MW  |  kW/m2"
            print "  --------------------|----------|-------------------|----------"

        for i in range(self.N):

            # Get fire granule for this particular day
            # ----------------------------------------
            t = self.tyme[i]
            p = mod14_path
            m = DOZIER(_getGran(t-dt,p) + _getGran(t,p) + _getGran(t+dt,p))

            # Select closest fires
            # --------------------
            x0 = cos(d2r*self.lat_f[i]) * cos(d2r*self.lon_f[i])
            y0 = cos(d2r*self.lat_f[i]) * sin(d2r*self.lon_f[i])
            z0 = sin(d2r*self.lat_f[i])
            dx = x0*cos(d2r*m.lat) * cos(d2r*m.lon)
            dy = y0*cos(d2r*m.lat) * sin(d2r*m.lon)
            dz = z0*sin(d2r*m.lat) 
            s  = a * arccos(dx+dy+dz) # great circle distance
            j = s.argmin()

            # Estimate fire heat flux
            # -----------------------
            if method=='bimodal':
                m.bimodal_u()   # bimodal Dozier, need to provide trasnmittance
                self.mod14.hflux[i] = m.h_F[j] # flaming Radiative Heat Flux
                self.mod14.fdist[i] = s[j]
                self.mod14.frp[i]  = m.pow_F[j]
                self.mod14.pixar[i] = m.pixar[j]
                self.mod14.farea[i] = m.r_F[j]
                self.mod14.qa[i] = m.m[j] # bolean
            else:
                m.classic_var() # classic Dozier, need to provide trasnmittance
                self.mod14.hflux[i] = m.hflux[j] # Radiative Heat Flux
                self.mod14.fdist[i] = s[j]
                self.mod14.frp[i]  = m.pow[j]
                self.mod14.pixar[i] = m.pixar[j]
                self.mod14.farea[i] = m.farea[j]
                self.mod14.qa[i] = m.m[j] # bolean
            if Verbose:
                print " ", t, "| %8.2f | %8.2f %8.2f | %8.2f"%\
                    (self.mod14.fdist[i], self.mod14.frp[i], self.mod14.frp[i], \
                     self.mod14.hflux[i] )

        # Save fire properties in NPZ file for later
        # ------------------------------------------
        if npzFile!=None:
            savez(npzFile,**self.mod14.__dict__)
                  
        if Verbose:
            print "  --------------------|----------|-------------------|----------"

#---
    def getOptBrute(self,I=None,Verbose=True,**kwopts):
 
        """
        Runs the Plume Rise extension to compute brute force optimal value of (hflux,area)
        to match MISR plume height.
	    """

        z_opt = __AMISS__ * ones(self.N)
        h_opt = __AMISS__ * ones(self.N)
        a_opt = __AMISS__ * ones(self.N)
        
        R = arange(self.N)
        if I is not None:
            R = R[I]

        Hflux_kW = linspace(1,100,10)
        Area = linspace(0.1e4,2e4,10)

        print Hflux_kW
        print Area

        if Verbose:
            print ""
            print "                    Plume Height Estimates"
            print ""
            print "  --------------------|----------|-------------------|----------"
            print "                      | Observed |    Opt Properties |  Optimal "
            print "    MINX Date/Time    |  Height  |   Area   Heat Flx |  Height"
            print "                      |    km    |    ha      kW     |    km"
            print "  --------------------|----------|-------------------|----------"
            

#       Loop over time
#       --------------
        for i in R:
            z = m.z[i]
            e = 1e20 # overestimate
            for hflux_kW in Hflux_kW:
                for area in Area:
                    z_ = self.getPlume1(i,Verbose=False,hflux_kW=hflux_kW,area=area)
                    if isnan(z)==False:
                        e_ = (z-z_)**2
                        if e_<e:
                            e = e_
                            z_opt[i] = z_
                            h_opt[i] = hflux_kW
                            a_opt[i] = area
            if e<1e20:
                if Verbose:
                    print " ", self.tyme[i], "| %8.2f | %8.2f %8.2f | %8.2f"%\
                          (z, a_opt[i]/1e4, h_opt[i], z_opt[i])
                                        
        if Verbose:
            print "  --------------------|----------|-------------------|----------"

        return (z_opt,h_opt,a_opt)

#---
    def getOptF(self,Verbose=True,area_m2=None):
        """
        Find optimal fire modified bstar to match observed plume height.
        """

        if Verbose:
            print ""
            print "    Fire Modified bstar Optimization"
            print ""
            print "Plume | f_opt  |    J     |  Ni | Nf"
            print "------|--------|----------|-----|-----"

        self.f_opt = ones(self.N)
        self.z_opt = ones(self.N)
        for i in range(self.N):
#            xmin, fval, iter, fcalls = opt.brent(CostFuncF,args=(self,i),brack=(1.,4),full_output=True)
            xmin, fval, iter, fcalls = opt.brent(CostFuncF,args=(self,i,area_m2),full_output=True)
            if isnan(fval):
                self.f_opt[i] = NaN
                self.z_opt[i] = NaN
            else:
                rad2conv = xmin**2
                self.f_opt[i] = rad2conv
                self.z_opt[i] = self.getPlume1(i,rad2conv=self.f_opt[i],area_m2=area_m2)

            if Verbose:
                print "%5d | %6.2f | %8.2f | %3d | %3d "%(i, rad2conv, fval, iter, fcalls)

        if Verbose:
            print "------|--------|----------|-----|-----"

#---
    def getOptA(self,Verbose=True,rad2conv=5):
        """
        Find optimal fire modified afac to match observed plume height.
        """

        if Verbose:
            print ""
            print "    Fire Modified afac Optimization"
            print ""
            print "Plume | f_opt  |    J     |  Ni | Nf"
            print "------|--------|----------|-----|-----"

        self.a_opt = ones(self.N)
        self.z_opt = ones(self.N)
        for i in range(self.N):
            xmin, fval, iter, fcalls = opt.brent(CostFuncA,args=(self,i,rad2conv),full_output=True) 
            afac = xmin**2
            self.a_opt[i] = afac
            self.z_opt[i] = self.getPlume1(i,afac=self.a_opt[i],rad2conv=rad2conv)
            if Verbose:
                print "%5d | %6.2f | %8.2f | %3d | %3d "%(i, afac, fval, iter, fcalls)

        if Verbose:
            print "------|--------|----------|-----|-----"

#---
    def getOptFanneal(self,Verbose=True):
        """
        Find optimal fire hflux scaling to match observed plume height.
        """

        if Verbose:
            print ""
            print "    Fire Modified bstar Optimization"
            print ""
            print "Plume | f_opt  |    J     |  Ni | Nf"
            print "------|--------|----------|-----|-----"

        self.f_opt = ones(self.N)
        self.z_opt = ones(self.N)
        for i in range(self.N):
            xmin, fval, T, fcalls, iter, accept, retval = opt.anneal(CostFuncF,1.0,args=(self,i),full_output=True)
            ffac = xmin**2
            self.f_opt[i] = ffac
            self.z_opt[i] = self.getPlume1(i,ffac=self.f_opt[i])
            if Verbose:
                print "%5d | %6.2f | %8.2f | %3d | %3d "%(i, ffac, fval, iter, fcalls)

        if Verbose:
            print "------|--------|----------|-----|-----"

#---
    def getOptFbnd(self,Verbose=True):
        """
        Find optimal fire hflux scaling to match observed plume height.
        """

        if Verbose:
            print ""
            print "    Fire Modified bstar Optimization"
            print ""
            print "Plume | f_opt  |    J     |  Ni | Nf"
            print "------|--------|----------|-----|-----"

        self.f_opt = ones(self.N)
        self.z_opt = ones(self.N)
        for i in range(self.N):
            xmin,fval,ier,fcalls  = opt.fminbound(CostFuncF,0.,2.,args=(self,i),full_output=True)
            ffac = xmin**2
            self.f_opt[i] = ffac
            self.z_opt[i] = self.getPlume1(i,rad2conv=self.f_opt[i])
            if Verbose:
                print "%5d | %6.2f | %8.2f | %3d | %3d "%(i, ffac, fval, fcalls, fcalls)

        if Verbose:
            print "------|--------|----------|-----|-----"

#----
def CostFuncF(f,m,i,area_m2):
    z = m.getPlume1(i,rad2conv=f**2,area_m2=area_m2)
    return 1e-6 * (m.z[i]-z)**2

def CostFuncA(f,m,i,rad2conv):
    z = m.getPlume1(i,afac=f**2,rad2conv=rad2conv)
    return 1e-6 * (m.z[i]-z)**2

#----------------------------------------------------------------------------------------------

class MOD14(object):
    def __init__(self,N):
        """
        Simple container class for fire properties.
        """
        self.hflux = __AMISS__ * ones(N)  # fire heat flux estimate
        self.fdist = __AMISS__ * ones(N)  # distance to detected fire
        self.frp   = __AMISS__ * ones(N)  # nearest MODIS FRP
        self.pixar = __AMISS__ * ones(N)  # pixel area
        self.farea = __AMISS__ * ones(N)  # fire area
        self.qa    = __AMISS__ * ones(N)  # quality flag
        return
    
#----------------------------------------------------------------------------------------------
class PLUME_L2(DOZIER):

    """
    Extension of the MxD14,IGBP and DOZIER classes, adding the
    Plume Rise functionality. This class handles non-gridded,
    observation location fires.
    """

    def getPlume1(self,i,Verbose=False,rad2conv=5,area_m2=None):
        """
        Compute plume height for the ith fire.
        """

        # Fire properties
        # ---------------
        if area_m2 is None:
            area = 1e6 * self.sample.farea[i] # km2 --> m2
        else:
            area = area_m2
        hflux_kW = 1e3 * rad2conv * self.sample.pow[i] / area
            
        # Meteorology
        # -----------        
        ptop = 1. # Pascal
        u = self.sample.u[i]
        v = self.sample.v[i]
        T = self.sample.t[i]
        q = self.sample.qv[i]
        delp = self.sample.delp[i]

        if delp.min()<=0 or T.min()<=0:
            print "out of range: ", delp.min(), T.min() 
            return None

        # Ensure arrays are top-down as in GEOS-5
        # ---------------------------------------
        if delp[0] > delp[-1]:
            u = u[::-1]
            v = v[::-1]
            T = T[::-1]
            q = q[::-1]
            delp = delp[::-1]

        # Run plume rise model
        # --------------------        
        z_i,z_d,z_a,z_f,z,w,rc = plumevmd(u,v,T,q,delp,ptop,hflux_kW,area)
        if rc:
            raise ValueError, "error on return from <plume>, rc = %d"%rc

        if Verbose:
            print "| %8.2f %8.2f %s | %8.2f %8.2f | %8.2f %8.2f %8.2f %8.2f | %03d"%\
                (self.lon[i],self.lat[i], str(self.tyme[i]), area/1e4, hflux_kW, z_i,z_d,z_a,z_f,i)
        
        return (z_i,z_d,z_a,z_f,z,w)

#---
    def getPlume(self,algo='dozier',Verbose=False,**kwargs):
 
        """
        Runs the Plume Rise extension to compute the extent of the plume for each fire.
        It is assumed that the necessary met fields have already been loaded in attribute
        *sample* by method sampleFP or alternative.

        """

        I = self.sample.I             # spatial subsetting
        if len(I) != self.lon.size:
                raise ValueError, 'sampling data appear inconsistent'
            
        N = len(self.lon)     # all obs
        n = len(self.lon[I])  # reduced set in case of regional subsetting
        self.sample.z_i = __AMISS__ * ones(n,dtype='float32')
        self.sample.z_d = __AMISS__ * ones(n,dtype='float32')
        self.sample.z_a = __AMISS__ * ones(n,dtype='float32')
        self.sample.z_f = __AMISS__ * ones(n,dtype='float32')

          
        # Gather fire properties
        # ----------------------
        self.sample.pow   = self.pow[I]
        self.sample.m     = self.m[I]
            
        # Assume that fire area has alread been estimated by classic_var()
        # ----------------------------------------------------------------
        if algo is not None:
            self.sample.farea = self.farea[I]
            if self.algo != algo:
                raise ValueError, 'only dozier algorithm supportted'

        if Verbose:
            print ""
            print "                         Plume Height Estimates"
            print ""
            print "  --------------------|-----------------------|-------------------------------------------"
            print "                      |    Fire Properties    |             Plume Height"
            print "    Lon  Lat  Time    |   Area     Conv Pwr   |    z_i      z_d      z_a      z_f"
            print "                      |    ha         kW      |    km       km       km       km"
            print "  --------------------|-----------------------|-------------------------------------------"
            

#       Loop over time
#       --------------
        s = self.sample # shorthand
        for i in range(n):

#           Interpolate met fields to fire location
#           ---------------------------------------
            if not self.sample.m[i]: continue   # skip over bad points

#           Compute plume rise for this fire
#           --------------------------------
            s.z_i[i],s.z_d[i],s.z_a[i],s.z_f[i],z,q = self.getPlume1(i,Verbose=Verbose,**kwargs)

#---
    def getPlumeDeprecated(self,met):
 
        """
        Runs the Plume Rise extension to compute the extent of the plume.
        On input,

        met --- MET object for computing Met fields

        Notice that the Dozier must have been run and produced the farea
        attribute. In addition, the veg attribute with the biome type must
        have been defined as well.
	"""

        if self.algo != 'dozier':
            raise ValueError, 'only Dozier algorithm supported for now'

        if self.veg is None:
            raise ValueError, 'veg attribute with biome type has not been defined'


#       Initialize Plume Rise ranges
#       ----------------------------
        ntd = met.ntd
        N = self.lon.size
        self.p_plume = zeros((ntd,N,2))
        self.k_plume = zeros((ntd,N,2))
        self.z_plume = zeros((ntd,N,2))
        yyyy = self.yyyy[N/2]
        jjj = self.jjj[N/2]

#       Loop over time
#       --------------
        for t in range(1,ntd+1):
    
#           Interpolate met fields to fire location
#           ---------------------------------------
            met.interp(yyyy,jjj,t,lon=self.lon,lat=self.lat)

#           Compute plume rise for this time
#           --------------------------------
            farea = self.r_F * self.farea # reduce area by flaming fraction
            p, z, k = getPlume(farea,self.veg,met,t,ntd,self.verb)  

#           Save in the approapriate containers
#           -----------------------------------
            self.p_plume[t-1,:,:] = p[:,:]
            self.z_plume[t-1,:,:] = z[:,:]
            self.k_plume[t-1,:,:] = k[:,:]

#...............................................................................

class PLUME_L3(object):

    """
    Extension of the MxD14,IGBP and DOZIER classes, adding the
    Plume Rise functionality. This class handles non-gridded,
    observation location fires.
    """

    def __init__(self,plume,refine=4,res=None):
        """
        Create a gridded Plume Rise object from a Level 2
        PLUME_L2 object *plume*. The grid resolution is
        specified by

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

        After initialization only the FRP weighted average area is
        set, for each gridbox/biome, along with the coordinates of the
        global grid.
        
        """

        N = plume.lon.size
        self.verb = plume.verb

#       Output grid resolution
#       ----------------------
        if res is not None:
            if res=='a': refine = 1 
            if res=='b': refine = 2
            if res=='c': refine = 4
            if res=='d': refine = 8
            if res=='e': refine = 16

#       Lat lon grid
#       ------------
        dx = 5. / refine
        dy = 4. / refine
        im = int(360. / dx)
        jm = int(180. / dy + 1)
        self.im = im
        self.jm = jm
        self.glon = linspace(-180.,180.,im,endpoint=False)
        self.glat = linspace(-90.,90.,jm)
        Lat, Lon  = meshgrid(self.glat,self.glon)  # shape should be (im,jm)

        self.yyyy = plume.yyyy[N/2]
        self.jjj  = plume.jjj[N/2]
        self.date = date((int(self.yyyy),1,1)) + (int(self.jjj) - 1)*DAY
        self.col  = plume.col
        
#       Supperobed fire attributes for each biome
#       These will have 1D arrays for each biome, using
#       a standard sparse matrix storage
#       -----------------------------------------------
        self.bioma = [1,2,3,4]
        NONE       = [None,None,None,None] 
        self.idx   = NONE[:]  # non-zero indices for each biome
        self.area  = NONE[:]  # non-zero areas   for 
        self.r_F   = NONE[:]  # corresponding flaming fraction
        self.lon   = NONE[:]  # corresponding lon
        self.lat   = NONE[:]  # corresponding lat

#       Plume extent to be filled later
#       -------------------------------
        self.p_plume = NONE[:]
        self.k_plume = NONE[:]
        self.z_plume = NONE[:]

#       Grid box average of fire flaming area, weighted by FRP
#       -----------------------------------------------------
        for b in self.bioma:

#           Compute average area in gridbox for this biome
#           ----------------------------------------------
            Frac = zeros((im,jm))
            Area = zeros((im,jm))
            FRP  = zeros((im,jm))
            i = (plume.veg==b)
            if any(i):
                blon = plume.lon[i]
                blat = plume.lat[i]
                bfrac = plume.r_F[i] * plume.pow[i] 
                barea = plume.r_F[i] * plume.farea[i] * plume.pow[i] # notice flaming fraction 
                bfrp  = plume.pow[i]
                Area +=  binareas(blon,blat,barea,im,jm) # to be normalized
                Frac +=  binareas(blon,blat,bfrac,im,jm) # to be normalized
                FRP  +=  binareas(blon,blat,bfrp, im,jm)
                I = (FRP>0.0)
                if any(I):
                    Area[I] = Area[I] / FRP[I]
                    Frac[I] = Frac[I] / FRP[I]

#           Use sparse matrix storage scheme
#           --------------------------------
            I = Area.nonzero()
            if any(I):
                self.area[b-1] = Area[I]
                self.r_F[b-1]  = Frac[I] # average flaming/total energy fraction
                self.lon[b-1]  = Lon[I]
                self.lat[b-1]  = Lat[I]
                self.idx[b-1]  = I       # save indices for going to global grid

    def getPlume(self):
 
        """
        Runs the Plume Rise extension to compute the extent of the plume.
        On input,

        met --- MET object for computing Met Fields at obs location.

	"""

        ntd = met.ntd # number of time steps per day
        self.ntd = ntd
        
#       Loop over bioma
#       ---------------
        for i in range(len(self.bioma)):

#           No data for this bioma, nothing to do
#           -------------------------------------
            if self.idx[i] is None:
                if self.verb>0:
                    print "[x] no data for %s"%Bioma[i] 
                continue

            lon = self.lon[i]
            lat = self.lat[i]
            area = self.area[i]
            N = lon.size
            veg = self.bioma[i] * ones(N)

            p_plume = zeros((ntd,2,N))
            z_plume = zeros((ntd,2,N))
            k_plume = zeros((ntd,2,N))
                            
            if self.verb>0:
                print "[ ] got %d burning gridboxes in %s"%(N,Bioma[i]) 

#           Loop over time 
#           --------------
            for t in range(1,ntd+1):

#               Interpolate met fields to fire locations
#               ----------------------------------------
                met.interp(self.yyyy,self.jjj,t,lon=lon,lat=lat)

#               Compute plume rise for this time, biome
#               ---------------------------------------
                p, z, k = getPlume(area,veg,met,t,ntd,self.verb)  

#               Save in the approapriate containers
#               -----------------------------------
                p_plume[t-1,:,:] = p.T[:,:]
                z_plume[t-1,:,:] = z.T[:,:]
                k_plume[t-1,:,:] = k.T[:,:]

#           Plume extent for this biome (sparse storage)
#           --------------------------------------------
            self.p_plume[i] = p_plume
            self.z_plume[i] = z_plume
            self.k_plume[i] = k_plume

#---
    def write(self,filename=None,dir='.',expid='qfed2',tag=None):
       """
       Writes gridded Area and FRP to file.
       """

       vtitle = {}
       vtitle['fa'] = 'Flaming Area'
       vtitle['ff'] = 'Fraction of Flaming Energy'
       vtitle['p2'] = 'Plume Bottom Pressure'
       vtitle['p1'] = 'Plume Top Pressure' 
       vtitle['z2'] = 'Plume Bottom Height' 
       vtitle['z1'] = 'Plume Top Height' 
       vtitle['k2'] = 'Plume Bottom Vertical Index' 
       vtitle['k1'] = 'Plume Top Vertical Index' 
                 
       vunits = {}
       vunits['fa'] = 'km2'
       vunits['ff'] = '1'
       vunits['p2'] = 'Pa'
       vunits['p1'] = 'Pa'
       vunits['z2'] = 'meter'
       vunits['z1'] = 'meter'
       vunits['k2'] = '1'
       vunits['k1'] = '1'
                 
       btitle = {}
       btitle['tf'] = 'Tropical Forest' 
       btitle['xf'] = 'Extra-Tropical Forest'
       btitle['sv'] = 'Savanna'
       btitle['gl'] = 'Grassland'

#      Create master variable list
#      ---------------------------
       Vname  = []
       Vtitle = []
       Vunits = []
       for v in vtitle.keys():
           vt = vtitle[v]
           vu = vunits[v]
           for b in btitle.keys(): 
               bt = btitle[b]
               Vname.append(v+'_'+b)
               Vtitle.append(vt+' ('+bt+')')
               Vunits.append(v)

#      Global metadata
#      ---------------
       title = 'QFED Level3c v%3.1f (%s) Gridded Plume Rise Estimates' % (__VERSION__, _getTagName(tag))
       source = 'NASA/GSFC/GMAO GEOS-5 Aerosol Group'
       contact = 'arlindo.dasilva@nasa.gov'

#      Time/date handling
#      ------------------
       if self.date is None:
           print "[x] did not find matching files, skipped writing an output file"
           return

       if 24%self.ntd != 0:
           raise ValueError,"invalid number of times per day (%d),"%self.ntd\
                 +"it must be a divisor of 24."
       else:
           dT = 240000/self.ntd # timestep in hhmmss format
           NHMS = range(0,240000,dT)

       nymd = 10000*self.date.year + 100*self.date.month + self.date.day
       nhms = NHMS[0]
       col = self.col

#      Create output file name
#      -----------------------
       if filename is None:
           filename = '%s/%s.plumerise.%s.%d.nc'%(dir,expid,col,nymd)
       self.filename = filename
       f = GFIO()
       f.create(filename,Vname, nymd, nhms,
                lon=self.glon, lat=self.glat,
                vtitle=Vtitle, vunits=Vunits,
                timinc=dT, amiss=__AMISS__,
                title=title, source=source, contact=contact)

#      Write out Plume Rise variables
#      ------------------------------
       d = (self.im,self.jm)
       for t in range(self.ntd):
           nhms = NHMS[t]
           b = 0
           for bn in btitle.keys():
               I = self.idx[b]
               f_area  = self.area[b]
               f_frac  = self.r_F[b]
               p_plume = self.p_plume[b]
               z_plume = self.z_plume[b]
               k_plume = self.k_plume[b]
               _writeOne(f,'fa_'+bn,nymd,nhms,I,f_area, t,0,d)
               _writeOne(f,'ff_'+bn,nymd,nhms,I,f_frac, t,0,d)
               _writeOne(f,'p1_'+bn,nymd,nhms,I,p_plume,t,0,d)
               _writeOne(f,'p1_'+bn,nymd,nhms,I,p_plume,t,0,d)
               _writeOne(f,'p2_'+bn,nymd,nhms,I,p_plume,t,1,d)
               _writeOne(f,'z1_'+bn,nymd,nhms,I,z_plume,t,0,d)
               _writeOne(f,'z2_'+bn,nymd,nhms,I,z_plume,t,1,d)
               _writeOne(f,'k1_'+bn,nymd,nhms,I,k_plume,t,0,d)
               _writeOne(f,'k2_'+bn,nymd,nhms,I,k_plume,t,1,d)
               b += 1

       try:
           f.close()
       except:
           pass

       if self.verb >=1:
           print "[w] Wrote file "+filename

#..............................................................................

#
#                                    Static Methods
#                                    --------------
#

def _writeOne(f,vname,nymd,nhms,I,S,t,k,d):
    """
    Write one sparse variable to a GFIO file.
    """
    A = zeros(d) + __AMISS__
    if I is not None:
        if len(S.shape)==3:
            A[I] = S[t,k,:]
        elif len(S.shape)==1:
            A[I] = S[:]
        else:
            raise ValueError, 'invalid S rank = %d'%len(S.shape)
            
    f.write(vname,nymd,nhms,A)


def getPlumeDeprecated(farea,veg,met,t,ntd,Verb=0):
    
    """
    Runs the Plume Rise extension to compute the extent of the plume.

          p, z, k = getPlume(farea,veg,met,t,ntd)

    where p, z and k are nd-arrays of shape(N,2), N being the
    number of observations (as in met.lon.size). On input,

    farea --- (flaming) fire area
    veg   --- biome type
    
    """

    N = met.lon.size
    nominal_area = 1.e6 # 1 km^2: pixar and farea are in units of km2
    km = met.lev.size
    ptop = met.ptop
    ktop = met.ktop # 1-offset

    if Verb:
        if N>100:
            Np = range(0,N,N/10)
        elif N>10:
            Np = range(0,N,N/10)
        else:
            Np = range(N)
        print ""
        print "                   Plume Rise Estimation for t=%d"%t
        print "                   ------------------------------"
        print ""
        print "  %  |    Lon    Lat  b |   p_bot    p_top  |  z_bot z_top  |  k   k"     
        print "     |    deg    deg    |    mb       mb    |   km     km   | bot top"
        print "---- |  ------ ------ - | -------- -------- | ------ ------ | --- ---"

#   Allocate space
#   --------------
    p_plume = zeros((N,2))
    k_plume = zeros((N,2))
    z_plume = zeros((N,2))
        
#   Compute plume extent, one fire at a time
#   ----------------------------------------
    for i in range(N):

        u = met.fields['u'][i]
        v = met.fields['v'][i]
        T = met.fields['t'][i]
        q = met.fields['qv'][i]
        delp = met.fields['delp'][i]

#       Ensure arrays are top-down as in GEOS-5
#       ---------------------------------------
        if delp[0] > delp[-1]:
            u = u[::-1]
            v = v[::-1]
            T = T[::-1]
            q = q[::-1]
            delp = delp[::-1]

#       Units:
#           farea - km2 (must multiply by nominal area for m2)
#            area  - m2 as required by plume rise model
#       ------------------------------------------------------   
        area = farea[i] * nominal_area
        veg_ = veg[i]

#       Run plume rise model
#       --------------------
        p1, p2, z1, z2, k1, k2, rc = \
            biome(u, v, T, q, delp, ptop, area, veg_)

        k1, k2 = (k1+ktop-1, k2+ktop-1)

        p_plume[i,:] = (p1, p2)
        k_plume[i,:] = (k1, k2)
        z_plume[i,:] = (z1, z2)

        if Verb:
            if i in Np:
                ip = int(0.5+100.*i/N)
                print "%3d%% | %7.2f %6.2f %d | %8.2f %8.2f | %6.2f %6.2f | %3d %3d "%\
                      (ip,met.lon[i],met.lat[i],veg[i], \
                       p2/100,p1/100,z2/1000,z1/1000,k2,k1)

    return (p_plume, z_plume, k_plume)


def _getTagName(tag):
    if tag != None:
        tag_name = tag
    else:    
        if __CVSTAG__ not in (None, ''):
            tag_name = __CVSTAG__
        else:
            tag_name = 'unknown'

    return tag_name

#----
def _getGran(t,mod14_path):
    doy = t.date().toordinal() - date(t.year-1,12,31).toordinal()
    hhmm = "%02d%02d"%(t.hour,5*(t.minute/5))
    patt = mod14_path+'/%4d/%03d/MOD14.A%4d%03d.%s.005.*.hdf'%(t.year,doy,t.year,doy,hhmm)
#    print '--> ', patt
    return glob(patt)

#---

#..............................................................................

def plotMINX(m,imfile=None,Title=None):

    from matplotlib.pyplot import plot, legend, xlabel, ylabel, title, grid, figure, savefig

    I_ = m.z>0
    I = m.z>(m.sample.pblh+250)
    J = m.zm>(m.sample.pblh+250)
    K = m.zt>(m.sample.pblh+250)

    figure(dpi=120)
    plot([0,5],[0,5],'k') # 1:1 line
    plot(m.sample.pblh[I_]/1000,m.z[I_]/1000, 'bo',label='Mode Height < PBL')
    plot(m.sample.pblh[I]/1000,m.z[I]/1000, 'co',label='Mode Height > PBL')
    plot(m.sample.pblh[K]/1000,m.zt[K]/1000, 'ro',label='95%-ile Height > PBL')

    print "Percent above PBL: ",100.*len(m.z[I])/len(m.z[I_])
    
    
    x, ya, yb = m.sample.pblh[K]/1000, m.zt[K]/1000, m.z[K]/1000
    for i in range(len(x)):
        plot([x[i],x[i]],[ya[i],yb[i]],'k')
    
    xlabel('GEOS-5 PBL Height AGL [km]')
    ylabel('MINX Plume Height AGL [km]')
    legend(loc='upper right',fontsize='medium')
    grid()

    if Title is not None:
        title(Title)
        
    if imfile is not None:
        savefig(imfile,bbox_inches='tight')

def plotAreaFRP(m,z_plume,imfile=None,Title=None):
    from matplotlib.pyplot import plot, legend, xlabel, ylabel, title, grid, figure, savefig

    figure(dpi=120)
    J = z_plume>(m.sample.pblh+250)
    area = m.mod14.farea * m.mod14.pixar * 1e6 / 1e4
    plot(area,m.mod14.frp,'bo')
    plot(area[J],m.mod14.frp[J],'ro')
    grid()
    xlabel(r'Fire Area [Ha]')
    ylabel('FRP/Area [MW]')
    
def plotPR(m,z_plume,imfile=None,Title=None):
    from matplotlib.pyplot import plot, legend, xlabel, ylabel, title, grid, figure, savefig
    
    I = m.z>(m.sample.pblh+250)
    J = z_plume>(m.sample.pblh+250)
    
    figure(dpi=120)
    plot([0.5,4],[0.5,4],'k')
    plot(m.sample.pblh/1000,z_plume/1000, 'bo',label='Modeled Below PBL')
    plot(m.sample.pblh[J]/1000,z_plume[J]/1000,'co',label='Modeled above PBL')
    #plot(m.sample.pblh[I]/1000,z_plume[I]/1000,'ro',label='Observed above PBL')
    plot(m.sample.pblh[I]/1000,m.z[I]/1000,'ro',label='Observed above PBL')

    x, ya, yb = m.sample.pblh[I]/1000, z_plume[I]/1000, m.z[I]/1000
    for i in range(len(x)):
        plot([x[i],x[i]],[ya[i],yb[i]],'k')
    
    xlabel('GEOS-5 PBL Height AGL [km]')
    ylabel('Plume Height AGL [km]')
    legend(loc='lower right')
    grid()

    if Title is not None:
        title(Title)
        
    if imfile is not None:
        savefig(imfile,bbox_inches='tight')
    
def scatOpt(m,imfile=None,Title=None,norm=True,xymax=1.9999):

    from matplotlib.pyplot import plot, legend, xlabel, ylabel, title, grid, figure, savefig

    z     = m.z/m.sample.pblh
    z_opt = m.z_opt/m.sample.pblh
    
    I = z>0
    J = z_opt>1
    K = (z>1)&(z_opt>1)

    if not norm:    
        z     = m.z/1000
        z_opt = m.z_opt/1000
        
    figure(dpi=120)
    if norm:
        plot([0,xymax],[0,xymax],'k')
        plot([0,xymax],[1,1],'k')
        plot([1,1],[0,xymax],'k')
    plot(z[I],z_opt[I],'yo',label=r'$z_{OBS}<z_{PBL}$, $z_{OPT}<z_{PBL}$')
    plot(z[J],z_opt[J],'ro',label=r'$z_{OBS}<z_{PBL}$, $z_{OPT}>z_{PBL}$')
    plot(z[K],z_opt[K],'go',label=r'$z_{OBS}>z_{PBL}$, $z_{OPT}>z_{PBL}$')

    if norm:
        xlabel(r'$z_{OBS}/z_{PBL}$')
        ylabel(r'$z_{OPT}/z_{PBL}$')
    else:
        xlabel(r'$z_{OBS}$')
        ylabel(r'$z_{OPT}$')
        
    legend(loc='lower right')
    grid()

    if Title is not None:
        title(Title)
        
    if imfile is not None:
        savefig(imfile,bbox_inches='tight')

def scatOptF(m,imfile=None,Title=None,norm=True,xymax=1.9999):

    from matplotlib.pyplot import plot, legend, xlabel, ylabel, title, grid, figure, savefig, loglog

    z     = m.z/m.sample.pblh
    z_opt = m.z_opt/m.sample.pblh
    
    I = (z>0)&(isnan(m.f_opt)==False)#&(m.f_opt<100)&(m.mod14.frp<100)
    J = I&(z_opt>1)
    K = I&(z>1)&(z_opt>1)

    if not norm:    
        z     = m.z/1000
        z_opt = m.z_opt/1000
        
    figure(dpi=120)
    loglog(m.mod14.frp[I],m.f_opt[I],'yo',label=r'$z_{OBS}<z_{PBL}$, $z_{OPT}<z_{PBL}$')
    loglog(m.mod14.frp[J],m.f_opt[J],'ro',label=r'$z_{OBS}<z_{PBL}$, $z_{OPT}>z_{PBL}$')
    loglog(m.mod14.frp[K],m.f_opt[K],'go',label=r'$z_{OBS}>z_{PBL}$, $z_{OPT}>z_{PBL}$')

    xlabel('FRP [MW]')
    ylabel(r'$\gamma$')     
    #legend(loc='upper right')
    grid()

    if Title is not None:
        title(Title)
        
    if imfile is not None:
        savefig(imfile,bbox_inches='tight')

def scatPR(m,z_plume,imfile=None,Title=None):

    from matplotlib.pyplot import plot, legend, xlabel, ylabel, title, grid, figure, savefig

    z     = m.z/m.sample.pblh
    z_pr = z_plume/m.sample.pblh
    
    I = z>0
    J = (z_pr>1)|((z_pr<1)&(z>1))
    K = (z>1)&(z_pr>1)

    figure(dpi=120)
    plot([0,1.5999],[0,1.5999],'k')
    plot([0,2.999],[1,1],'k')
    plot([1,1],[0,1.5999],'k')
    plot(z[I],z_pr[I],'yo',label=r'$z_{OBS}<z_{PBL}$, $z_{PR}<z_{PBL}$')
    plot(z[J],z_pr[J],'ro',label=r'$z_{OBS}<z_{PBL}$, $z_{PR}>z_{PBL}$')
    plot(z[K],z_pr[K],'go',label=r'$z_{OBS}>z_{PBL}$, $z_{PR}>z_{PBL}$')

    xlabel(r'$z_{OBS}/z_{PBL}$')
    ylabel(r'$z_{PR}/z_{PBL}$')
        
    legend(loc='lower right')
    grid()

    if Title is not None:
        title(Title)
        
    if imfile is not None:
        savefig(imfile,bbox_inches='tight')

#--

def plotKDE(s,s2=None,Title='',imfile=None):

    from matplotlib.pyplot import plot, legend, xlabel, ylabel, title, grid, figure, savefig

    J = (s.pblh>100)&(s.z_a>0)
    r = s.z_a/s.pblh
    if s2 is not None:
        J2 = (s2.pblh>100)&(s2.z_a>0)
        r2 = s2.z_a/s2.pblh
    
    # KDE 1d
    # ------
    figure(dpi=120)
    bins, P = kde.calc_kde1d(r[J],range=(0,5))
    plot(bins,P,'b',linewidth=2,label='Terra')
    if s2 is not None:
        bins, P2 = kde.calc_kde1d(r2[J2],range=(0,5))
        plot(bins,P2,'r',linewidth=2,label='Aqua')
        legend(loc='upper right')
    grid()
    xlabel(r'$z_{PR}/z_{PBL}$')
    ylabel('p.d.f.')
    if Title is not None:
        title(Title)
    if imfile is not None:
        savefig('kde1d.'+imfile,bbox_inches='tight')

    if s2 is not None: return
    
    # KDE 2d
    # ------
    x, y, P = kde.calc_kde2d(s.pblh[J]/1000,s.z_a[J]/1000,x_range=(0,3),y_range=(0,3))
    kde.plot_kde2d(x,y,P,dpi=120,Title=Title,
                   xLabel=r'$z_{PBL}$', yLabel=r'$z_{PR}$')
    grid()
    if imfile is not None:
        savefig('kde2d.'+imfile,bbox_inches='tight')
    
def aveVMD(s):
    """
    Compute average vertical mass distribution, weighted by FRP.
    """
    J = (s.pblh>100)&(s.z_a>0)

    z_c   = s.z_d[J]/1000
    delta = (s.z_f[J] - s.z_d[J])/1000
    pow = s.pow[J]
    pblh = s.pblh[J]/1000
    
    nz = 100
    nf = len(pblh)
    z = linspace(0,5,nz)


    # PBL mass distribution
    # ---------------------
    V = zeros(nz)
    for i in range(nf):
        v = zeros(nz)
        v[z<pblh[i]] = 1
        v = pow[i] * v / sum(v)
        V += v
    v_pbl = V / pow.sum()

    return v_pbl

    v = getvmd(z,z_c,delta)

    

    
#-------------------------------
        
if __name__ == "__main__":

    import crtmmodis_

    m = MINXs_PR('/Users/adasilva/workspace/misrPlumes/western_fires_2013/*.txt')
    m.sampleLoadz('/Users/adasilva/workspace/Data_Analysis/AGU-2014/seac4rs_01.npz')
    m.mod14 = NPZ(('/Users/adasilva/workspace/Data_Analysis/AGU-2014/mod14_seac4rs.npz',))

def _allfires():
    
    print "ok"
    topdir = '/Volumes/ArlindoSD' # SDXC card

    t1 = datetime(2013,8,2)
    t2 = datetime(2013,8,31)

    #for p in ( 'MOD14', 'MYD14'):
    for p in ( 'MYD14',):

        print "Loading fires"
        f = PLUME_L2(None)
        f.restart('%s.fires_nam.%4d-%02d.npz'%(p,t1.year,t1.month))

        f.m = f.qc>0
        
        print "Loading Meteorology"
        f.sampleLoadz('%s/AGU-2014/%s.sample_nam.%4d-%02d.npz'%(topdir,p,t1.year,t1.month))

        I_na = (f.lon>-170)&(f.lon<-50)&(f.lat>15)&(f.lat<80)
        f.sample.I = I_na

        f.getPlume(algo=None,Verbose=True,area_m2=10e4,rad2conv=5)

        f.sample.pblh = f.sample.pblh[I_na]
        savez('%s.sample_a1_r10_r5.2013-08.npz'%p,**f.sample.__dict__)
        
def xxxxx():
    
#--
                
        Files = granules(t1,t2,product=p,rootdir=topdir+'/MODIS')

        print 'Number of files: ', len(Files)
        f = PLUME_L2(Files,Verb=1)

        print "Computing fire properties"
        ###f.classic_var()

        print "Checkpointing"
        f.checkpoint('%s.fires_nam.%4d-%02d.npz'%(p,t1.year,t1.month))

def minx_test():
    
    m = MINXs_PR('/Users/adasilva/workspace/misrPlumes/western_fires_2013/*.txt')
    m.sampleLoadz('/Users/adasilva/workspace/Data_Analysis/AGU-2014/seac4rs_01.npz')
    m.mod14 = NPZ(('/Users/adasilva/workspace/Data_Analysis/AGU-2014/mod14_seac4rs.npz',))
    #plotMINX(m)

    #z_plume = m.getPlume()

    
def arctas():
    # m.getFires(npzFile='mod14_seac4rs.npz')
#    m = MINXs_PR('/Users/adasilva/workspace/misrPlumes/canada2008/Plumes*.txt')

    m = MINXs_PR('/Users/adasilva/workspace/misrPlumes/canada2008/Plumes*.txt')
    m.sampleLoadz('merraero.npz')
    m.mod14 = NPZ('mod14.npz')
    
#    z_plume = m.getPlume()

def Mac():
    m = MINXs_PR('/Users/adasilva/workspace.local/misrPlumes/canada2008/Plumes*.txt')
#    m.sampleLoadz('/Users/adasilva/workspace.local/misrPlumes/canada2008/merra.npz')
#    m.sampleLoadz('merra.npz')
#    m.mod14 = NPZ('mod14.npz')
#    z_plume = m.getPlume()
