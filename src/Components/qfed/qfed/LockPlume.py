"""
A Python interface to the simpler Lock Plume model.

"""

from numpy import zeros, ones, meshgrid, linspace, any, array, dot, arange, \
                  pi, sin, cos, arccos
from scipy import optimize as opt

from datetime import date, timedelta
from glob     import glob

from MAPL.constants import *
from LockPlume_     import plume, plumeg  # f2py extension
from pyobs.minx     import MINXs

from pyobs import NPZ

import eta

__VERSION__ = 1.0
__CVSTAG__  = '@CVSTAG'

#----
class MINXs_PR(MINXs):

    """
    Extension of the MINXs class adding the Lock Plume Rise
    functionality. This class handles non-gridded, observation
    location fires.
    """
    
    def getPlume1(self,i,rad2conv=5,bfac=15.,ffac=None,hfac=1):
 
        """
        Runs the Plume Rise extension to compute the extent of the plume.
        On input,
           i     ---  index to operate on
           bfac  ---  factor to increase bstar
           hfac  ---  factor to increase fire heat flux
           
        Recall:
            bstar = (g/(rho*ustar)) * (hflx+(1-eps)*eflx)/(cp * Ta)
	"""

        ptop = 1. # Pascal
        ustar = self.sample.ustar[i]

        # Convective Heat Flux
        # --------------------
        area = 1e6 * self.mod14.farea[i]           # km2 --> m2
        frp_MW = self.mod14.frp[i]                 # MW
        hflux_W = 1e6 * rad2conv * frp_MW / area  # kW/m2
        
        # Derive (scaled) fire Heat Flux contribution
        # -------------------------------------------
        bfire = (MAPL_GRAV / (self.sample.rhoa[i]*ustar)) \
              * (hflux_W / (MAPL_CP * self.sample.tsh[i]))
        bstar = bfac * self.sample.bstar[i] + hfac * bfire
        
        u = self.sample.u[i]
        v = self.sample.v[i]
        T = self.sample.t[i]
        q = self.sample.qv[i]
        delp = self.sample.delp[i]

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
        return plume(u,v,T,q,delp,ptop,bstar,ustar)

    def getPlume(self,I=None,**kwopts):
 
        """
        Runs the Plume Rise extension to compute the extent of the plume.
	"""

        z_lock = -99. * ones(self.N)

        R = arange(self.N)
        if I is not None:
            R = R[I]

#       Loop over time
#       --------------
        for i in R:
            z_lock[i] = self.getPlume1(i,**kwopts)
    
        return array(z_lock)

#---
    def getOptB(self,Verbose=True):
        """
        Find optimal bstar factor to match observed plume height.
        """

        if Verbose:
            print ""
            print "    Scaled bstar Optimization"
            print ""
            print "Plume | b_opt  |    J     |  Ni | Nf"
            print "------|--------|----------|-----|-----"

        self.b_opt = ones(self.N)
        self.z_opt = ones(self.N)
        for i in range(self.N):
            xmin, fval, iter, fcalls = opt.brent(CostFuncB,args=(self,i),brack=(1.,5.),full_output=True)
            bfac = xmin**2
            if Verbose:
               print "%5d | %6.2f | %8.2f | %3d | %3d "%(i, bfac, fval, iter, fcalls)
            self.b_opt[i] = bfac
            self.z_opt[i] = self.getPlume1(i,bfac=self.b_opt[i])

        if Verbose:
            print "------|--------|----------|-----|-----"

#---
    def getOptF(self,Verbose=True):
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
            if self.frp[i]>0:
                xmin, fval, iter, fcalls = opt.brent(CostFuncF,args=(self,i),brack=(0.,1.),full_output=True)
                ffac = xmin**2
                self.f_opt[i] = ffac
                self.z_opt[i] = self.getPlume1(i,ffac=self.f_opt[i])
                if Verbose:
                    print "%5d | %6.2f | %8.2f | %3d | %3d "%(i, ffac, fval, iter, fcalls)
            else:
                self.f_opt[i] = -99.
                self.z_opt[i] = -99.

        if Verbose:
            print "------|--------|----------|-----|-----"

#---
    def getOptGfrp(self,AbovePBL=False,Verbose=True):
        """
        Find optimal fire modified bstar to match observed plume height.
        """

        if Verbose:
            print ""
            print "    FRP Modified bstar Global Optimization"
            print ""
            print "Plume | f_opt  |    J     |  Ni | Nf"
            print "------|--------|----------|-----|-----"

        self.ffac_opt = -99. * ones(self.N)
        self.zfrp_opt = -99. * ones(self.N)

        if AbovePBL:
            I = (self.frp>0)&(m.z>m.sample.pblh)
        else:
            I = self.frp>0

        xmin, fval, iter, fcalls = opt.brent(CostFuncG,args=(self,I),brack=(0.,1.),full_output=True)
        ffac = xmin**2
        self.ffac_opt = ffac
        self.zfrp_opt = self.getPlume(ffac=ffac,I=I)
        if Verbose:
            print "%5d | %6.2f | %8.2f | %3d | %3d "%(0, ffac, fval, iter, fcalls)

        if Verbose:
            print "------|--------|----------|-----|-----"

#---
    def getOptGhf(self,AbovePBL=False,dtol=30.,Verbose=True):
        """
        Find optimal fire modified bstar (by HF) to match observed plume height.

        dtol = distance (in km) tolerance to detected fire.
        
        """

        if Verbose:
            print ""
            print "    Fire Heat Flux Modified bstar Global Optimization"
            print ""
            print "Plume | f_opt  |    J     |  Ni | Nf"
            print "------|--------|----------|-----|-----"

        self.hfac_opt = -99. * ones(self.N)
        self.zhf_opt = -99. * ones(self.N)

        if AbovePBL:
            I = (self.fdist<dtol)&(m.z>m.sample.pblh)
        else:
            I = (self.fdist<dtol)

        xmin, fval, iter, fcalls = opt.brent(CostFuncGhf,args=(self,I),brack=(0.,1.),full_output=True)
        hfac = xmin**2
        self.hfac_opt = hfac
        self.zhf_opt = self.getPlume(hfac=hfac,I=I)
        if Verbose:
            print "%5d | %6.2f | %8.2f | %3d | %3d "%(0, hfac, fval, iter, fcalls)

        if Verbose:
            print "------|--------|----------|-----|-----"

#----
    def getFires(self,mod14_path='/home/adasilva/iesa/aerosol/data/MODIS/Level2/MOD14',
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


    def getFiresOld(self,mod14_path='/nobackup/MODIS/Level2/MOD14', Verbose=True):
        """
        Retrieves Level2 MOD14 fire data for each MINX fire.
        """
        from dozier import DOZIER

        self.mod14 = []
        d2r = pi / 180.
        a = MAPL_RADIUS/1000. # Earth radius in km
        dt = timedelta(seconds=5*60)

        self.hflux = -99. * ones(self.N) # fire heat flux estimate
        self.fdist = -99. * ones(self.N) # distance to detected fire
        self.mfrp = -99. * ones(self.N)  # nearest MODIS FRP
        self.pixar = -99. * ones(self.N)  # nearest MODIS FRP

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
            x0 = cos(d2r*self.lat[i]) * cos(d2r*self.lon[i])
            y0 = cos(d2r*self.lat[i]) * sin(d2r*self.lon[i])
            z0 = sin(d2r*self.lat[i])
            dx = x0*cos(d2r*m.lat) * cos(d2r*m.lon)
            dy = y0*cos(d2r*m.lat) * sin(d2r*m.lon)
            dz = z0*sin(d2r*m.lat) 
            s  = a * arccos(dx+dy+dz) # great circle distance
            j = s.argmin()

            # Estimate fire heat flux
            # -----------------------
            m.classic_var() # classic Dozier
            self.hflux[i] = m.hflux[j] # radiative heat flux
            self.fdist[i] = s[j]
            self.mfrp[i] = m.pow[j]
            self.pixar[i] = m.pixar[j]
            if Verbose:
                print " ", t, "| %8.2f | %8.2f %8.2f | %8.2f"%(s[j], self.frp[i], m.pow[j], m.hflux[j] )

        if Verbose:
            print "  --------------------|----------|-------------------|----------"
#----
def _getGran(t,mod14_path):
    doy = t.date().toordinal() - date(t.year-1,12,31).toordinal()
    hhmm = "%2d%02d"%(t.hour,5*(t.minute/5))
    patt = mod14_path+'/%4d/%03d/MOD14.A%4d%03d.%s.005.*.hdf'%(t.year,doy,t.year,doy,hhmm)
    return glob(patt)

#----

def CostFuncB(b,m,i):
    z_lock = m.getPlume1(i,bfac=b**2)
    return (m.z[i]-z_lock)**2

#----
def CostFuncF(f,m,i):
    z_lock = m.getPlume1(i,ffac=f**2)
    return (m.z[i]-z_lock)**2

#----
def CostFuncG(f,m,I):
    z_lock = m.getPlume(ffac=f**2,I=I)
    d = m.z[I] - z_lock[I]
    return dot(d,d)

#----
def CostFuncGhf(f,m,I):
    z_lock = m.getPlume(hfac=f**2,I=I)
    d = m.z[I] - z_lock[I]
    return dot(d,d)

#..............................................................................
class Gridded_PR(object):

    def __init__(self,lon_range=None,lat_range=None,Verbose=True):
        self.lon_range = lon_range
        self.lat_range = lat_range
        self.Verbose = Verbose

 #---
    def _setdim(self,ga,gatime,fh=None):
        ga('set time '+gatime)
        if self.lon_range is not None:
            ga('set lon %f %f'%self.lon_range)
        if fh is not None:
            ga('set x 1 %d'%fh.nx)
        if self.lat_range is not None:
            ga('set lat %f %f'%self.lat_range)
 #---
    def getMERRA(self,gatime):
        """
        Get Met fields from MERRA.
        """
        from grads import GrADS
        ga = GrADS(Window=False,Echo=False)

        fh = ga.open('http://goldsmr2.sci.gsfc.nasa.gov:80/dods/MAT1NXFLX') # 2D fluxes
        self._setdim(ga,gatime)
        for q in ('ustar', 'bstar', 'pblh', 'hflux', 'eflux', 'rhoa', 'tsh'):
            x = ga.exp('re(%s,1.25,1)'%q) 
            self.__dict__[q] = x.data[:,:]
            if self.Verbose:
                print '          [2d] Read ', q, x.shape
        ga('close 1')

        fh = ga.open('http://goldsmr1.sci.gsfc.nasa.gov:80/dods/MAT3FVCHM') # 3D fields
        self._setdim(ga,gatime,fh)
        ga('set z 1 72')
        for q in ('u','v','t','qv','delp'):
            x = ga.exp(q)
            self.__dict__[q] = x.data[:,:,:]
            if self.Verbose:
                print '          <3d> Read ', q, x.shape
        ga('close 1')
        
 #---
    def getPlume1(self,j,i,bfac=15.):
 
        """
        Runs the Plume Rise extension to compute the extent of the plume, for
        a single (lat,lon).

        On input,
           i     ---  longitude index to operate on
           j     ---  latitude index to operate on
           bfac  ---  factor to scale bstar

        Recall:
            bstar = (g/(rho*ustar)) * (hflx+(1-eps)*eflx)/(cp * Ta)
	"""

        ptop = 1. # Pascal
        ustar = self.ustar[j,i]

        # Simply scale model bstar
        # ------------------------
        bstar = bfac * self.bstar[j,i]
        
        u = self.u[:,j,i]
        v = self.v[:,j,i]
        T = self.t[:,j,i]
        q = self.qv[:,j,i]
        delp = self.delp[:,j,i]

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
        return plume(u,v,T,q,delp,ptop,bstar,ustar)

 #---

    def getPlume(self,j,i,bfac=15.):
 
        """
        Runs the Plume Rise extension to compute the extent of the plume.
        This works on gridded arrays.

        On input,
           bfac  ---  factor to scale bstar

        Recall:
            bstar = (g/(rho*ustar)) * (hflx+(1-eps)*eflx)/(cp * Ta)
	"""

        ptop = 1. # Pascal
        ustar = self.ustar

        # Simply scale model bstar
        # ------------------------
        bstar = bfac * self.bstar
        
        u = self.u
        v = self.v
        T = self.t
        q = self.qv
        delp = self.delp

        # Ensure arrays are top-down as in GEOS-5
        # ---------------------------------------
        if delp[0,0,0] > delp[-1,0,0]:
            u = u[::-1,:,:]
            v = v[::-1,:,:]
            T = T[::-1,:,:]
            q = q[::-1,:,:]
            delp = delp[::-1,:,:]

        # Run plume rise model
        # --------------------
        return plumeg(u,v,T,q,delp,ptop,bstar,ustar)

#..............................................................................

if __name__ == "__main__":

    print "ok"
    
    #def hold():
    n = MINXs_PR('/home/adasilva/workspace/misrPlumes/western_fires_2013/*.txt')
    n.sampleLoadz('/home/adasilva/workspace/Data_Analysis/AGU-2014/seac4rs_01.npz')
    n.mod14 = NPZ(('/home/adasilva/workspace/Data_Analysis/AGU-2014/mod14_seac4rs.npz',))
    z_plume = n.getPlume()

def hold():
    
#    m = MINXs_PR('/Users/adasilva/workspace.local/misrPlumes/canada2008/Plumes_O450*.txt')
#    m.sampleLoadz('/Users/adasilva/workspace.local/misrPlumes/canada2008/merra_O450.npz')
    m = MINXs_PR('/Users/adasilva/workspace.local/misrPlumes/canada2008/Plumes*.txt')
    m.sampleLoadz('/Users/adasilva/workspace.local/misrPlumes/canada2008/merra.npz')
#    m.getOptB()    
#    m.getOptF()    
#    m.getOptG()    

