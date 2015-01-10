from pyobs import LIDAR_L2, NPZ
from numpy import savez, load, tile, zeros, size, transpose, pi, interp, sqrt
from numpy import cos, sin, arcsin, where, array 
from datetime import datetime

#import VLIDORT_      f77
import VLIDORT_OMI_  #f90
from mie    import getTopo, getAlbedo
from mieobs import *  
from time import clock 
from netCDF4 import Dataset
from gfio   import GFIO, GFIOctl
from math import copysign
MISSING = -1.2676506E+030
pi = 3.1415926

def szangle(c):

    # find the sza given the lat, lon and python time 
    #------------------------------------------------         
   
    a0 = 0.006918
    a1 = 0.399912
    a2 = 0.006758
    a3 = 0.002697
    b1 = 0.070257
    b2 = 0.000907
    b3 = 0.001480

    doy, xhour = [], []

    for t in c.tyme:
        t0 = datetime(t.year,1,1)
        doy.append((t - t0).days +1 )                       # calendar day
        xhour.append(t.hour + t.minute/60.+ t.second/3600.) # TU decimal hour
    doy, xhour = array(doy), array(xhour)
  
    # solar declination in radians
    #------------------------------
    r  = 2.*pi*(doy-1) / 365.
    
    dec=a0 - a1*cos(r) + b1*sin(r) - a2*cos(2.*r) + b2*sin(2.*r)- a3*cos(3.*r) + b3*sin(3.*r)  
    

    # local time in hours
    # --------------------
    timloc = xhour + c.lon[:]/15.    

#    time equation (in mn.dec)
#    -------------------------
    c1=.000075
    c2=.001868
    c3=.032077
    c4=.014615
    c5=.040849    
    et=c1+c2*cos(r)-c3*sin(r)-c4*cos(2.*r)-c5*sin(2.*r)
    et=et*12.*60./pi
    
#   true solar time
#   ---------------- 
    tsv=timloc+et/60.
    tsv=(tsv-12.)
    
#   hour angle
#   ----------- 
    ah=tsv*15.*pi/180.

#   latitude in radians
#   ------------------
    latRad = c.lat[:]*2*pi/360.              

#   Elevation and Solar Zenith and Azimuth Angles in degrees
#   --------------------------------
    x= sin(latRad)*sin(dec)+cos(latRad)*cos(dec)*cos(ah)
    elev=arcsin(x)
         
    az = cos(dec)*sin(ah)/cos(elev)
    s =cos(dec)*sin(ah)/cos(elev)
    for i in range(len(az)) :
        if ((abs(az[i])-1)> 0.00) :
           az[i] = copysign(1.,az[i])
    
    caz=(-cos(latRad)*sin(dec)+sin(latRad)*cos(dec)*cos(ah))/cos(elev)
    azim=arcsin(az)
   
    for i in range(len(caz)) :
        if(caz[i] < 0.) :
           azim[i]=pi-azim[i]
        if(caz[i] > 0. and az[i] < 0.) :
           azim[i]=2*pi+azim[i]
   
    azim=azim+pi
    
    pi2=2*pi
    for i in range(len(azim)) :
        if(azim[i] > pi2) :
           azim[i]=azim[i]-pi2
    
    elev=elev*180./pi        
    asol=90.-elev     # conversion in degree
    phi=180.*azim/pi

    return  asol,phi


def vlidort_scalar (c, channels, sza, albedo, tau, ssa, g, pe, he, te, verbose=0 ):
        """
        
        Uses VLIDORT f90 scalar to compute radiances. On Input,
 
           channels       ---  wavelengths [nm]
             

           solar_zenith   ---  solar zenith angle   
           relat_azymuth  ---  relative azymuth angle         
           sensor_zenith  ---  sensor zenith angle          
        
           albedo   --- surface albedo
           tau      --- aerosol optical depth
           ssa      --- aerosol single scattering albedo
           g        --- aerosol asymmetry factor
           pe       --- pressure at edges [Pa]
           te       --- temperature at edges [K]
           he       --- height at edges above sea-level  [m]
          
        Outputs :
           radiance      ---  radiance (normalized)
           reflectance   ---  reflectance
        """

        nch = len(channels)               # wavelength number
        N = len(c.lon)                    # nobs  
        
        solar_zenith = sza                # solar zenith angle
        relat_azymuth = 60* ones(N)       # relative azimuth angle
#        sensor_zenith = zeros(N)         # sensor zenith angle
        sensor_zenith = 0.03* ones(N)     # sensor zenith angle
                                          # 0.03 deg like CALIPSO (keep the same option in VLIDORT_Lambmod.F)
                                          #(DO_SSCORR_OUTGOING -> True 
                                          # if sensor_zenith = 0 -> DO_SSCORR_OUTGOING -> False 
                                          # and DO_SSCORR_NADIR - > True )
        radiance_, reflectance_, rc = VLIDORT_OMI_.scalar(channels,
                                                     tau, ssa, g, pe, he, te, 
                                                     albedo,solar_zenith,
                                                     relat_azymuth,
                                                     sensor_zenith, MISSING,verbose)
        print 'rad', radiance_
        if rc != 0:
            raise ValueError, "on return from VLIDORT_OMI_.scalar, rc = "+str(rc)

        return (radiance_, reflectance_)

#----
def vlidort_vector (c, channels, sza, albedo, tau, ssa, g,pmom, pe, he, te,I=None, verbose=0 ):
        """

        
        Uses VLIDORT to compute radiances. On Input,
 
           channels       ---  wavelengths [nm]

           solar_zenith   ---  solar zenith angle   
           relat_azymuth  ---  relative azymuth angle         
           sensor_zenith  ---  sensor zenith angle          
        
           albedo         --- surface albedo
           tau      --- aerosol optical depth
           ssa      --- aerosol single scattering albedo
           g        --- aerosol asymmetry factor
           pmom     --- scatter phase matrix
           pe       --- pressure at edges [Pa]
           te       --- temperature at edges [K]
           he       --- height at edges above sea-level  [m]
          
        Outputs :
           radiance       ---  radiance (normalized)
           reflectance    ---  reflectance
           Q              ---  Q stokes vector (normalized)
           U              ---  U stokes vector (normalized)
        """

        nch = len(channels)               # wavelength number
        N = len(I)                        # nobs
        solar_zenith = sza                # solar zenith angle
        relat_azymuth = zeros(N)          # relative azimuth angle (sensor azim - sun azim)
#        sensor_zenith = 0.03* ones(N)    # sensor zenith angle
        sensor_zenith = zeros(N)          # if 0.03 deg like CALIPSO (keep the same option in VLIDORT_Lambmod.F)
                                          #(DO_SSCORR_OUTGOING -> True 
                                          # if sensor_zenith = 0 -> DO_SSCORR_OUTGOING -> False 
                                          # and DO_SSCORR_NADIR - > True )

        radiance_, reflectance_,Q, U, rc = VLIDORT_OMI_.vector(channels,
                                                     tau, ssa, g,pmom, pe, he, te, 
                                                     albedo,solar_zenith,
                                                     relat_azymuth,
                                                     sensor_zenith,MISSING, verbose)
        print 'rad', radiance_, Q, U
        if rc != 0:
            raise ValueError, "on return from VLIDORT_OMI.vector/vector, rc = "+str(rc)

        return (radiance_, Q, U)
#----------------------------
def get_Ocean_Albedo(albe_fn,channels,u10m,v10m):
    """
    Returns ocean surafce albedo fct of wind speed.
    """

# Get precomputed albedo LUT
    # --------------------------
    npz = NPZ(albe_fn)

    # Trimmed wind speed
    # ------------------
    w10m = sqrt(u10m*u10m+v10m*v10m)
    w10m[w10m<0] = 0
    w10m[w10m>50.] = 50.
    print 'len wind', len(w10m)
    # Interpolate albedo
    # ------------------
    print 'len lambda',len(channels)
    albedo = zeros((len(w10m),len(channels)))
    for ch,i in zip(channels,range(len(channels))):
        j = list(npz.channels).index(ch)
        albedo[:,i] = interp(w10m,npz.speed,npz.albedo[:,j])

    return albedo

#----------------------------
def get_IGBP(lat,lon,Path):
    """
    Returns IGBP surface type.
    """
    from pyobs import igbp
        
    vtype = igbp.getDetailedVeg(lon,lat,Path)
        
    return vtype

#----------------------------
def get_ALBEDO(lat,lon,time,vtype,channels,month,path):
        
        """
        Returns albedo surface for 355, 532, 1064 nm.
        """
        albe_fn = 'ocean_albedo_lidar.npz'
        wind_fn = '/nobackup/MERRA/opendap/slv_Nx.ctl'
       

        nch = len(channels)
        nobs = len(lat)
        albedo=ones((nobs))
        
        # Albedo over ocean fct wind speed
        #----------------
        w = GFIOctl(wind_fn).sampleVars(lon,lat,time,
                                              onlyVars=['U10M','V10M'],Verbose=True)
              
        albedo_ocean = get_Ocean_Albedo(albe_fn,channels,w.U10M,w.V10M)
        print albedo_ocean.shape

        # Albedo for each wavelength (355, 532, 1064 nm)
        #------------------
        for ch in channels :

          if (ch == 355.) :
              print 'ok channels'
              TOMS = load(path+'TOMS_clim_360nm.npz') # shape lat, lon, month 
              
              for l in range(nobs): 
                 if  lon[l] < 179.5: 
                   i=where(TOMS['lats']<lat[l])
                   j=where(TOMS['lons']>lon[l])                 
                   ilat=i[0][0]
                   ilon=j[0][0]
                 else:
                   i=where(TOMS['lats']<lat[l])
                   ilat=i[0][0]
                   ilon=359
                                   
                 if(abs(lat[l]-TOMS['lats'][ilat])<0.5):
                   ilat_=ilat
                 else: 
                    ilat_=ilat - 1
                 if(abs(lon[l]-TOMS['lons'][ilon])<0.5):
                    ilon_=ilon
                 else: 
                    ilon_=ilon - 1
                 albedo[l] = TOMS['ref360all'][ilat_,ilon_,month]/100.
                  
          if (ch == 532.) :
              albedo_532 = ones(nobs)

              albedo_MODIS= load(path+'albedo_MODIS_555nm_20090715.npz')['albedo']
              
              for l in range(nobs): 
                  albedo_532[l] =  albedo_MODIS[l]                    
                  if (vtype[l] ==17) :  # if water replace undef by ocean albedo
                      albedo_532[l] = albedo_ocean[l]
              undef = where(albedo_532 < 0)
              for x in undef[0]:  # 1) replace undef by albedo of obs+-1 if same type 
                  if (vtype[x] == vtype[x-1]):
                      if (albedo_532[x-1] > 0) :
                          albedo_532[x] = albedo_532[x-1]
                  elif (vtype[x]  == vtype[x+1]) :
                      if (albedo_532[x+1] > 0) :
                          albedo_532[x] = albedo_532[x+1]
              undef = where(albedo_532 < 0)        
              for x in undef[0]:  # 2) replace undef by albedo value over antartica
                  if (vtype[x]==15) and lat[x] < -65.:
                      albedo_532[x] = 0.9 
              undef = where(albedo_532 < 0)  # 3) points that stay
              for x in undef[0]: 
                   if (vtype[x]==15) :
                      albedo_532[x] = 0.8
                   else :
                      albedo_532[x] = 0.07
              albedo[:] = albedo_532



          if (ch == 1064.) :

              albedo_1064 = -999*ones(nobs)
              MODIS858 = load(path+'albedo_MODIS_858nm_20090715.npz')['albedo']
              MODIS1240 = load(path+'albedo_MODIS_1240nm_20090715.npz')['albedo']
              def_858 = where(MODIS858 > 0)
              def_1240 = where(MODIS1240 > 0)
              MODIS_858 = MODIS858[def_858[0]]
              MODIS_1240 = MODIS1240[def_858[0]]
              a = (MODIS_1240 - MODIS_858)/(1240. - 858.)
              b = MODIS_1240 - a * 1240.
              MODIS_1064 = a *1064. +b

              albedo_1064[def_858[0]] =  MODIS_1064 
             
              for l in range(nobs):
                   if (vtype[l] ==17) :   # if water replace undef by ocean albedo at 1064 nm
                       albedo_1064[l] = albedo_ocean[l]
              undef = where(albedo_1064 < 0)
              for x in undef[0]:  # 1) replace undef by albedo of obs+-1 if same type 
                  if (vtype[x] == vtype[x-1]):
                      if (albedo_1064[x-1] > 0) :
                          albedo_1064[x] = albedo_1064[x-1]
                  elif (vtype[x]  == vtype[x+1]) :
                      if (albedo_1064[x+1] > 0) :
                          albedo_1064[x] = albedo_1064[x+1]
              undef = where(albedo_1064 < 0)        
              for x in undef[0]:  # 2) replace undef by albedo value over antartica
                  if (vtype[x]==15) and lat[x] < -65.:
                      albedo_1064[x] =  0.60  # albedo = S. Wuttke et al.,2006 
              undef = where(albedo_1064 < 0)  # 3) points that stay
              for x in undef[0]: 
                  if (vtype[x]==15) :
                     albedo_1064[x] = 0.60
                  else :
                     albedo_1064[x] = 0.05
              albedo[:] = albedo_1064

        return albedo

#-------------------------------------------------------------------------------
if __name__ == "__main__":

    channels = [355.]
    month = [07]    

    topo_fn  = '/nobackup/1/VLIDORT/topo/topography.1152x721.nc'
    albedo_dir = '/nobackup/2/vbuchard/LIDAR/ALBEDO/'
    
    c_dn = '/nobackup/2/vbuchard/LIDAR/'         # Pete's Calipso interp
    a_dn = '/nobackup/2/vbuchard/LIDAR/'         # aer_Nv + interp

#    c_fn = 'dR_Fortuna-2-4-b4/dR_Fortuna-2-4-b4.calipso_532nm.20090715.nc' # Pete's interp
#    a_fn = 'dR_Fortuna-2-4-b4/dR_Fortuna-2-4-b4.calipso_aer.20090715.npz' 
#    c_fn = 'dR_Fortuna-M-1-1/orig/dR_Fortuna-M-1-1.glas_532nm.20031009.nc' # Pete's interp
#    a_fn = 'dR_Fortuna-M-1-1/dR_Fortuna-M-1-1.glas_aer.20031009.npz'       # aer collocation

    c_fn = 'dR_MERRA-AA-r2/dR_MERRA-AA-r2.calipso_355nm.20090715.nc' # Pete's interp
    a_fn = 'dR_MERRA-AA-r2/dR_MERRA-AA-r2.calipso_aer.20090715.npz'  # aer collocation
    
    # Read relevant data
    # ------------------
    c = LIDAR_L2(c_dn+c_fn)   # Pete's interp with CALIPSO coordinates
    a = NPZ(a_dn+a_fn)        # aer_v interpolated to obs location
    aerToUpper(a)
    km = a.delp.shape[1]
    
    # Get IGBP surface type
    # ---------------------
    IGBP_dir = '/nobackup/Emissions/Vegetation/GL_IGBP_INPE'  # on calculon
    vtype = get_IGBP(c.lat,c.lon,IGBP_dir)
    
    # Get albedo at lat lon
    # ------------------------
    alb_dir = '/nobackup/2/vbuchard/LIDAR/ALBEDO/'  
    albedo=get_ALBEDO(c.lat,c.lon,c.tyme,vtype,channels,month,alb_dir)
  
    # Get the VLIDORT inputs scalar or vector
    # ---------------------------------------
    zs = getTopo(topo_fn,c.lon,c.lat)
    pe, ze, te = getEdgeVars(a)
    
#    tau, ssa, g = getAOPscalar(a,channels,rcfile='Aod_EOS.rc')
  
    # Height above sea level
    # ---------------------- 
    he = ze + tile(zs,(km+1,1))
       
    # Get the Solar Zenith and Azimuth Angles
    # -------------------------
    asol, phi = szangle(c)
    
    # call VLIDORT scalar or vector (if vector -> check mie tables in Aod_EOS.rc)
    #-------------------------------
    nch = len(channels)               # wavelength number
    N = c.lon.shape[0]
#    N = 1000
    rad_ = ones((N,nch)) 
    Q_ = ones((N,nch)) 
    U_ = ones((N,nch)) 
    nMom = 300
    mobs = 200
    for i in range(0,N,mobs):
        Irange = range(i,min(i+mobs,N))
        
        tau_, ssa_, g_, pmom_ = getAOPvector(a,channels, I= Irange,nMom=nMom,rcfile='Aod_EOS.rc')
        asol_=asol[Irange]    
        albedo_=albedo[Irange]  
        pe_=pe[:,Irange].copy()
        he_=he[:,Irange].copy()
        te_=te[:,Irange].copy()
        rad, Q, U = vlidort_vector(c,channels,asol_,albedo_,tau_, ssa_, g_,pmom_, pe_, he_, te_,I=Irange)
        
        rad_[Irange,:] = rad
        Q_[Irange,:] = Q
        U_[Irange,:]= U
    print 'save in a file'
#    savez(c_dn+'radiances_1064nm_20090715_vector_U_Q_0.npz',\
#          lon=c.lon,lat=c.lat,rads=rad_,Q=Q_,U=U_,sza=asol,phi=phi)
#def HOLD():
    # Add the radiances and SZA values in the same netcdf file 
    #--------------------------------------------------------
    file=Dataset(c_dn+c_fn,'a',format='NETCDF4')    
    sza=file.createVariable('sza','f4','nt')
    azim=file.createVariable('azimuth','f4','nt')
    radiances=file.createVariable('radiances','f4','nt')
    Q=file.createVariable('Q_STOKES','f4','nt')
    U=file.createVariable('U_STOKES','f4','nt')
    IGBP=file.createVariable('IGBP_type','i4','nt')
    Alb=file.createVariable('Surf_Alb','f4','nt')
    radiances.desciption = 'TOA sun-normalized radiances'
    sza.desciption = 'Solar Zenith Angle'
    azim.desciption = 'Solar azimuth Angle'
    Q.desciption = 'Q Stokes vector'
    U.desciption = 'U Stokes vector'
    IGBP.description = 'IGBP surface type'
    Alb.description = 'Surface Albedo'
    radiances[:]=rad_
    Alb[:] = albedo
    sza[:]=asol
    azim[:]=phi
    IGBP[:]=vtype
    Q[:]=Q_
    U[:]=U_
    file.close()


