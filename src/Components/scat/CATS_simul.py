#!/usr/bin/env python

from numpy import savez, load, tile,zeros,interp, transpose
from numpy import loadtxt, ones, array
from datetime import datetime
from netCDF4 import Dataset

import VLIDORT_
import VLIDORT_OMI_
import pdb
from math import *
from time import clock 

MISSING = -1.2676506E+030
pi = 3.1415926

SDS_dR = ['h', 't', 'latitude','longitude','ps','delp']

def vlidort_scalar (channels, sza, albedo, tau, ssa, g, pe, he, te, verbose=0 ):
        """
        
        Uses VLIDORT to compute radiances. On Input,
 
           channels       ---  wavelengths [nm]
           ps             ---  surface pressure in Pascal   

           solar_zenith   ---  solar zenith angle   
           relat_azymuth  ---  relative azymuth angle         
           sensor_zenith  ---  sensor zenith angle          
        
           albedo         --- surface albedo
           tau      --- aerosol optical depth
           ssa      --- aerosol single scattering albedo
           g        --- aerosol asymmetry factor
           pe       --- pressure at edges [Pa]
           te       --- temperature at edges [K]
           he       --- height at edges above sea-level  [m]
          
        Outputs :
           radiance       ---  radiance (normalized)
           ai             ---  aerosol index
           reflectivity   ---  reflectivity
        """

        nch = len(channels)               # wavelength number
#                                          # nobs  
        N=len(sza)
        solar_zenith = sza                # solar zenith angle
        relat_azymuth = zeros(N)          # relative azimuth angle
        sensor_zenith = zeros(N)         # sensor zenith angle
#        sensor_zenith = 0.03 * ones(N)    # 0.03 deg like CALIPSO (keep the same option (DO_SSCORR_OUTGOING -> True in VLIDORT_Lambmod.F if sensor_zenith = 0 -> DO_SSCORR_OUTGOING -> False and DO_SSCORR_NADIR - > True )

        print ' essai', channels, sza.shape,albedo.shape,tau.shape,ssa.shape, g.shape,pe.shape,he.shape, te.shape  

        radiance_, reflectance_, rc = VLIDORT_OMI_.scalar(channels,
                                                     tau, ssa, g, pe, he, te, 
                                                     albedo,solar_zenith,
                                                     relat_azymuth,
                                                     sensor_zenith, MISSING,verbose)
#        print 'rad', radiance_
        if rc != 0:
            raise ValueError, "on return from VLIDORT_OMI_.scalar, rc = "+str(rc)

        return (radiance_, reflectance_)

#---------------------------------------------------------------------
def radiances_scalar(ssa_du,ssa_ss,g_du,g_ss,albedo_):


  # Inputs to run the code (wavelength, cloud phase fct file name)
  # ----------------------------
 
    channels = [532.,]
 
    cloud_fn = '/home/vbuchard/workspace/Phase_function_cloud/P11_MODIS_ice_de060_560nm_170mom.dat'   

    filename1='/nobackup/2/vbuchard/LIDAR/dR_Fortuna-2-4-b4/dR_Fortuna-2-4-b4.calipso_532nm.20090715.nc'
      
  # Read variables delp and P from filename1
  # ---------------------
    f1 = Dataset(filename1)
    self={}
    for sds_dr in SDS_dR:
       v = f1.variables[sds_dr]
       self[sds_dr] = v

    delp=self['delp'][1000] 

    P = ones(72)
    P[71] = self['ps'][1000]
    for K in range (70,-1,-1):                   # Pressure for each layer
         P[K] = P[K+1]- 0.5*(delp[K]+delp[K+1])  # Pressure in Pa
    P = P * 0.01                                 # Pressure in hPa
    P_=P[::-1]

    z=self['h'][1000]  # at lat in tropics
    t=self['t'][1000]  

    z_=z[::-1]         # array croissant
    t_=t[::-1]
    new_z = range(0,20500,500)
          
    t_interp = ones(41)
    p_interp = ones(41)
    te_interp = interp(new_z,z_,t_)               # interpolation T in the grid 0-20 km (500m grid)
    pe_interp = interp(new_z,z_,P_)               # interpolation P in the grid 0-20 km (500m grid)
     
    pe=array(pe_interp[::-1])
    te=array(te_interp[::-1])   
    he=array(new_z[::-1])

    pe_=ones((20,41))
    pe_[:,:]=pe
    pe_=transpose(pe_)
    te_=ones((20,41))
    te_[:,:]=te
    te_=transpose(te_)
    he_=ones((20,41))
    he_[:,:]=he
    he_=transpose(he_)             
 
    km= 40   # km = 40 and nobs = 20
    
#    albedo = ones((20,1))
#    albedo[:,0]=0.07
    albedo = ones((20,1))
    albedo[:,0]=albedo_ 
 

    # Get the VLIDORT inputs
    # ----------------------

    sza=array(20*[20.])   # nobs = 20 and SZA = 20 deg

    tau=zeros((40,1,20))
    tau[1:5,0,0]=0.0005
    tau[1:5,0,1]=0.001
    tau[1:5,0,2]=0.002
    tau[1:5,0,3]=0.005
    tau[1:5,0,4]=0.01
    tau[1:5,0,5]=0.02
    tau[1:5,0,6]=0.05
    tau[1:5,0,7]=0.1
    tau[1:5,0,8]=0.2
    tau[1:5,0,9]=0.5
    tau[20:24,0,10]=0.0005
    tau[20:24,0,11]=0.001
    tau[20:24,0,12]=0.002
    tau[20:24,0,13]=0.005
    tau[20:24,0,14]=0.01
    tau[20:24,0,15]=0.02
    tau[20:24,0,16]=0.05
    tau[20:24,0,17]=0.1
    tau[20:24,0,18]=0.2
    tau[20:24,0,19]=0.5

    ssa = zeros((40,1,20))
    ssa[1:5,0,0:10]= ssa_ss       # ssa sea salt 
    ssa[20:24,0,10:20]= ssa_du    # ssa dust 
     
    g = zeros((40,1,20))
    g[1:5,0,0:10]= g_ss       # g sea salt 
    g[20:24,0,10:20]=g_du     # g dust 
     
    tau_=tau[::-1]
    ssa_=ssa[::-1]
    g_=g[::-1]

    # read cloud phase fction P11
    # -----------------------------
    
    cdsmom = loadtxt(cloud_fn,skiprows=2)  # 170 moments
    nMom_cd = len(cdsmom)

#    pdb.set_trace()
    # call vlidort
    radiance, reflectance=vlidort_scalar(channels, sza, albedo, tau_, ssa_, g_, pe_, he_, te_)
    
#    t = clock() - t0
#    print 'time', t
#    print 'radiance', radiance
    savez('radiances_532nm_CATS_simul.npz',sza=sza,tau=tau_,ssa=ssa_,g=g_,rads=radiance,albedo=albedo)
    return (radiance)


#-------------------------------------------------------------------------------
if __name__ == "__main__":

#    SZA=[20.]

    ALBEDO = [0.07,0.1,0.15,0.35]

    TAU=[0.0005,0.001,0.002,0.005,0.01,0.02,0.05,0.1,0.2,0.5]

#    SSA_DU=[0.85,0.90,0.95]
    SSA_DU=[1.,1.,1.]  
    SSA_SS=[1.]     

#    G_DU=[0.70,0.75]
#    G_SS=[0.80,0.85]
    G_DU=[0.60,0.8]
    G_SS=[0.60,0.8]
    rad_ss=ones((4,10,1,2))
    rad_dust=ones((4,10,3,2))

    for j in range(4):
        albedo_=ALBEDO[j]
        for i in range(len(SSA_DU)):
            ssa_du=SSA_DU[i]
            ssa_ss=1.                 # same value for all
            for k in range(len(G_DU)):
                g_du=G_DU[k]
                g_ss=G_SS[k]
                radiance=radiances_scalar(ssa_du,ssa_ss,g_du,g_ss,albedo_)
                rad_ss[j,0:10,0,k]=radiance[0:10,0]
                rad_dust[j,0:10,i,k]=radiance[10:20,0]

# write in a netCDF file

    file=Dataset('TEST_2_6/Solar_background_seasalt_dust_sza_20deg_f90_old_2.5_test_alt_aerosols_g.nc4','w',format='NETCDF4')
    file.createDimension('albedo',4)
#    file.createDimension('ssa_dust',3)
    file.createDimension('ssa_dust',3)
    file.createDimension('ssa_ss',1)
    file.createDimension('g_dust',2)
    file.createDimension('g_ss',2)
    file.createDimension('AOD',10)
    
    albedo=file.createVariable('surface_albedo','f8',('albedo',))
    ssa_dust=file.createVariable('ssa_dust','f8',('ssa_dust',))
    ssa_ss=file.createVariable('ssa_ss','f8',('ssa_ss',))
    g_dust=file.createVariable('g_dust','f8',('g_dust',))
    g_ss=file.createVariable('g_ss','f8',('g_ss',))
    AOD=file.createVariable('AOD','f8',('AOD',))
    
    radiances_ss=file.createVariable('radiances_ss','f8',('albedo','AOD','ssa_ss','g_ss',))
    radiances_dust=file.createVariable('radiances_dust','f8',('albedo','AOD','ssa_dust','g_dust',))
    radiances_ss.desciption = ' TOA radiances at nadir for the 10 profiles with a sea salt layer between 0.5-2.5 km at SZA = 20deg'
    radiances_dust.desciption = ' TOA radiances at nadir for the 10 profiles with a dust layer between 10-12 km at SZA = 20deg'
    albedo.desciption='Surface albedo'
    ssa_dust.desciption='Single Scattering Albedo for Dust'
    ssa_ss.desciption='Single Scattering Albedo for Sea Salt'
    g_dust.desciption='Asymetry factor for Dust'
    g_ss.desciption='Asymetry factor for Sea Salt'
    AOD.desciption='Aerosol Optical Depth for a 0.5km layer'
    albedo[:]=ALBEDO
    ssa_dust[:]=SSA_DU
    ssa_ss[:]=1.
    g_dust[:]=G_DU
    g_ss[:]=G_SS
    AOD[:]=TAU
    radiances_ss[:]=rad_ss
    radiances_dust[:]=rad_dust

    file.close



