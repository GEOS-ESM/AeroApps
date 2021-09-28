#!/usr/bin/env python

"""
    Does some simple calculations to test trace gas  PACE calculations

    Adapted from benchmark4Amir.py
    Patricia Castellanos, April 2020

"""

import os
import sys
from   netCDF4 import Dataset
from   netCDF4 import Dataset as ncread
import numpy   as np
from MAPL.constants import *
from py_leo_vlidort import VLIDORT_POLAR_
from scipy.interpolate import interp1d
import scipy.integrate as integrate
from pyhdf.SD import SD, SDC
from multiprocessing import Pool

format   = 'NETCDF4_CLASSIC'
plane_parallel = True
MISSING = -1.e+20


WrapperFuncs = {'MODIS_BRDF'                : VLIDORT_POLAR_.vector_brdf_modis,
                'MODIS_BRDF_BPDF'           : VLIDORT_POLAR_.vector_brdf_modis_bpdf,
                'BPDF'                      : VLIDORT_POLAR_.vector_bpdf,
                'LAMBERTIAN'                : VLIDORT_POLAR_.vector_lambert,
                'LAMBERTIAN_BPDF'           : VLIDORT_POLAR_.vector_lambert_bpdf,
                'GissCX'                    : VLIDORT_POLAR_.vector_gisscx,
                'CX'                        : VLIDORT_POLAR_.vector_cx,
                'OCICX'                     : VLIDORT_POLAR_.vector_ocicx,
                'OCIGissCX'                 : VLIDORT_POLAR_.vector_ocigisscx,
                'OCIGissCX_NOBM_CLOUD'      : VLIDORT_POLAR_.vector_ocigisscx_nobm_cloud,
                'ROT_CALC'                  : VLIDORT_POLAR_.rot_calc}



def get_geom(Iscan,Icross):
    inFile = '/nobackup/PACE/L1B/Y2020/M03/D24/OCI2020084005000.L1B_PACE.nc'
    nc = Dataset(inFile)
    grp = nc.groups['geolocation_data']
    sza = grp.variables['solar_zenith'][:]
    saa = grp.variables['solar_azimuth'][:]
    vza = grp.variables['sensor_zenith'][:]
    vaa = grp.variables['sensor_azimuth'][:]

    # make azimuths clockwise from north
    I = saa < 0
    saa[I] = 360. + saa[I]

    I = vaa < 0
    vaa[I] = 360. + vaa[I]

    # define SAA according to photon travel direction
    saa = saa + 180.0
    I = saa >= 360.
    saa[I] = saa[I] - 360.

    raa = vaa - saa   

    I = raa < 0
    raa[I] = raa[I] + 360.0

    saa = grp.variables['solar_azimuth'][:]
    vaa = grp.variables['sensor_azimuth'][:]


    sza = np.array([sza[Iscan,Icross]])
    vza = np.array([vza[Iscan,Icross]])
    raa = np.array([raa[Iscan,Icross]])
    saa = np.array([saa[Iscan,Icross]])
    vaa = np.array([vaa[Iscan,Icross]])

    return sza,vza,raa,saa,vaa

def get_ROT(ch,pe,te,ze,verbose=True):
    """
    calculate ROT and depol ratio
    """

    args = [ch, pe, ze, te, MISSING, verbose]
    ROT, depol_ratio, rc = VLIDORT_POLAR_.rot_calc(*args)


    return ROT,depol_ratio

def get_TOA_unpack(args):
    return get_TOA(*args)

def get_TOA(channel,
            F0,
            ROT,depol_ratio,
            tau,ssa,pmom,
            alpha,
            SZA,VZA,RAA,            
            km,pe,te,ze,
            nstreams,            
            albedoType,
            sleave,
            U10m,V10m,mr,
            verbose):
    """
    Do RT calculation to get TOA
    """

    # wrapper function based on albedo
    vlidortWrapper = WrapperFuncs[albedoType]
    
    if albedoType == 'OCIGissCX':
        args = [channel, nstreams, plane_parallel, ROT, depol_ratio, alpha, tau, ssa, pmom,
                pe, ze, te,
                U10m, V10m, mr,
                SZA, RAA, VZA,
                F0,
                MISSING,
                verbose]
        I, reflectance, surf_reflectance, Q, U, BR_Q, BR_U, rc = vlidortWrapper(*args)
    else:
        args = [channel, nstreams, plane_parallel, ROT, depol_ratio, alpha, tau, ssa, pmom,
                tau, ssa, pmom,
                tau, ssa, pmom,
                pe, ze, te,
                U10m, V10m, mr,
                sleave, False,
                SZA, RAA, VZA,
                F0,
                MISSING,
                verbose]
        I, reflectance, surf_reflectance, Q, U, BR_Q, BR_U, rc, adjusted_sleave = vlidortWrapper(*args)


    return I,reflectance,surf_reflectance


def writenc(outFile,channels,
            F0,
            SZA,SAA,VZA,VAA,
            I,BR,reflectance,
            ROD,depol_ratio,
            alphaD,
            pe,te,ze,
            U10m,V10m,mr):



    nch = len(channels)
    km  = len(pe)

    # Write data to netcdf file
    nc = Dataset(outFile,'w',format='NETCDF4_CLASSIC')
    nc.title = 'Line-by-Line TOA radiance calculation for one OCI pixel. Model 5 - US Standard 1962'

    dc = nc.createDimension('channels',nch)
    dn = nc.createDimension('npixel',1)
    dk = nc.createDimension('leve',km)

    ch = nc.createVariable('channels','f4',('channels',))
    ch.long_name = "Wavelength in nm"
    ch[:] = channels

    f  = nc.createVariable('solar_irradiance','f4',('channels',))
    f.long_name = 'Thuillier solar irradiance spectrum'
    f.units = 'uW/cm^2/nm'
    f[:] = F0

    rad = nc.createVariable('I','f4',('channels',))
    rad.long_name = "sun normalized TOA radiance"
    rad[:] = I

    ref = nc.createVariable('R','f4',('channels',))
    ref.long_name = "TOA reflectance"
    ref[:] = reflectance

    sr = nc.createVariable('surface_reflectance','f4',('channels',))
    sr.long_name = "surface bidirectional reflectance"
    sr[:] = BR  

    a = nc.createVariable('sza','f4',('npixel',))
    a.long_name = "solar zenith angle"
    a.units = 'degrees'
    a[:] = SZA    

    a = nc.createVariable('saa','f4',('npixel',))
    a.long_name = "solar azimiuth angle"
    a.units = 'degrees 0-360 clockwise from north'
    a[:] = SAA

    a = nc.createVariable('vza','f4',('npixel',))
    a.long_name = "sensor zenith angle"
    a.units = 'degrees'
    a[:] = VZA

    a = nc.createVariable('vaa','f4',('npixel',))
    a.long_name = "sensor azimiuth angle"
    a.units = 'degrees 0-360 clockwise from north'
    a[:] = SAA 

    rod = nc.createVariable('ROD','f4',('channels',))
    rod.long_name = "Rayleigh Optical Depth"
    rod[:] = ROD

    depol = nc.createVariable('depol_ratio','f4',('channels',))
    depol.long_name = 'Rayleigh depolarization ratio'
    depol[:] = depol_ratio

    alpha = nc.createVariable('ALPHA','f4',('channels',))
    alpha.long_name = "Trace Gas Absorption Optical Depth"
    alpha[:] = alphaD

    ke = nc.createVariable('PE','f4',('leve',))
    ke.long_name = 'Pressure at layer edge'
    ke.units     = 'Pa'
    ke[:] = pe

    ke = nc.createVariable('TE','f4',('leve',))
    ke.long_name = 'Temperature at layer edge'
    ke.units     = 'K'
    ke[:] = te

    ke = nc.createVariable('ZE','f4',('leve',))
    ke.long_name = 'height above surface'
    ke.units = 'm'
    ke[:] = ze

    wi = nc.createVariable('U10M','f4',('npixel',))
    wi.long_name = 'U10M wind speed'
    wi.units = 'm/s'
    wi[:] = U10m

    wi = nc.createVariable('V10M','f4',('npixel',))
    wi.long_name = 'V10M wind speed'
    wi.units = 'm/s'
    wi[:] = V10m

    wi = nc.createVariable('mr','f4',('npixel',))
    wi.long_name = 'ocean water refractive index'
    wi[:] = mr

    nc.close()    

def get_PTWV_profile(inFile,model=5):
    """
    Read in height [km],pressure [mb], temperature [K], water vapor vmr profile [ppm]
    for selected model:
    0:Tropical
    1:Mid Latitude Summer
    2:Mid Latitude Winter
    3:Subarctic Summer
    4:Subarctic Winter
    5:US Standard 1962
    6:User Defined Model
    """
    nc = Dataset(inFile)
    # make array because interp doesn't take masked arrays
    pe = np.array(nc.variables['p'][model,:])  
    te = np.array(nc.variables['t'][model,:])
    ze = np.array(nc.variables['h'][model,:])    
    vmre =np.array( nc.variables['vmr'][model,:])

    km   = len(ze) - 1

    # get DP, from Amir's code
    DP= np.empty([km])
    for i in range(km):
        DP[i] = (pe[i] - pe[i+1])/1013

    # convert mb to Pascal
    pe = pe *100
    DELP = pe[:-1] - pe[1:]

    # get T middles
    DELP = pe[:-1] - pe[1:]
    pm = pe[:-1] - 0.5*DELP
    f = interp1d(pe,te)
    tm = f(pm)


    # convert km to m
    ze = ze*1000.
    dz = ze[1:] - ze[:-1]

    # calculate air number density [molecules/m3]
    g = 9.80616 # gravity m/s2
    Na = 6.022e23 # avogadro's number
    MW_AIR =  28.964*1e-3 # kg/mole

    AIRDENS = DELP/(dz*g) # kg/m3
    rho = AIRDENS/MW_AIR      # moles/m3
    rho = rho*Na  #[molecules/m3]

    return km, pe, te, tm, ze, dz, rho, vmre, DP, AIRDENS 


def get_abs(inFile):
    """
    Read aborption coefficients from Amir's HITRAN calculations
    """
    nc = Dataset(inFile)
    # wavenumber [cm-1]
    waveno = np.array(nc.variables['waveno'][:])
    # abosrption coefficients [not sure about units]
    abs_o2  = nc.variables['abscf_o2'][:]
    abs_h2o = nc.variables['abscf_h2o'][:]
    abs_co  = nc.variables['abscf_co'][:]
    abs_co2 = nc.variables['abscf_co2'][:]
    abs_ch4 = nc.variables['abscf_ch4'][:]
    abs_n2o = nc.variables['abscf_n2o'][:]

    nc.close()

    return abs_o2, abs_h2o, abs_co, abs_co2, abs_ch4, abs_n2o, waveno


def get_rsr(inFile):
    """
    Read in OCI RSR File
    """
    hdf = SD(inFile, SDC.READ)
    rsr = hdf.select('RSR')[:]
    wav_rsr = hdf.select('rsrwave')[:]
    wav_oci = hdf.select('wave')[:]
    hdf.end()

    return rsr, wav_rsr, wav_oci

def get_alpha(A,VMR,rho,dz):
    """
    Calculate Absorption optical depth profile
    A - absorption coefficient [m2/molecule]
    VMR - trace gas mixing ratio [vol/vol, dimensionless]
    rho - air number density [molecules/m3]
    ze  - profile layer thickness [m]
    """

    # convert vmr to molecules/m3
    nxe = VMR*rho

    # get the optical depth subcolumns
    km, nch = A.shape
    alpha = np.zeros([km,nch])
    for i in range(km):
        alpha[i,:] = nxe[i]*dz[i]*A[i,:]

    return alpha


# ---
def read_o3(inFile,te):
    f = open(inFile)
    nhead = 10
    for i in range(nhead):
        hh = f.readline()

    wav_o3, c0, c1, c2 = [],[],[],[]
    for l in f:
        a = np.array(l.split()).astype(float)
        wav_o3.append(a[0])
        c0.append(a[1])
        c1.append(a[2])
        c2.append(a[3])
    f.close()

    wav_o3 = np.array(wav_o3)
    c0     = np.array(c0)
    c1     = np.array(c1)
    c2     = np.array(c2)
    T0     = 273.15

    # calculate xsec for te
    nwav = len(wav_o3)
    km   = len(te)
    xsec_o3 = np.zeros([km,nwav])

    for i,t in enumerate(te):
        xsec_o3[i,:] = c0 + c1*(t-T0) + c2*(t-T0)**2

    xsec_o3 = xsec_o3*1e-20

    # convert from cm2/molecule to m2/molecule
    xsec_o3 = xsec_o3*1e-4

    return wav_o3, xsec_o3

# --
def read_ROD_table(inFile):
    f = open(inFile)

    for i in range(16):
        f.readline() #header

    wav = []
    rod = []
    depol = []
    for l in f:
        w, r, d = l.split()
        wav.append(w)
        rod.append(r)
        depol.append(d)

    f.close()

    wav = np.array(wav).astype('float')
    rod = np.array(rod).astype('float')
    depol = np.array(depol).astype('float')

    return wav, rod, depol

#------------------------------------ M A I N ------------------------------------

if __name__ == "__main__":

    outRoot = 'hyperTest/'
    outFile   = '{}/outputs/hyperTest_CK_Thuillier_g5nr.nc4'.format(outRoot)


    # Pressure [Pa], temperature [K], height [m], water vapor [ppm]  profile - standard atmosphere
    # used to make OCI look up tables
    # rhoe = air number density [molecules/m3]
    inFile    = '{}/atrem_tpvmr.nc'.format(outRoot)
    km, pe, te, tm, ze, dz, rho, h2oe, DP, AIRDENS =  get_PTWV_profile(inFile)

    # Read in G5NR CO, CO2, O3, RH
    Iscan  = 600
    Icross = 1000
    LevelB = '/nobackup/PACE/LevelB/Y2006/M03/D24'
    inFile = '{}/pace-g5nr-std.lb.chm_Nv.20060324_005000.nc4'.format(LevelB)
    nc = Dataset(inFile)
    co = nc.variables['CO'][0,:,Iscan,Icross]    # VMR
    co2 = nc.variables['CO2'][0,:,Iscan,Icross]  # VMR
    o3 = nc.variables['O3'][0,:,Iscan,Icross]  #kg/kg
    nc.close()
    inFile = '{}/pace-g5nr-std.lb.aer_Nv.20060324_005000.nc4'.format(LevelB)
    nc = Dataset(inFile)
    h2o = nc.variables['WV_VMR'][0,:,Iscan,Icross]  #ppm
    nc.close()
    
    # flip from bottom to top
    co  = co[::-1]
    co2 = co2[::-1]
    o3  = o3[::-1]
    h2o = h2o[::-1]

    # Read in alpha_table
    inFile = 'alphaTable_v0/alpha_CK_Thuillier_o3.nc4'
    nc = Dataset(inFile)
    all_wl = np.array(nc.variables['channels'][:])
    g_bins   = nc.variables['g_bins'][:]
    alpha_o2 = nc.variables['alpha_o2'][:]
    ROD      = nc.variables['ROD'][:]
    depol    = nc.variables['depol_ratio'][:]
    F0_int   = nc.variables['solar_irradiance'][:]

    nc.close()

    alpha = alpha_o2

#    # integrate air density in each layer
#    rhoint = rho*dz


    # get absorption optical depth with new aborption coefficient
#    co_vmr = 0.1*1e-6
#    co_vmr = co
#    alpha_co = get_alpha(abs_co_z,co_vmr,rho,dz)
    
#    o2_vmr = 0.21
#    alpha_o2 = get_alpha(abs_o2_z,o2_vmr,rho,dz)
    
#    co2_vmr = 400.*1.0E-06
#    co2_vmr = co2
#    alpha_co2 = get_alpha(abs_co2_z,co2_vmr,rho,dz)

#    ch4_vmr = 1.8*1.0E-06
#    alpha_ch4 = get_alpha(abs_ch4_z,ch4_vmr,rho,dz)

#    n2o_vmr = 0.3*1.0E-06
#    alpha_n2o = get_alpha(abs_n2o_z,n2o_vmr,rho,dz)

#    h2o_vmr = h2o*1.0e-6
#    alpha_h2o = get_alpha(abs_h2o_z,h2o_vmr,rho,dz)

    # add up all the alphas
#    alpha = alpha_h2o + alpha_n2o + alpha_ch4 + alpha_co2 + alpha_o2 + alpha_co

#    # ----
#    # OZONE Stuff
#    # ---
#    # read xsec
#    inFile = 'hyperTest/o3_bremen/Ozone_abs_x_wTemperatureFit.dat'
#    # wav [nm], C0, C1(T), C2(T^2)
#    # xsec is in m2/molecule
#    wav_o3,abs_o3 = read_o3(inFile,tm)
#    # reverse so going from max to min wavelength
#    wav_o3 = wav_o3[::-1]
#    abs_o3 = abs_o3[:,::-1]

#    # interpolate to LBL wavelengths
#    # append zeros to max lbl wavelength
#    wav_new = np.arange(wav_o3.max()+0.1,wav_abs.max()+0.1,0.1)
#    # reverse
#    wav_new = wav_new[::-1]
#    wav_o3  = np.append(wav_new,wav_o3)
#    nnew = len(wav_new)
#    abs_o3 = np.append(np.zeros([km,nnew]),abs_o3,axis=1)

#    abs_o3_lbl = np.zeros(abs_h2o_z.shape)
#    for k in range(km):
#       # xsec_f = interp1d(wav_o3,abs_o3[k,:],kind='linear')
#        abs_o3_lbl[k,:] = xsec_f(wav_abs)
    
#    # append UV-Vis to LBL that stops at 555
#    i = wav_o3 < wav_abs.min()
#    all_wl = np.append(wav_abs,wav_o3[i])
#    abs_o3_lbl = np.append(abs_o3_lbl,abs_o3[:,i],axis=1)
    
#    # convert mass mixing ratio to molecules/m3
#    Na = 6.022e23 # avogadro's number
#    O3_MW  = 48.0*1e-3     # kg/mole

#    o3_conc = o3*AIRDENS   # kg/m3
#    o3_conc = o3_conc*Na/O3_MW   # molecules/m3

#    # get the optical depth subcolumns
#    alpha_o3 = np.zeros(abs_o3_lbl.shape)
#    for i in range(km):
#        alpha_o3[i,:] = o3_conc[i]*dz[i]*abs_o3_lbl[i,:]

    nwav = len(all_wl)
    
#    # add ozone to total alpha
#    # extend array down to uv
#    nnew = nwav - len(wav_abs)
#    alpha = np.append(alpha,np.zeros([km,nnew]),axis=1)
#    alpha = alpha + alpha_o3
    
#    # limit to wavelengths covered by RSR
#    I = all_wl <= wav_rsr.max()+1.0
#    all_wl = all_wl[I]
#    alpha   = alpha[:,I]
#    alpha[:] = 0.0


#    # get ROD from oci_tables
#    inFile = 'oci_tables/rayleigh_bodhaine.txt'
#    wav,rod,depol = read_ROD_table(inFile)

#    # interpolate to lbl wavelengths
#    rod_f = interp1d(wav,rod,kind='linear')
#    rod_lbl = rod_f(all_wl)
#    depol_f = interp1d(wav,depol,kind='linear')
#    depol_lbl = depol_f(all_wl)
    
    # flip everything vertically so going from top of atmosphere to surface
    pe = pe[-1::-1]
    te = te[-1::-1]
    ze = ze[-1::-1]
#    alpha = alpha[-1::-1,:]

    # add dimension to be in km+1,nobs
    pe.shape = (km+1,1)
    te.shape = (km+1,1)
    ze.shape = (km+1,1)


    # read in granule geometry
    SZA,VZA,RAA,SAA,VAA = get_geom(Iscan,Icross)
    csza = np.cos(np.radians(SZA))

#    # Read in solar irradiance spectrum
#    # second dim = wavelength, irradiance
#    # units=nm, uW/cm^2/nm
#    inFile = '{}/Thuillier_F0.npy'.format(outRoot)
#    F0 = np.load(inFile)
#    # interpolate to wavelengths 
#    F0_f = interp1d(F0[:,0],F0[:,1],kind='linear',fill_value="extrapolate")
#    F0_int = F0_f(all_wl)

    # --------------
    # surface stuff
    # --------------

    # SLEAVE
    LevelB = '/nobackup/PACE/LevelB/surface/SLEAVE/NOBM/Y2006/M03/D24'
    inFile = '{}/pace-g5nr.lb.sleave.20060324_005000.nc4'.format(LevelB)
    nc     = Dataset(inFile)
    nobm_wav = nc.variables['wavelength'][:]
    rrs      = nc.variables['rrs'][0,:,Iscan,Icross]
    nc.close()
    rrs_f    = interp1d(nobm_wav,rrs,kind='linear',fill_value=0.0,bounds_error=False)
    sleave   = rrs_f(all_wl)*csza
    
    # wind speed
    LevelB = '/nobackup/PACE/LevelB/Y2006/M03/D24'
    inFile = '{}/pace-g5nr-std.lb.met_Nv.20060324_005000.nc4'.format(LevelB)
    nc     = Dataset(inFile)
    U10m   = np.array([nc.variables['U10M'][0,Iscan,Icross]])
    V10m   = np.array([nc.variables['V10M'][0,Iscan,Icross]])

    # water refractive index
    mr   = 1.334

    # RT stuff
    nstreams = 12
    albedoType = 'OCIGissCX_NOBM_CLOUD'
    albedoType = 'OCIGissCX'

    # loop through channels
    nproc = 50
    nwl   = len(all_wl)
#    nwl   = 100
#    wlstep = int(nwl/nproc)
    wlstep  = 10
    sys.exit()
    args = []
    ROD  = []
    depol = []
    sys.exit()
    for ich in np.arange(0,nwl,wlstep):
        endch = ich + wlstep
        if endch > nwl:
            endch = nwl
        nch = endch - ich
        ch = all_wl[ich:endch]

        # Get Rayleigh
        # ROT shape is nlev,1,nch
        ROT, depol_ratio = get_ROT(ch,pe,te,ze,verbose=False)    
       
        ROT = ROT*rod_lbl[ich:endch]/np.squeeze(ROT.sum(axis=0))
        ROD.append(np.squeeze(ROT.sum(axis=0)))

        depol_ratio = depol_lbl[ich:endch]
        depol.append(depol_ratio)

        # trace gas
        alpha_ch = alpha[:,ich:endch]
        alpha_ch.shape = (km,1,nch)

        # AOP vectors [km,nch,nobs]
        tau = np.zeros([km,nch,1])
        ssa = np.zeros([km,nch,1])
        pmom = np.zeros([km,nch,1,30,6])
 
        # water refractive index
        mr_in = np.ones(nch)
        mr_in = mr_in*mr

        # solar irradiance
        F0_in = F0_int[ich:endch]
        F0_in.shape = (nch,1)

        # sleave
        sleave_ch = sleave[ich:endch]

        args.append([ch,
                     F0_in,
                     ROT,depol_ratio,
                     tau,ssa,pmom,
                     alpha_ch,
                     SZA,VZA,RAA,
                     km,pe,te,ze,
                     nstreams,
                     albedoType,
                     sleave_ch,
                     U10m,V10m,mr_in,
                     False])

#       I,reflectance,BR =  get_TOA(ch,
#                                    ROT,depol_ratio,
#                                    tau,ssa,pmom,
#                                    alpha_ch,
#                                    SZA,VZA,RAA,
#                                    km,pe,te,ze,
#                                    nstreams,
#                                    albedoType,
#                                    U10m=U10m,V10m=V10m,mr=np.array([1.334]),
#                                    verbose=False)


    
    # use multiprocessing
    p = Pool(nproc)
    result = p.map(get_TOA_unpack,args)
    I = []
    reflectance = []
    BR = []
    for r in result:    
        I_r,reflectance_r,BR_r = r
        I.append(np.squeeze(I_r))
        reflectance.append(np.squeeze(reflectance_r))
        BR.append(np.squeeze(BR_r))

    p.close()
    p.join()
    
    # concatenate arrays
    ROD = np.concatenate(ROD)
    depol_ratio = np.concatenate(depol)

    I  = np.concatenate(I)
    reflectance = np.concatenate(reflectance)
    BR = np.concatenate(BR)
                
    alphaD = alpha.sum(axis=0)            
    # write to outFile
    writenc(outFile,all_wl,
            F0_int,
            SZA,SAA,VZA,VAA,
            I,BR,reflectance,
            ROD,depol_ratio,
            alphaD,
            pe,te,ze,
            U10m,V10m,mr)    

