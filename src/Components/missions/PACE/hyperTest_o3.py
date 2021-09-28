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


def read_ROD_table():
    f = open('amir/OCI_ROD_Table_adjusted.txt','r')
    f.readline() #header
    wav = f.readline().split()
    f.readline() #wav center
    f.readline() #wav width
    f.readline() #F0
    rod = f.readline().split()
    depol = f.readline().split() 

    f.close()

    rod = rod[2:]
    rod = np.array(rod).astype('float')

    depol = depol[2:]
    depol = np.array(depol).astype('float')

    wav = wav[3:]
    wav = np.array(wav).astype('float')

    return wav,rod,depol


def get_channels():
    inFile = '/nobackup/PACE/L1B/Y2020/M03/D24/OCI2020084005000.L1B_PACE.nc'
    nc = Dataset(inFile)
    grp = nc.groups['sensor_band_parameters']
    blue = grp.variables['blue_wavelength'][:]
    red  = grp.variables['red_wavelength'][:]
    swir = grp.variables['SWIR_wavelength'][:]

    return blue,red,swir

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
            U10m,V10m,mr,
            verbose):
    """
    Do RT calculation to get TOA
    """

    # wrapper function based on albedo
    vlidortWrapper = WrapperFuncs[albedoType]

    if 'CX' in albedoType:
        args = [channel, nstreams, plane_parallel, ROT, depol_ratio, alpha, tau, ssa, pmom,
                pe, ze, te,
                U10m, V10m, mr,
                SZA, RAA, VZA,
                F0,
                MISSING,
                verbose]

    I, reflectance, surf_reflectance, Q, U, BR_Q, BR_U, rc = vlidortWrapper(*args)    


    return I,reflectance,surf_reflectance


def writenc(outFile,channels,
            F0,
            SZA,SAA,VZA,VAA,
            I,BR,reflectance,
            ROD,depol_ratio,
            alphaD,
            pe,te,ze,
            U10m,V10m,mr):

    nch = len(I)
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

    # convert km to m
    ze = ze*1000.

    # calculate air number density [molecules/m3]
    R = 8.3145  # gas constant [J mol-1 K-1]
    Na = 6.022e23 # avogadro's number
    airdens = pe/(R*te) # moles/m3
    airdens = airdens*Na

    return km, pe, te, ze, airdens, vmre, DP 


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

def get_abs_no2(inFile):
    """
    Read in NO2 cross section
    """
    f = open(inFile)
    # header
    for i in range(8):
        ll = f.readline()

    scale_factor = float(f.readline())
    xsec = []
    wav  = []
    for l in f:
        l = l.split()
        wav.append(l[0])
        xsec.append(l[1])

    # nm
    wav = np.array(wav).astype(float)
    xsec = np.array(xsec).astype(float)

    # mutiply by scale factor
    # units are  [cm2/molecule]
    xsec = xsec*scale_factor

    # convert to m2/molecule
    xsec = xsec*1e-4

    return xsec, wav


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

def get_alpha(A,VMR,rhoe,ze):
    """
    Calculate Absorption optical depth profile
    A - absorption coefficient [m2/molecule]
    VMR - trace gas mixing ratio [vol/vol, dimensionless]
    rhoe - air number density [molecules/m3]
    ze  - profile altitude [m]
    """

    # convert vmr to molecules/m3
    nxe = VMR*rhoe

    # integrate to get the optical depth subcolumns
    km, nch = A.shape
    alpha = np.zeros([km,nch])
    for i in range(km):
        c1 = A[i,:]*nxe[i]
        c2 = A[i,:]*nxe[i+1]
        c1.shape = (1, nch)
        c2.shape = (1, nch)
        c  = np.append(c1,c2,axis=0)
        alpha[i,:] = np.trapz(c,ze[i:i+2],axis=0)

    #alpha = np.trapz(A*nxe,ze)
    
    return alpha

def get_o3_profile(inFile):
    """
    read o3 profile (molecules/cm3) from mcclams.dat
    p = hPa
    t = K
    airdens = molecules/cm3
    """
    f = open(inFile)
    # header
    f.readline()
    f.readline()
    
    z = []
    p = []
    t = []
    airdens = []
    o3_conc = []
    for l in f:
        aa = np.array(l.split()).astype(float)
        z.append(aa[0])
        p.append(aa[1])
        t.append(aa[2])
        airdens.append(aa[3])
        o3_conc.append(aa[4])
    
    f.close()
    z = np.array(z)
    p = np.array(p)
    t = np.array(t)
    airdens = np.array(airdens)
    o3_conc = np.array(o3_conc)

    # o3 mixing ratio, unitless
    o3_vmr = o3_conc/airdens

    # get altitude in m
    z = z*1e3

    # get conc in molecules/m3
    o3_conc = o3_conc*1e6

    return z,o3_conc

def get_abs_o3(inFile,wav_lbl):
    """
    Read Amir's Ozone Cross Sections
    Interpolate to LBL wavelengths
    Keep resolution down to UV
    """
    f = open(inFile)

    # header 
    nhead = 37
    for i in range(nhead):
        f.readline()

    wav = []
    xsec = []
    for l in f:
        w,x = l.split()
        wav.append(float(w))
        xsec.append(float(x))

    f.close()
    wav = np.array(wav)
    xsec = np.array(xsec)

    # reverse so going from max to min wavelength
    wav = wav[::-1]
    xsec = xsec[::-1]

    # cross sections are for 1 DU = 2.6967e19 molecules/cm3 @ STP
    # divide by 1DU to get to units cm2/molecule
    xsec = xsec/2.6967e19

    # multipy by 1e-4 to get m2/molecule
    xsec = xsec*1e-4

    # append zeros to max lbl wavelength
    wav_new = np.arange(wav.max()+0.1,wav_lbl.max()+0.1,0.1)
    # reverse
    wav_new = wav_new[::-1]
    wav_o3  = np.append(wav_new,wav)
    nnew = len(wav_new)
    xsec_o3 = np.append(np.zeros(nnew),xsec)

    # interpolate ozone xsec to lbl wavelengths
    xsec_f = interp1d(wav_o3,xsec_o3,kind='linear')
    xsec_lbl = xsec_f(wav_lbl)

    # append UV-Vis to lbl
    i = wav_o3 < wav_lbl.min()
    wav_total = np.append(wav_lbl,wav_o3[i])
    xsec_total = np.append(xsec_lbl,xsec_o3[i])

    return wav_total,xsec_total

#------------------------------------ M A I N ------------------------------------

if __name__ == "__main__":

    outRoot = 'hyperTest/'
    outFile   = '{}/outputs/hyperTest_LBL_Thuillier_o3.nc4'.format(outRoot)

    # Read in OCI RSR
    inFile = '{}/OCI_RSR_v0.hdf'.format(outRoot)
    rsr, wav_rsr, wav_oci = get_rsr(inFile)

    # Pressure [Pa], temperature [K], height [m], water vapor [ppm]  profile - standard atmosphere
    # used to make OCI look up tables
    # rhoe = air number density [molecules/m3]
    inFile    = '{}/atrem_tpvmr.nc'.format(outRoot)
    km, pe, te, ze, rhoe, h2oe, DP =  get_PTWV_profile(inFile)
    zm = (ze[1:] - ze[0:-1])*0.5 + ze[0:-1]

    # Read in Amir's absorption cross sections and wavelengths
    # covers 555 nm  - 3.3 um
    inFile = '{}/abscf_gas.nc'.format(outRoot)
    abs_o2, abs_h2o, abs_co, abs_co2, abs_ch4, abs_n2o, wavno_abs = get_abs(inFile)

    # from Amir's code:
    abs_o2[:,9001:106600] = abs_o2[:,9001:106600]*2.6

    # convert wavenumber [1/cm] to wavelength [nm]
    wav_abs = 1e7*(1./wavno_abs)

    # convert gas absorption coefficients to m2/molecule
    # ---------------------------------------------------
    Q = 2.15199993E+25
    C = (Q*28.966) / 6.0225e23 / 1e-6

    # integrate air density in each layer
    rhoint = np.zeros(km)
    for i in range(km):
        rhoint[i] = np.trapz(rhoe[i:i+2],ze[i:i+2])

    DP.shape = (km,1)
    rhoint.shape = (km,1) 
    abs_o2_z  = abs_o2*C*DP/rhoint
    abs_co_z  = abs_co*C*DP/rhoint
    abs_co2_z = abs_co2*C*DP/rhoint
    abs_ch4_z = abs_ch4*C*DP/rhoint
    abs_n2o_z = abs_n2o*C*DP/rhoint
    abs_h2o_z = abs_h2o*C*DP/rhoint

    # get absorption optical depth with new aborption coefficient
    co_vmr = 0.1*1e-6
    alpha_co = get_alpha(abs_co_z,co_vmr,rhoe,ze)
    
    o2_vmr = 0.21
    alpha_o2 = get_alpha(abs_o2_z,o2_vmr,rhoe,ze)

    co2_vmr = 400.*1.0E-06
    alpha_co2 = get_alpha(abs_co2_z,co2_vmr,rhoe,ze)

    ch4_vmr = 1.8*1.0E-06
    alpha_ch4 = get_alpha(abs_ch4_z,ch4_vmr,rhoe,ze)

    n2o_vmr = 0.3*1.0E-06
    alpha_n2o = get_alpha(abs_n2o_z,n2o_vmr,rhoe,ze)

    
    # scale water vapor so total precipitable water 
    # is equal to 1 cm
    # ----------------

    # first calculate water vapor column [molecules/m2]
    h2ocol = np.trapz(h2oe*1e-6*rhoe,ze)
    # normalize profile so water vapor column expressed as total precipitable water is equal to 1 cm 
    # use 1 mm of rainfall = 1 kg/m2
    # or 1 cm = 10 kg/m2
    # 10 kg/m2 is equal to 3.34e22 water molecules/cm2
    h2ocolnew = 3.34e22
    # convert to meters
    h2ocolnew = h2ocolnew*1e4
  
    h2oenew = h2oe*(h2ocolnew/h2ocol)
    # get in vmr units
    h2oe_vmr = 1e-6*h2oenew
    alpha_h2o = get_alpha(abs_h2o_z,h2oe_vmr,rhoe,ze)

    # add up all the alphas
    alpha = alpha_h2o + alpha_n2o + alpha_ch4 + alpha_co2 + alpha_o2 + alpha_co

    # ----
    # OZONE Stuff
    # ---
    # Read in mcclams ozone concentration profile [molecules/m3]
    inFile = '{}/mcclams.dat'.format(outRoot)
    z_o3, o3_conc = get_o3_profile(inFile)

    # interpolate o3_conc to PT profile
    o3_f = interp1d(z_o3,o3_conc,kind='linear',fill_value="extrapolate")
    o3_conce = o3_f(ze)

    # Read in OCI ozone cross section in m2/molecule 
    # interpolate to LBL wavelengths
    inFile = '{}/k_o3.txt'.format(outRoot)
    wav_o3, abs_o3 = get_abs_o3(inFile,wav_abs)
    abs_o3.shape = (1,len(wav_o3))
    for k in range(km-1):
        abs_o3 = np.append(abs_o3,abs_o3[0:1,:],axis=0)

    # integrate conc_o3*xsec_o3 to get o3  alpha
    nwav = len(wav_o3)
    alpha_o3 = np.zeros([km,nwav])
    for i in range(km):
        o3 = o3_conce[i:i+2]
        o3.shape = (2,1)
        alpha_o3[i,:] = np.trapz(o3*abs_o3[i:i+2,:],ze[i:i+2],axis=0)

    # add ozone to total alpha
    # extend array down to uv
    nnew = nwav - len(wav_abs)
    alpha = np.append(alpha,np.zeros([km,nnew]),axis=1)
    alpha = alpha + alpha_o3
    all_wl = wav_o3
    
    # limit to wavelengths covered by RSR
    I = all_wl <= wav_rsr.max()+1.0
    all_wl = all_wl[I]
    alpha   = alpha[:,I]
    
    # flip everything vertically so going from top of atmosphere to surface
    pe = pe[-1::-1]
    te = te[-1::-1]
    ze = ze[-1::-1]
    alpha = alpha[-1::-1,:]

    # add dimension to be in km+1,nobs
    pe.shape = (km+1,1)
    te.shape = (km+1,1)
    ze.shape = (km+1,1)


    # read in granule geometry
    Iscan  = 600
    Icross = 1000
    SZA,VZA,RAA,SAA,VAA = get_geom(Iscan,Icross)

    # Read in solar irradiance spectrum
    # second dim = wavelength, irradiance
    # units=nm, uW/cm^2/nm
    inFile = '{}/Thuillier_F0.npy'.format(outRoot)
    F0 = np.load(inFile)
    # interpolate to c-k wavelength bins
    F0_f = interp1d(F0[:,0],F0[:,1],kind='linear',fill_value="extrapolate")
    F0_int = F0_f(all_wl)

    # set up some vlidort stuff
    albedoType = 'OCIGissCX'
    U10m = np.array([3.0])
    V10m = np.array([4.0])
    mr   = 1.334
    nstreams = 12

    # loop through channels
    nproc = 50
    nwl   = len(all_wl)
#    nwl   = 100
#    wlstep = int(nwl/nproc)
    wlstep  = 10

    args = []
    ROD  = []
    depol = []
    for ich in np.arange(0,nwl,wlstep):
        endch = ich + wlstep
        if endch > nwl:
            endch = nwl
        nch = endch - ich
        ch = all_wl[ich:endch]

        # Get Rayleigh
        ROT, depol_ratio = get_ROT(ch,pe,te,ze,verbose=False)           
        ROD.append(np.squeeze(ROT.sum(axis=0)))
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


        args.append([ch,
                     F0_in,
                     ROT,depol_ratio,
                     tau,ssa,pmom,
                     alpha_ch,
                     SZA,VZA,RAA,
                     km,pe,te,ze,
                     nstreams,
                     albedoType,
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


