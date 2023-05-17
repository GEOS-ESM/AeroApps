#!/usr/bin/env python3

"""
    compared o3 cross sections

    Patricia Castellanos, April 2020

"""

import os
import sys
from   netCDF4 import Dataset
from   netCDF4 import Dataset as ncread
import numpy   as np
from MAPL.constants import *
from scipy.interpolate import interp1d
import scipy.integrate as integrate
from pyhdf.SD import SD, SDC

format   = 'NETCDF4_CLASSIC'

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


#------------------------------------ M A I N ------------------------------------

if __name__ == "__main__":


    # Read in OCI RSR
    inFile = '../OCI_RSR_v0.hdf'
    rsr, wav_rsr, wav_oci = get_rsr(inFile)

    # Pressure [Pa], temperature [K], height [m], water vapor [ppm]  profile - standard atmosphere
    # used to make OCI look up tables
    # rhoe = air number density [molecules/m3]
    inFile    = '../atrem_tpvmr.nc'
    km, pe, te, tm, ze, dz, rho, h2oe, DP, AIRDENS =  get_PTWV_profile(inFile)

    # Read in G5NR CO, CO2, O3, RH
    Iscan  = 600
    Icross = 1000
    LevelB = '/nobackup/PACE/LevelB/Y2006/M03/D24'
    inFile = '{}/pace-g5nr-std.lb.chm_Nv.20060324_005000.nc4'.format(LevelB)
    nc = Dataset(inFile)
    o3 = nc.variables['O3'][0,:,Iscan,Icross]  #kg/kg
    nc.close()
    
    # flip from bottom to top
    o3  = o3[::-1]

    # Read in Amir's absorption cross sections and wavelengths
    # covers 555 nm  - 3.3 um
    inFile = '../abscf_gas.nc'
    abs_o2, abs_h2o, abs_co, abs_co2, abs_ch4, abs_n2o, wavno_abs = get_abs(inFile)

    # convert wavenumber [1/cm] to wavelength [nm]
    wav_abs = 1e7*(1./wavno_abs)

    # integrate air density in each layer
    rhoint = rho*dz

    DP.shape = (km,1)
    rhoint.shape = (km,1) 
    
    # ----
    # OZONE Stuff
    # ---
    # read xsec
    inFile = '../o3_bremen/Ozone_abs_x_wTemperatureFit.dat'
    # wav [nm], C0, C1(T), C2(T^2)
    # xsec is in m2/molecule
    wav_o3,abs_o3 = read_o3(inFile,tm)
    # reverse so going from max to min wavelength
    wav_o3 = wav_o3[::-1]
    abs_o3 = abs_o3[:,::-1]

    # interpolate to LBL wavelengths
    # append zeros to max lbl wavelength
    wav_new = np.arange(wav_o3.max()+0.1,wav_abs.max()+0.1,0.1)
    # reverse
    wav_new = wav_new[::-1]
    wav_o3  = np.append(wav_new,wav_o3)
    nnew = len(wav_new)
    abs_o3 = np.append(np.zeros([km,nnew]),abs_o3,axis=1)

    abs_o3_lbl = np.zeros([km,len(wav_abs)])
    for k in range(km):
        xsec_f = interp1d(wav_o3,abs_o3[k,:],kind='linear')
        abs_o3_lbl[k,:] = xsec_f(wav_abs)

    # append UV-Vis to LBL that stops at 555
    i = wav_o3 < wav_abs.min()
    all_wl = np.append(wav_abs,wav_o3[i])
    abs_o3 = np.append(abs_o3_lbl,abs_o3[:,i],axis=1)

    # Read in alphaTable
    inFile = '../../alphaTable_v0/alpha_CK_Thuillier_o3.nc4'
    nc = Dataset(inFile)
    xsec_o3 = nc.variables['xsec_o3'][:]
    o3_ch   = nc.variables['channels'][:]
    # flip so bottom to top
    xsec_o3 = xsec_o3[::-1,:]

    # convert mass mixing ratio to molecules/m3
    Na = 6.022e23 # avogadro's number
    O3_MW  = 48.0*1e-3     # kg/mole

    o3_conc = o3*AIRDENS   # kg/m3
    o3_conc = o3_conc*Na/O3_MW   # molecules/m3

    # get the optical depth subcolumns
    alpha_o3 = np.zeros(abs_o3.shape)
    for i in range(km):
        alpha_o3[i,:] = o3_conc[i]*dz[i]*abs_o3[i,:]


    sys.exit()

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



