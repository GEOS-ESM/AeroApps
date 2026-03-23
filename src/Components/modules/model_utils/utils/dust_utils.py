#!/usr/bin/env python3
"""
Generic dust functions (e.g., PSD calculations)
"""

import numpy as np
import os
import sys
from scipy.special import erf
from carma_utils import findrmin, carmabins
from plot_utils import plotdifference, plotone
import xarray as xr


def vwettogwet(vwet,sandfrac=0.8):
    """
    Conversion of volumetric soil moisture fraction (0-1) to gravimetric
    following Zender et al. 2003 equations 7 - 9
    """
    rhowater = 1000.  # water density kg m-3
    rhosoil  = 2500.  # soil density kg m-3
    w_s      = 0.489 - 0.126*sandfrac
    gwet     = vwet*rhowater/(rhosoil*(1.-w_s))
    return gwet


def fecanmoisture(gwet,clayfrac=0.2):
    """
    Enhancement in soil friction threshold velocity due to soil moisture
    Parameterization follows Fe'can et al. Ann. Geophysicae 17 (1999)
    equation 15
    Input gwet = soil gravimetric moisture fraction (0-1)
    """
    wgt  = clayfrac * (14.*clayfrac + 17.) # %-gravimetric wetness threshold
    wetp = gwet*100.                       # %-gravimetric wetness from model
    fw   = np.full_like(wetp, 1.)
    fw[np.where(wetp>wgt)] = np.sqrt(1.+1.21*(wetp[np.where(wetp>wgt)]-wgt)**0.68)
    return fw


def thresholdM97(diam,rhop=2650.,rhoa=1.25):
    """
    Saltation threshold friction velocity [m s-1] for dry soils
    Equation follow Marticorena et al. JGR, 1997 Eq. 1
    Input:  diam = diameter of soil aggregate [m]
    Output: u_thresh0 = threshold friction velocity [m s-1]
    """
    grav = 9.81  # acceleration of gravity in m s-2
    u_thresh0 = 0.129 * np.sqrt(rhop*grav*diam/rhoa) \
                 * np.sqrt(1.+6.e-7/(rhop*grav*diam**2.5)) \
                 / np.sqrt(1.928*(1331.*(100.*diam)**1.56+0.38)**0.092 - 1.)
    return u_thresh0


def paggregate(diam):
    """
    Aggregate size distribution after meng2021 below
    """
    massMedianDiameterAgg = 127.0     # aggregate diameter [um]
    geometricStdDevAgg = 2.95         # geometric standard deviation [1]
    pagg = 1. / (diam * np.log(geometricStdDevAgg) * np.sqrt(2.*np.pi)) \
               *  np.exp(- 1. * ((np.log(diam)-np.log(massMedianDiameterAgg)) ** 2.) \
                         / (2.*np.log(geometricStdDevAgg)**2.))
    return pagg


def meng2021(diam, diamLow, diamUp):
    """
    Calculate probability of dust volume size distribution
    given a central diameter and range (i.e., bin width)
    based on equations 3-4 in Meng et al. 2021: 
    https://doi.org/10.1029/2021GL097287
    """
#   These are what is stated in the paper under equation 3
#    medianMassDiameter = 1.13         # median mass diameter [um]
#    geometricStdDev = 1.92            # geometric standard deviation [1]
#   These are the number from Kok 2011 and reproduce Figure 1 in Meng 2021
    medianMassDiameter = 3.4          # median mass diameter [um]
    geometricStdDev = 3.0             # geometric standard deviation [1]
    flambda = 0.15                    # Ratio of lambda to Daggregate
    geometricStdDevAgg = 2.95         # geometric standard deviation [1]
    factor = 1. / (np.sqrt(2.) * np.log(geometricStdDev))  # auxiliary constant

#   Get a PSD for the aggregates
    rlow = .01
    rup  = 50000.
    nbin = 1000
    rmin, rmrat = findrmin(nbin, rlow, rup)
    r, dr, rlow, rup, rmassup = carmabins(nbin,rmrat,rmin)
#   Integrate over the aggregate distribution (right side of eqn. 3)
    fact = 0.
    for i in range(0,nbin):
        dAggregate = 2.*r[i]
        lam = flambda*dAggregate
        fact += np.exp(-diam/lam)**3 * paggregate(dAggregate) *2.*dr[i]
    dvdlnd = diam * (1.0 + erf(factor * np.log(diam / medianMassDiameter))) * fact
    return dvdlnd


def kok2011(diam, diamLow, diamUp):
    """
    Calculate probability of dust volume size distribution
    given a central diameter and range (i.e., bin width)
    based on equation 6 in Kok et al. 2011: 
    https://doi.org/10.1073/pnas.1014798108
    """
    medianMassDiameter = 3.4          # median mass diameter [um]
    geometricStdDev = 3.0             # geometric standard deviation [1]
    crackPropagationLength = 12.0     # crack propagation length [um]
    factor = 1. / (np.sqrt(2.) * np.log(geometricStdDev))  # auxiliary constant

    
    dlam   = diam / crackPropagationLength
    dvdlnd = diam * (1.0 + erf(factor * np.log(diam / medianMassDiameter))) * np.exp(-dlam**3) * np.log(diamUp / diamLow)
    return dvdlnd

def kok2021(diam, diamLow, diamUp):
    """
    This is as in kok2011 above but the parameters are modified according
    to GOCART2G_Process.F90 code in kok2021 subroutine. The numbers are
    claimed to be from Meng 2021 but I don't see them there. This is a 
    small difference when integrated across the GOCART bins.
    """
    medianMassDiameter = 3.5          # median mass diameter [um]
    geometricStdDev = 2.8             # geometric standard deviation [1]
    crackPropagationLength = 11.0     # crack propagation length [um]
    factor = 1. / (np.sqrt(2.) * np.log(geometricStdDev))  # auxiliary constant

    
    dlam   = diam / crackPropagationLength
    dvdlnd = diam * (1.0 + erf(factor * np.log(diam / medianMassDiameter))) * np.exp(-dlam**3) * np.log(diamUp / diamLow)
    return dvdlnd

def plotsource(src,extent=[-180,180,-90,90]):
    if src == "topo":
        fname = "/home/pcolarco/ExtData/chemistry/DUST/v0.0.0/sfc/gocart.dust_source.v5a.x1152_y721.nc"
        vname = "du_src"
        ds    = xr.open_mfdataset(fname)
        src_  = ds[vname]
    elif src == "flex":
        fname = "/discover/nobackup/projects/ARCSIX/kamoore6/Useful_Files/flexdust.nc"
        vname = "soil"
        ds    = xr.open_mfdataset(fname)
        src_  = ds[vname][0]
    elif src == "ssm":
        fname = "/discover/nobackup/pcolarco/fvInput/chemistry/DUST/v0.0.0/sfc/ssm_global_2m.x10800_y5400.prc.nc"
        vname = "ssm"
        ds_   = xr.open_mfdataset(fname)
        ds    = ds_.rename({'lon':'longitude','lat':'latitude'})
        src_  = ds[vname][0]
    else:
        print("incorrect source")
        sys.exit()
    try:
        lon    = ds["lon"].values
        lat    = ds["lat"].values
    except:
        lon    = ds["longitude"].values
        lat    = ds["latitude"].values
     
    fname = f"DU_SRC.{src}.png"
    title = f"Dust Sources"
    varv = np.ma.array(np.squeeze(src_.values))
    varv = np.ma.masked_values(varv,0.)
    print(varv)
    plotone(lon, lat, varv, src,
            extent=extent,box=None,
            cbartitle="Dust Source", title=title,
            prange=[0,1],fname=fname,proj="PlateCarree")
         

def plotsourcediff(src1,src2,extent=[-180,180,-90,90],drange=[-0.2,0.2]):
#   Read source file
    for src in [src1,src2]:
        if src == "topo":
            fname = "/home/pcolarco/ExtData/chemistry/DUST/v0.0.0/sfc/gocart.dust_source.v5a.x1152_y721.nc"
            vname = "du_src"
            ds    = xr.open_mfdataset(fname)
            src_  = ds[vname]
        elif src == "flex":
            fname = "/discover/nobackup/projects/ARCSIX/kamoore6/Useful_Files/flexdust.nc"
            vname = "soil"
            ds    = xr.open_mfdataset(fname)
            src_  = ds[vname][0]
        elif src == "ssm":
            fname = "/discover/nobackup/pcolarco/fvInput/chemistry/DUST/v0.0.0/sfc/ssm_global_2m.x10800_y5400.prc.nc"
            vname = "ssm"
            ds_   = xr.open_mfdataset(fname)
            ds    = ds_.rename({'lon':'longitude','lat':'latitude'})
            src_  = ds[vname][0]
        else:
            print("incorrect source")
            sys.exit()
        if src == src1:
            src1__ = np.squeeze(src_)
            src1v  = np.squeeze(src1__.values)
            try:
                lon    = ds["lon"].values
                lat    = ds["lat"].values
            except:
                lon    = ds["longitude"].values
                lat    = ds["latitude"].values
        else:
            src2__ = src_.interp_like(src1__,method="nearest",assume_sorted=True)
            src2v  = np.squeeze(src2__.values)
    fname = f"DU_SRC.{src1}_{src2}.png"
    title = f"Dust Sources"
    plotdifference(lon, lat, src1v, src2v, src1, src2,
                   extent=extent,box=None,
                   cbartitle="Dust Source", title=title,
                   prange=[0,1],drange=drange,fname=fname,proj="PlateCarree")
        
            


if __name__ == "__main__":
    rlow = 0.1
    rup  = 10.
    nbin = 1000
    
    rmin, rmrat = findrmin(nbin, rlow, rup)
    r, dr, rlow, rup, rmassup = carmabins(nbin,rmrat,rmin)

    dvdlnd = []
    for i in np.arange(0,1000):
        d    = r[i]*2.
        dup  = rup[i]*2.
        dlow = rlow[i]*2.
#        dvdlnd.append(kok2011(d,dlow,dup))
#        dvdlnd.append(kok2021(d,dlow,dup))
        dvdlnd.append(meng2021(d,dlow,dup))

#   Normalize
    dvdlnd = dvdlnd/np.sum(dvdlnd*dr/r)

#   Calculate volume fractions at GOCART bins
    rl = [0.1,1.0,1.8,3.0,6.0]
    ru = [1.0,1.8,3.0,6.0,10.]
    frac = []
    for i in np.arange(0,5):
        dv = dvdlnd[np.where((r>=rl[i]) & (r<ru[i]))]
        f  = r[np.where((r>=rl[i]) & (r<ru[i]))]/dr[np.where((r>=rl[i]) & (r<ru[i]))]
        frac.append(np.sum(dv*f))

    print(frac)

#    print(thresholdM97(100.e-6))
#    plotsourcediff("topo","flex",extent=[-20,60,0,40])
#    plotsourcediff("topo","ssm",extent=[-20,60,0,40],drange=[-0.5,0.5])
#    plotsource("topo",extent=[-20,60,0,40])
#    plotsource("ssm",extent=[-180,180,-60,60])
#    plotsource("flex",extent=[-20,60,0,40])
