#!/usr/bin/env python3
import numpy as np
import os
import sys
import xarray as xr
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from closeness import read_model, read_sat, plotdifference
from carma_utils import findrmin, carmabins
from dust_utils import thresholdM97, vwettogwet, fecanmoisture, kok2011, meng2021


def plotthreshold():
    rlow = .1e-6
    rup  = 500.e-6
    nbin = 1000
    rmin, rmrat = findrmin(nbin, rlow, rup)
    r, dr, rlow, rup, rmassup = carmabins(nbin,rmrat,rmin)

    ut = []
    for i in range(0,nbin):
        ut.append(thresholdM97(2.*r[i]))
    fig, ax = plt.subplots(figsize=(8,6))
    ax.set_xscale("log")
    ax.set_xlim(1,1000)
    ax.set_ylim(0,4)
    ax.set_xlabel("Diameter [um]")
    ax.set_ylabel("Threshold Friction Speed [m s-1]")
    #ax.set_yscale("log")
    ax.plot(2.*r*1e6,ut)
    plt.show()
    return


def plotmoisture():
    vwet = np.arange(0,1,.01)
    gwet = vwettogwet(vwet)
    print(gwet)
    fw   = fecanmoisture(gwet)
    print(fw)
    fig, ax = plt.subplots(figsize=(8,6))
    ax.set_xlim(0,20)
    ax.set_ylim(0,4)
    ax.set_xlabel("gravimetric soil moisture (%)")
    ax.set_ylabel("erosion threshold velocity ratio")
    ax.plot(100.*gwet,fw)
    plt.show()
    return

def plotemittedsize():
    rlow = 0.05
    rup  = 20.
    nbin = 200
    
    rmin, rmrat = findrmin(nbin, rlow, rup)
    r, dr, rlow, rup, rmassup = carmabins(nbin,rmrat,rmin)

    dvdlndk = []
    dvdlndm = []
    print(r*2.)
    for i in np.arange(0,nbin):
        d    = r[i]*2.
        dup  = rup[i]*2.
        dlow = rlow[i]*2.
        dvdlndk.append(kok2011(d,dlow,dup))
        dvdlndm.append(meng2021(d,dlow,dup))

#   Normalize
    dvdlndk = dvdlndk/np.max(dvdlndk)
    dvdlndk[np.where(r>10.)] = 0.
    dvdlndm = dvdlndm/np.max(dvdlndm)
    a = np.sum(dvdlndk*dr/r)
    b = np.sum(dvdlndm[np.where(r<10)]*dr[np.where(r<10)]/r[np.where(r<10)])
    dvdlndm = dvdlndm/b
    print(np.sum(dvdlndk*dr/r), np.sum(dvdlndm*dr/r))
    print(np.sum(dvdlndm[np.where(r<10)]*dr[np.where(r<10)]/r[np.where(r<10)]))
#    dvdlndm = dvdlndm/np.sum(dvdlndm*dr/r)
    
    fig, ax = plt.subplots(figsize=(8,6))
    ax.set_xlim(.2,40)
    ax.set_ylim(.001,3)
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("diameter [um]")
    ax.set_ylabel("normalized dvdlnd")
    ax.plot(2.*r[np.where(r<10)],1.3*dvdlndk[np.where(r<10.)])
    ax.plot(2.*r,1.3*dvdlndm)
    plt.show()
    return


if __name__ == "__main__":
#    plotthreshold()
#    plotmoisture()
    plotemittedsize()
    sys.exit()


    baseline = "c180R_v11.8.0_develop"
    perturb  = "c180R_v11.8.0_kok"
    sat      = "MYD04"
    yy       = 2019
    mm       = [1,2,3,4,5,6,7,8,9,10,11,12]
    basedir  = f"/home/pcolarco/geos_aerosols/pcolarco/{baseline}/"
    pertdir  = f"/home/pcolarco/geos_aerosols/pcolarco/{perturb}/"

    baseaod, lon, lat  = read_model(baseline,basedir,yy,mm,sat)
    pertaod, lon, lat  = read_model(perturb,pertdir,yy,mm,sat)
    baseaoddu, lon, lat = read_model(baseline,basedir,yy,mm,sat,varn="DUEXTTAU550")
    obsaod, lono, lato = read_sat(sat,yy,mm)

#    baseaod = baseaod-0.3*baseaoddu
    
    obsaodi = obsaod.interp_like(baseaod,method="nearest",assume_sorted=True)
    extent=[-40,40,0,40]
#    plotdifference(yy, mm, lon, lat, baseaod.values,baseline,np.squeeze(obsaodi.values),sat,extent=extent)
    plotdifference(yy, mm, lon, lat, pertaod.values,perturb,np.squeeze(obsaodi.values),sat,extent=extent)
