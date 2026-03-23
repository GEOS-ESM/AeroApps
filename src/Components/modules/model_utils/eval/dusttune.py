#!/usr/bin/env python3
import numpy as np
import os
import sys
import xarray as xr
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from closeness import read_model, read_sat, plotdifference, plotcloseness, plotone

def tuneup(baseline, yy=2019, mm=[1,2,3,4,5,6,7,8,9,10,11,12], extent=[-20,40,10,30]):
    basedir  = f"/home/pcolarco/geos_aerosols/pcolarco/{baseline}/"
#   Find mean in a region
    latmx = extent[3]
    latmn = extent[2]
    lonmx = extent[1]
    lonmn = extent[0]

    obsmean = []
    minf    = []

    for sat in ["MYD04","MOD04"]:
        for m in mm:
            baseaod, lon, lat  = read_model(baseline,basedir,yy,[m],sat)
            baseaoddu, lon, lat = read_model(baseline,basedir,yy,[m],sat,varn="DUEXTTAU550")
            obsaod, lono, lato = read_sat(sat,yy,[m])

            baseafr = baseaod.where((baseaod.lat >= latmn) & (baseaod.lat <= latmx) &
                                    (baseaod.lon >=lonmn) & (baseaod.lon <= lonmx))
            baseduafr = baseaoddu.where((baseaoddu.lat >= latmn) & (baseaoddu.lat <= latmx) &
                                        (baseaoddu.lon >=lonmn) & (baseaoddu.lon <= lonmx))
            obsafr = obsaod.where((obsaod.lat >= latmn) & (obsaod.lat <= latmx) &
                                  (obsaod.lon >=lonmn) & (obsaod.lon <= lonmx))

            obsafr_ = obsafr.mean()
            obsmean.append(obsafr_.values)

            totmean = []
            minscale = 0.4
            for f in np.arange(minscale,2.0,0.01):
                tot = baseafr+(f-1.)*baseduafr
                tot_ = obsafr_ - tot.mean()
                totmean.append(np.abs(tot_.values))
            minf.append(np.argmin(totmean))

    for j in mm:
        i = j-mm[0]+1
        lmm = len(mm)
#        print(f"{j:02d}  {obsmean[i-1]:.3f}  {minf[i-1]*0.01+minscale:.2f}  {obsmean[i-1+lmm]:.3f}  {minf[i-1+lmm]*0.01+minscale:.2f}")
        print(f"{j:02d}  {obsmean[i-1]:.3f}  {minf[i-1]*0.01+minscale:.2f}  {obsmean[i-1+lmm]:.3f}  {minf[i-1+lmm]*0.01+minscale:.2f}")
    return f

def plottuneddust(baseline,yy,mm,sat="MYD04",fscale=1.0,extent=[-180,180,-90,90],
                  box=None,diffrange=[-0.2,0.2]):
    basedir  = f"/home/pcolarco/geos_aerosols/pcolarco/{baseline}/"
    baseaod, lon, lat  = read_model(baseline,basedir,yy,mm,sat)
    baseaoddu, lon, lat  = read_model(baseline,basedir,yy,mm,sat,varn="DUEXTTAU550")
    obsaod, lono, lato = read_sat(sat,yy,mm)
    obsaodi = obsaod.interp_like(baseaod,method="nearest",assume_sorted=True)
    plotdifference(yy,mm, lon, lat, baseaod+(fscale-1.)*baseaoddu,baseline,
                   np.squeeze(obsaodi),sat,varn="TOTEXTTAU550",extent=extent,
                   box=box,diffrange=diffrange)

def plotmodis(sat,yy,mm,extent=[-180,180,-90,90],
              box=None,weighted=False):
    obsaod, lon, lat = read_sat(sat,yy,mm,weighted=weighted)
    obscnt, lon, lat = read_sat(sat,yy,mm, varn="count_tau_",weighted=weighted)
    satn = sat
    if weighted:
        satn = f"{sat}w"
    plotone(yy,mm, lon, lat, np.squeeze(obsaod.values),satn,
            varn="TOTEXTTAU550",extent=extent,
            box=[-5,3,20,28], var2=np.squeeze(obscnt.values))

def plotdiffdustmod(baseline,perturb,yy,mm,sat="MYD04",fscale1=1.0,fscale2=1.0,
                    extent=[-180,180,-90,90],box=None,diffrange=[-0.1,0.1],
                    varn="TOTEXTTAU550",proj="PlateCarree"):
    basedir  = f"/home/pcolarco/geos_aerosols/pcolarco/{baseline}/"
    baseaod, lon, lat  = read_model(baseline,basedir,yy,mm,sat,varn=varn)
    pertdir  = f"/home/pcolarco/geos_aerosols/pcolarco/{perturb}/"
    pertaod, lon, lat  = read_model(perturb,pertdir,yy,mm,sat,varn=varn)
    if varn == "TOTEXTTAU550":
        baseaoddu, lon, lat  = read_model(baseline,basedir,yy,mm,sat,varn="DUEXTTAU550")
        pertaoddu, lon, lat  = read_model(perturb,pertdir,yy,mm,sat,varn="DUEXTTAU550")
        var1 = baseaod+(fscale1-1.)*baseaoddu
        var2 = pertaod+(fscale2-1.)*pertaoddu
    else:
        var1 = fscale1*baseaod
        var2 = fscale2*pertaod
    plotdifference(yy,mm, lon, lat, var1,baseline,
                   var2,perturb,proj=proj,
                   varn=varn,extent=extent,box=box,diffrange=diffrange)

def plotscalecloseness(baseline,perturb,yy,mm,sat="MYD04",fscale1=1.0,fscale2=1.0,
                       extent=[-180,180,-90,90],box=None):
    basedir  = f"/home/pcolarco/geos_aerosols/pcolarco/{baseline}/"
    baseaod, lon, lat  = read_model(baseline,basedir,yy,mm,sat)
    baseaoddu, lon, lat  = read_model(baseline,basedir,yy,mm,sat,varn="DUEXTTAU550")
    pertdir  = f"/home/pcolarco/geos_aerosols/pcolarco/{perturb}/"
    pertaod, loni, lati  = read_model(perturb,pertdir,yy,mm,sat)
    pertaoddu, loni, lati  = read_model(perturb,pertdir,yy,mm,sat,varn="DUEXTTAU550")
    obsaod, lono, lato = read_sat(sat,yy,mm)
    obsaodi = obsaod.interp_like(baseaod,method="nearest",assume_sorted=True)
    a = baseaod+(fscale1-1.)*baseaoddu
    b_ = pertaod+(fscale2-1.)*pertaoddu
    if len(lon) == len(loni):
        b = b_
    else:
        b = b_.interp_like(a,method="nearest",assume_sorted=True)
    
    plotcloseness(yy, mm, lon, lat, a, b,
                  np.squeeze(obsaodi), baseline, perturb, sat,extent=extent,box=box)


def plotverticalpsd(expid,yy=2019,mm=6,dd=5,hh=0,fscale=1.0,lonw=23.21,latw=19.23):
    filen = f"/home/pcolarco/geos_aerosols/pcolarco/{expid}/holding/inst3d_aer_v/{yy}{mm:02d}/{expid}.inst3d_aer_v.{yy}{mm:02d}{dd:02d}_{hh:02d}00z.nc4"
    ds  = xr.open_mfdataset(filen)
    rho = np.squeeze(ds["AIRDENS"].sel(lon=lonw,lat=latw,method="nearest").values)
    h   = np.squeeze(ds["H"].sel(lon=lonw,lat=latw,method="nearest").values)/1000.
    du1 = np.squeeze(ds["DU001"].sel(lon=lonw,lat=latw,method="nearest").values)*rho
    du2 = np.squeeze(ds["DU002"].sel(lon=lonw,lat=latw,method="nearest").values)*rho
    du3 = np.squeeze(ds["DU003"].sel(lon=lonw,lat=latw,method="nearest").values)*rho
    du4 = np.squeeze(ds["DU004"].sel(lon=lonw,lat=latw,method="nearest").values)*rho
    du5 = np.squeeze(ds["DU005"].sel(lon=lonw,lat=latw,method="nearest").values)*rho
    psd = np.transpose(np.squeeze(np.array([[du1],[du2],[du3],[du4],[du5]])))*fscale*1e9
    fig, ax = plt.subplots(figsize=(8,8))
    clevs=np.arange(0,200,10)
    ax.set_xlim(0,6)
    ax.set_ylim(0,6)
    ax.set_xlabel("bin #")
    ax.set_ylabel("altitude [km]")
    im = ax.contourf(np.arange(1,6),h,psd,levels=clevs,cmap='Spectral_r',extend='max')
    cbar = fig.colorbar(im, orientation='horizontal', extend='max')
    cbar.set_label("Concentration [ug m-3]")
    plt.suptitle(f"{expid}: {yy:04d}{mm:02d}{dd:02d}_{hh:02d}", fontsize=14, y=0.95)
    plt.savefig(f"DUCONC.{expid}.{yy:04d}{mm:02d}{dd:02d}_{hh:02d}z.png")
    return


if __name__ == "__main__":
    """
#   Outputs a table of monthly scaling factors f that minimizes
#   the regional AOD difference in the defined box between
#   monthly MYD04/MOD04 and the model. The f factor is the 
#   amount to scale the dust part of the total AOD to achieve
#   this minimization.
#   Erg Chech
    box = [-5,3,20,28]
#   Mauritania
#    box = [-15,-10,15,20]
#   Bodele
#    box = [10,22,12,20]
#    fscale = tuneup("c180R_v11.8.0_develop",extent=box)
#    fscale = tuneup("c180R_v11.8.0_kok",extent=box)
#    fscale = tuneup("c180R_v11.8.0_newdust",extent=box)
#    fscale = tuneup("c180R_v11.8.0_dead",extent=box)
#    fscale = tuneup("c180R_v11.8.0_deadssm",extent=box)
#    fscale = tuneup("c180R_v11.8.0_deadflex",extent=box)
#    fscale = tuneup("c180R_v11.8.0_k14",extent=box)
#    fscale = tuneup("c180R_v11.8.0_kok_p2650",extent=box)
#    fscale = tuneup("c180R_v11.8.0_kok_p1000",extent=box)
#    fscale = tuneup("c180R_v11.8.0_kok_p250",extent=box)
#    fscale = tuneup("c180R_v11.8.0_develop_ssm",extent=box)
#    fscale = tuneup("c180R_v11.8.0_develop_flex",extent=box)
#    fscale = tuneup("c360R_v11.8.0_kok_p2650",extent=box,mm=[5,6])
#    sys.exit()
    """

    """
#   A series of plots that compare model to MYD04 or model to model
#   as tuned (by hand) from the scaling factors determined above.
    box = [-5,3,20,28]
    yy       = 2019
    mm       = [6]
    plotdiffdustmod("c180R_v11.8.0_develop","c180R_v11.8.0_kok",yy,mm,
                    fscale1=1.01,fscale2=1.00,extent=[-40,40,0,40],box=box)
    plotdiffdustmod("c180R_v11.8.0_develop","c180R_v11.8.0_dead",yy,mm,
                    fscale1=1.01,fscale2=0.94,extent=[-40,40,0,40],box=box)
    plotdiffdustmod("c180R_v11.8.0_develop","c180R_v11.8.0_deadssm",yy,mm,
                    fscale1=1.01,fscale2=1.08,extent=[-40,40,0,40],box=box)
    plotdiffdustmod("c180R_v11.8.0_develop","c180R_v11.8.0_deadflex",yy,mm,
                    fscale1=1.01,fscale2=1.01,extent=[-40,40,0,40],box=box)
    plotdiffdustmod("c180R_v11.8.0_develop","c180R_v11.8.0_k14",yy,mm,
                    fscale1=1.01,fscale2=1.03,extent=[-40,40,0,40],box=box)
    plotdiffdustmod("c180R_v11.8.0_kok","c180R_v11.8.0_kok_p2650",yy,mm,
                    fscale1=1.0,fscale2=1.0,extent=[-40,40,0,40],box=box,diffrange=[-0.05,0.05])
    plotdiffdustmod("c180R_v11.8.0_kok","c180R_v11.8.0_kok_p1000",yy,mm,
                    fscale1=1.0,fscale2=0.88,extent=[-40,40,0,40],box=box,diffrange=[-0.05,0.05])
    plotdiffdustmod("c180R_v11.8.0_kok","c180R_v11.8.0_kok_p250",yy,mm,
                    fscale1=1.0,fscale2=0.86,extent=[-40,40,0,40],box=box,diffrange=[-0.05,0.05])
    plottuneddust("c180R_v11.8.0_develop",yy,mm,fscale=1.01,extent=[-40,40,0,40],box=box)
    plottuneddust("c180R_v11.8.0_kok",yy,mm,fscale=1.0,extent=[-40,40,0,40],box=box)
    plottuneddust("c180R_v11.8.0_dead",yy,mm,fscale=0.94,extent=[-40,40,0,40],box=box)
    plottuneddust("c180R_v11.8.0_deadssm",yy,mm,fscale=1.08,extent=[-40,40,0,40],box=box)
    plottuneddust("c180R_v11.8.0_deadflex",yy,mm,fscale=1.01,extent=[-40,40,0,40],box=box)
    plottuneddust("c180R_v11.8.0_k14",yy,mm,fscale=1.03,extent=[-40,40,0,40],box=box)
    plottuneddust("c360R_v11.8.0_kok_p2650",yy,mm,fscale=0.94,extent=[-40,40,0,40],box=box)
    plottuneddust("c180R_v11.8.0_kok_p2650",yy,mm,fscale=1.00,extent=[-40,40,0,40],box=box)
    plottuneddust("c180R_v11.8.0_kok_p1000",yy,mm,fscale=0.88,extent=[-40,40,0,40],box=box)
    plottuneddust("c180R_v11.8.0_kok_p250",yy,mm,fscale=0.86,extent=[-40,40,0,40],box=box)
    sys.exit()
    """

    """
#   Plot some fun maps of the MODIS AOD/Counts and MODIS AOD/Dust source functions
    yy       = 2019
    mm       = [6]

    plotmodis("MYD04",yy,mm,extent=[-40,40,0,40])#,weighted=True)
    
#   Topo source
    filename = "/home/pcolarco/ExtData/chemistry/DUST/v0.0.0/sfc/gocart.dust_source.v5a.x1152_y721.nc"
    ds = xr.open_mfdataset(filename)
    du_src = ds["du_src"].values
    sat = "MYD04"
    obsaod, lon, lat = read_sat(sat,yy,mm)
    satn = f"{sat}.du_src"
    extent=[-40,40,0,40]
    plotone(yy,mm, lon, lat, np.squeeze(obsaod.values),satn,
            varn="TOTEXTTAU550",extent=extent,
            box=[-5,3,20,28], var2=np.squeeze(du_src))
#   SSM
    filename = "/discover/nobackup/pcolarco/fvInput/chemistry/DUST/v0.0.0/sfc/ssm_global_2m.x10800_y5400.prc.nc"
    ds = xr.open_mfdataset(filename)
    sat = "MYD04"
    obsaod, lon, lat = read_sat(sat,yy,mm)
    ssmi = ds["ssm"].interp_like(obsaod,method="nearest", assume_sorted=True)
    ssm = ssmi.values
    satn = f"{sat}.ssm"
    extent=[-40,40,0,40]
    plotone(yy,mm, lon, lat, np.squeeze(obsaod.values),satn,
            varn="TOTEXTTAU550",extent=extent,
            box=[-5,3,20,28], var2=np.squeeze(ssm))

#   FLEX
    filename = "/discover/nobackup/projects/ARCSIX/kamoore6/Useful_Files/flexdust.nc"
    ds = xr.open_mfdataset(filename)
    sat = "MYD04"
    obsaod, lon, lat = read_sat(sat,yy,mm)
    flex = ds["soil"].values
    satn = f"{sat}.flex"
    extent=[-40,40,0,40]
    plotone(yy,mm, lon, lat, np.squeeze(obsaod.values),satn,
            varn="TOTEXTTAU550",extent=extent,
            box=[-5,3,20,28], var2=np.squeeze(flex[0,:,:]))
    """

    """
#   Plot some Closeness maps of the pairs of model experiments
#   using scaled dust adjustments as above
    yy = 2019
    mm = [1,2,3,4,5,6,7,8,9,10,11,12]
    extent=[-20,120,-10,60]
    extent=[-180,180,-90,90]
    box=[-5,3,20,28]
    plotscalecloseness("c180R_v11.8.0_develop","c180R_v11.8.0_kok",yy,mm,
                       fscale1=1.01,fscale2=1.0,extent=extent,box=box)
    plotscalecloseness("c180R_v11.8.0_develop","c180R_v11.8.0_dead",yy,mm,
                       fscale1=1.01,fscale2=0.94,extent=extent,box=box)
    plotscalecloseness("c180R_v11.8.0_develop","c180R_v11.8.0_deadssm",yy,mm,
                       fscale1=1.01,fscale2=1.08,extent=extent,box=box)
    plotscalecloseness("c180R_v11.8.0_develop","c180R_v11.8.0_deadflex",yy,mm,
                       fscale1=1.01,fscale2=1.01,extent=extent,box=box)
    plotscalecloseness("c180R_v11.8.0_develop","c180R_v11.8.0_k14",yy,mm,
                       fscale1=1.01,fscale2=1.03,extent=extent,box=box)
    plotscalecloseness("c180R_v11.8.0_kok_p2650","c180R_v11.8.0_kok_p1000",yy,mm,
                       fscale1=1.0,fscale2=0.88,extent=extent,box=box)
    plotscalecloseness("c180R_v11.8.0_kok_p2650","c180R_v11.8.0_kok_p250",yy,mm,
                       fscale1=1.0,fscale2=0.86,extent=extent,box=box)

    mm = [6]
    extent=[-40,40,0,40]
    plotscalecloseness("c180R_v11.8.0_develop","c180R_v11.8.0_kok",yy,mm,
                       fscale1=1.01,fscale2=1.0,extent=extent,box=box)
    plotscalecloseness("c180R_v11.8.0_develop","c180R_v11.8.0_dead",yy,mm,
                       fscale1=1.01,fscale2=0.94,extent=extent,box=box)
    plotscalecloseness("c180R_v11.8.0_develop","c180R_v11.8.0_deadssm",yy,mm,
                       fscale1=1.01,fscale2=1.08,extent=extent,box=box)
    plotscalecloseness("c180R_v11.8.0_develop","c180R_v11.8.0_deadflex",yy,mm,
                       fscale1=1.01,fscale2=1.01,extent=extent,box=box)
    plotscalecloseness("c180R_v11.8.0_develop","c180R_v11.8.0_k14",yy,mm,
                       fscale1=1.01,fscale2=1.03,extent=extent,box=box)
    plotscalecloseness("c180R_v11.8.0_kok_p2650","c180R_v11.8.0_kok_p1000",yy,mm,
                       fscale1=1.0,fscale2=0.88,extent=extent,box=box)
    plotscalecloseness("c180R_v11.8.0_kok_p2650","c180R_v11.8.0_kok_p250",yy,mm,
                       fscale1=1.0,fscale2=0.86,extent=extent,box=box)
    plotscalecloseness("c180R_v11.8.0_kok_p2650","c360R_v11.8.0_kok_p2650",yy,mm,
                       fscale1=1.0,fscale2=0.86,extent=extent,box=box)
    
    """

    """
    plotverticalpsd("c360R_v11.8.0_kok_p2650",fscale=0.86)
    plotverticalpsd("c180R_v11.8.0_develop",fscale=1.01)
    plotverticalpsd("c180R_v11.8.0_kok",fscale=1.0)
    plotverticalpsd("c180R_v11.8.0_kok_p2650",fscale=1.0)
    plotverticalpsd("c180R_v11.8.0_kok_p1000",fscale=0.88)
    plotverticalpsd("c180R_v11.8.0_kok_p250",fscale=0.86)
    plotverticalpsd("c180R_v11.8.0_dead",fscale=0.94)
    plotverticalpsd("c180R_v11.8.0_deadflex",fscale=1.01)
    plotverticalpsd("c180R_v11.8.0_deadssm",fscale=1.08)
    plotverticalpsd("c180R_v11.8.0_k14",fscale=1.03)
    """


    """
    plotverticalpsd("c180R_v11.8.0_develop",fscale=1.01,hh=6,lonw=-52.07,latw=13.08)
    plotverticalpsd("c180R_v11.8.0_kok",fscale=1.0,hh=6,lonw=-52.07,latw=13.08)
    plotverticalpsd("c180R_v11.8.0_kok_p2650",fscale=1.0,hh=6,lonw=-52.07,latw=13.08)
    plotverticalpsd("c180R_v11.8.0_kok_p1000",fscale=0.88,hh=6,lonw=-52.07,latw=13.08)
    plotverticalpsd("c180R_v11.8.0_kok_p250",fscale=0.86,hh=6,lonw=-52.07,latw=13.08)
    """

    """
#   A series of plots that compare model to model
#   as tuned (by hand) from the scaling factors determined above.
    extent=[-90,-10,0,40]
    box = [-5,3,20,28]
    diffrange=[-0.2,0.2]
    yy       = 2019
    mm       = [6]
    plotdiffdustmod("c180R_v11.8.0_kok","c180R_v11.8.0_kok_p2650",yy,mm,varn="DUCMASS",
                    fscale1=1.0*1e3,fscale2=1.0*1e3,extent=extent,box=box,diffrange=diffrange)
    plotdiffdustmod("c180R_v11.8.0_kok_p2650","c180R_v11.8.0_kok_p1000",yy,mm,varn="DUCMASS",
                    fscale1=1.0*1e3,fscale2=0.88*1e3,extent=extent,box=box,diffrange=diffrange)
    plotdiffdustmod("c180R_v11.8.0_kok_p2650","c180R_v11.8.0_kok_p250",yy,mm,varn="DUCMASS",
                    fscale1=1.0*1e3,fscale2=0.86*1e3,extent=extent,box=box,diffrange=diffrange)
    diffrange=[-0.05,0.05]
    plotdiffdustmod("c180R_v11.8.0_kok","c180R_v11.8.0_kok_p2650",yy,mm,
                    fscale1=1.0,fscale2=1.0,extent=extent,box=box,diffrange=diffrange)
    plotdiffdustmod("c180R_v11.8.0_kok_p2650","c180R_v11.8.0_kok_p1000",yy,mm,
                    fscale1=1.0,fscale2=0.88,extent=extent,box=box,diffrange=diffrange)
    plotdiffdustmod("c180R_v11.8.0_kok_p2650","c180R_v11.8.0_kok_p250",yy,mm,
                    fscale1=1.0,fscale2=0.86,extent=extent,box=box,diffrange=diffrange)
     
    """
    extent=[-90,-10,0,40]
    box = [-5,3,20,28]
    diffrange=[-0.005,0.005]
    yy       = 2019
    mm       = [8]
#    plotdiffdustmod("c180R_v11.8.0_kok_p2650","c180R_v11.8.0_k14",yy,mm,proj="nps",
#                    fscale1=1.0,fscale2=1.03,extent=extent,box=box,diffrange=diffrange)
#    plotdiffdustmod("c180R_v11.8.0_kok_p2650","c180R_v11.8.0_deadssm",yy,mm,proj="nps",
#                    fscale1=1.0,fscale2=1.08,extent=extent,box=box,diffrange=diffrange)
#    plotdiffdustmod("c180R_v11.8.0_kok_p2650","c180R_v11.8.0_deadflex",yy,mm,proj="nps",
#                    fscale1=1.0,fscale2=1.01,extent=extent,box=box,diffrange=diffrange)

#    plotdiffdustmod("c180R_v11.8.0_kok_p2650","c180R_v11.8.0_k14",yy,mm,varn="DUCMASS",proj="nps",
#                    fscale1=1.0*1e3,fscale2=1.03*1e3,extent=extent,box=box,diffrange=diffrange)
#    plotdiffdustmod("c180R_v11.8.0_kok_p2650","c180R_v11.8.0_deadssm",yy,mm,varn="DUCMASS",proj="nps",
#                    fscale1=1.0*1e3,fscale2=1.08*1e3,extent=extent,box=box,diffrange=diffrange)
#    plotdiffdustmod("c180R_v11.8.0_kok_p2650","c180R_v11.8.0_deadflex",yy,mm,varn="DUCMASS",proj="nps",
#                    fscale1=1.0*1e3,fscale2=1.01*1e3,extent=extent,box=box,diffrange=diffrange)

    plotdiffdustmod("c180R_v11.8.0_develop","c180R_v11.8.0_develop_ssm",yy,mm,
                    fscale1=1.01,fscale2=0.46,extent=[-40,40,0,40],box=box)
    plotdiffdustmod("c180R_v11.8.0_develop","c180R_v11.8.0_develop_flex",yy,mm,
                    fscale1=1.01,fscale2=1.43,extent=[-40,40,0,40],box=box)
    plotscalecloseness("c180R_v11.8.0_develop","c180R_v11.8.0_develop_ssm",yy,mm,
                       fscale1=1.01,fscale2=0.46,extent=extent,box=box)
    plotscalecloseness("c180R_v11.8.0_develop","c180R_v11.8.0_develop_flex",yy,mm,
                       fscale1=1.01,fscale2=1.43,extent=extent,box=box)
