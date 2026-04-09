#!/usr/bin/env python3

from model_utils.utils import plot_utils
import numpy as np

def closeness(baseline, perturb, sat, yy, mm,
              expdir="/home/pcolarco/geos_aerosols/pcolarco/"):
    
    basedir  = f"{expdir}/{baseline}/"
    pertdir  = f"{expdir}/{perturb}/"
    baseaod, lon, lat  = plot_utils.read_model(baseline,basedir,yy,mm,sat)
    print(perturb, pertdir)
    pertaod, lon, lat  = plot_utils.read_model(perturb ,pertdir,yy,mm,sat)
    obsaod, lono, lato = plot_utils.read_sat(sat,yy,mm)
    obsaodi = obsaod.interp_like(baseaod,method="nearest",assume_sorted=True)
    plot_utils.plotcloseness(yy, mm, lon, lat, baseaod.values, pertaod.values,
                  np.squeeze(obsaodi.values),baseline,perturb,sat)
    return

def difference(src1, src2, yy, mm,
               expdir="/home/pcolarco/geos_aerosols/pcolarco/"):

#   There's some logic here to allow src1 or src2 to be either
#   a model or satellite data

#   Get field 1
    sat1 = 0
    expdir1  = f"{expdir}/{src1}/"
    if((src1 == "MOD04") or (src1 == "MYD04")):
        sat1 = 1
        var1_, lono, lato = plot_utils.read_sat(src1,yy,mm)
    else:
        var1_, lon, lat   = plot_utils.read_model(src1,expdir1,yy,mm,sat)

#   Get field 2
    sat2 = 0
    expdir2  = f"{expdir}/{src2}/"
    if((src2 == "MOD04") or (src2 == "MYD04")):
        sat2 = 2
        var2_, lono, lato = plot_utils.read_sat(src2,yy,mm)
    else:
        var2_, lon, lat   = plot_utils.read_model(src2,expdir2,yy,mm,sat)

#   Interpolate?
    if(sat1 and sat2):
        lon = lono
        lat = lato
        var1 = var1_
        var2 = var2_
    elif(not(sat1) and not(sat2)):
        var1 = var1_
        var2 = var2_
    else:
        if(sat1):
            var1 = var1_.interp_like(var2_,method="nearest",assume_sorted=True)
            var2 = var2_
        else:
            var2 = var2_.interp_like(var1_,method="nearest",assume_sorted=True)
            var1 = var1_
    plot_utils.plotdifference(yy,mm,lon,lat,
                              np.squeeze(var1.values),np.squeeze(var2.values)
                              ,src1,src2)
    return
            

if __name__ == "__main__":
    baseline = "c180R_v11.8.1"
    perturb  = "c180R_v11.8.1_v2xx"
    sat      = "MYD04"
    yy       = 2019
    mm       = [1,2,3,4,5,6,7,8,9,10,11,12]

#   Make a closeness plot
    closeness(baseline, perturb, sat, yy, mm)

#   Make difference plots
#    difference(baseline,perturb, yy, mm)
#    difference(baseline,sat, yy, mm)
#    difference(perturb,sat, yy, mm)
