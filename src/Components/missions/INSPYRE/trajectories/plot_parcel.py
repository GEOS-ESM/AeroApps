#!/usr/bin/env python3.11

import xarray as xr
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
from matplotlib.dates import MonthLocator, WeekdayLocator, DateFormatter
import matplotlib.ticker as ticker
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from shapely.geometry import Polygon
from cartopy.geodesic import Geodesic  # from geographiclib

filen = 'parcel_traj.dryrun_20250801_2100z.nc'
#filen = 'parcel_traj.dryrun_20250728_2030z.nc'
ds = xr.open_mfdataset(filen)
lat = ds["lat"].values
lon = ds["lon"].values
alt = ds["PAlt"].values
print(lon.shape)
ds.close()

# Make a map
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature

fig = plt.figure(figsize=(14,9))
ax = plt.axes(projection=ccrs.PlateCarree())
ax.coastlines(resolution="50m")
ax.gridlines()
ax.add_feature(cfeature.STATES, linestyle=':', edgecolor='black', linewidth=0.5)
#ax.set_xlim([-110,-20])
ax.set_xlim([-120,-60])
ax.set_ylim([20,70])
ax.stock_img()

#cs = plt.scatter(x=lon[24,:],y=lat[24,:],c=alt[24,:],vmin=100,vmax=160,marker="X",cmap="plasma")
#cs = plt.scatter(x=lon[48,:],y=lat[48,:],c=alt[24,:],vmin=100,vmax=160,marker="o",cmap="plasma")
for i in range(0,2500,25):
    cs = plt.scatter(x=lon[:,i],y=lat[:,i],c=alt[:,i],vmin=12,vmax=18,marker=".",cmap="plasma")
cbar = fig.colorbar(cs)
plt.show()
