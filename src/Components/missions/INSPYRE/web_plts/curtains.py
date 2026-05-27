import numpy as np
from netCDF4 import Dataset
import matplotlib.pyplot as plt
from matplotlib import cm #cm is colormap
import matplotlib.colors



def plot(x, xtitle, lev, var, vmax, title, figure, yyyymmdd, hh, yyyymmddt):
    fig = plt.figure(1, figsize=(12,7))
    ax1 = fig.add_axes([0.1, 0.1, 0.9, 0.8])
    #ax1.set_yscale("log")
    plt.ylim(1000,100)
    oc_lev = np.arange(0,vmax,vmax/50)
    im2 = ax1.contourf(x, lev, var*1e9,levels=oc_lev,extend = 'both',cmap='Spectral_r')
    #im2.cmap.set_over('m')
    im2.cmap.set_under('w')
    plt.xlabel(xtitle)
    plt.ylabel('Pressure [hPa]')
    plt.yticks(ticks=[1000,900,800,700,600,500,400,300,200,100],labels=['1000','900','800','700','600','500','400','300','200','100'])
    cbar2=plt.colorbar(im2, ax=ax1, orientation='vertical', ticks=[0,1,2,3,4,5])
    #cbar2.ax.set_yticks([0,1,2,3,4,5])
    #cbar2.ax.set_yticklabels(['0','1','2','3','4','>5'])
    ax1.set_title(title)
    plt.text(76,150,'Forecast: %s_%sz'%(yyyymmdd,hh),color='w')
    plt.text(76,175,'Valid:       %s_12z'%(yyyymmddt),color='w')
    fig.savefig(figure)
    plt.figure().clear()
    plt.close()
    plt.cla()
    plt.clf()

if __name__ == "__main__":
   yyyy = '2025'
   mm   = '07'
   dd   = '26'
   hh   = '00'
   yyyymmdd = '%s%s%s'%(yyyy,mm,dd)

   yyyyt = '2025'
   mmt   = '07'
   ddt   = '27'
   yyyymmddt = '%s%s%s'%(yyyyt,mmt,ddt)

# READ THE DATA

   filen = '/home/pcolarco/fp/forecast/Y%s/M%s/D%s/H%s/GEOS.fp.fcst.inst3_3d_aer_Np.%s_%s+%s_1200.V01.nc4'%(yyyy,mm,dd,hh,yyyymmdd,hh,yyyymmddt)
   with Dataset(filen,'r') as ncid:
    lon  = ncid.variables['lon'][:] # longitude points
    lat  = ncid.variables['lat'][:] # latitude points
    lev  = ncid.variables['lev'][:] # Altitude in meters
    time = ncid.variables['time'][:] #time inseconds and UTC on the day of flight
    oc   = np.squeeze(ncid.variables['OC'][:,:,:,:])
    du   = np.squeeze(ncid.variables['DU'][:,:,:,:])

# SUBSETS
# -------
# Alert to NP
   j = np.where((lat >= 30.) & (lat <= 70))[0]
   i = np.where((lon >= -110.) & (lon <= -110.))
   oc2  = np.squeeze(oc[:,j,i])
   du2  = np.squeeze(du[:,j,i])
   lat2 = lat[j]

# Greenland to NP
   j = np.where((lat >= 30.) & (lat <= 70))[0]
   i = np.where((lon >= -120.) & (lon <= -120.))
   oc3  = np.squeeze(oc[:,j,i])
   du3  = np.squeeze(du[:,j,i])
   lat3 = lat[j]

# Svalbard to NP
   j = np.where((lat >= 30.) & (lat <= 70))[0]
   i = np.where((lon >= -100.) & (lon <= -100.))
   oc4  = np.squeeze(oc[:,j,i])
   du4  = np.squeeze(du[:,j,i])
   lat4 = lat[j]

# EW Alert to Svalbard
   j = np.where((lat >= 40.) & (lat <= 40.1))[0]
   i = np.where((lon >= -120.) & (lon <= -80))
   oc5  = np.squeeze(oc[:,j,i])
   du5  = np.squeeze(du[:,j,i])
   lon5 = lon[i]

# FIGURES
# -------
# Alert to NP
   plot(lat2, 'Latitude', lev, oc2, 5, 
        'Forecast Organic Carbon Mixing Ratio [ppbm] @ 62.3$^\circ$ W (Alert to North Pole)',
	'alert.forecast_OC.%s_%sz.png'%(yyyymmdd,hh),
	yyyymmdd, hh, yyyymmddt)

# Greenland to NP
   plot(lat3, 'Latitude', lev, oc3, 5, 
        'Forecast Organic Carbon Mixing Ratio [ppbm] @ 60$^\circ$ W (Greenland to North Pole)',
	'greenland.forecast_OC.%s_%sz.png'%(yyyymmdd,hh),
	yyyymmdd, hh, yyyymmddt)

# Svalbard to NP
   plot(lat4, 'Latitude', lev, oc4, 5, 
	'Forecast Organic Carbon Mixing Ratio [ppbm] @ 15.5$^\circ$ E (Svalbard to North Pole)',
	'svalbard.forecast_OC.%s_%sz.png'%(yyyymmdd,hh),
	yyyymmdd, hh, yyyymmddt)

# Alert to Svalbard
   plot(lon5, 'Longitude', lev, oc5, 5, 
	'Forecast Organic Carbon Mixing Ratio [ppbm] @ 80$^\circ$ N (Alert to Svalbard)',
	'alert_svalbard.forecast_OC.%s_%sz.png'%(yyyymmdd,hh),
	yyyymmdd, hh, yyyymmddt)

# DUST
# Alert to NP
   plot(lat2, 'Latitude', lev, du2, 50, 
        'Forecast Dust Mixing Ratio [ppbm] @ 62.3$^\circ$ W (Alert to North Pole)',
	'alert.forecast_DU.%s_%sz.png'%(yyyymmdd,hh),
	yyyymmdd, hh, yyyymmddt)

# Greenland to NP
   plot(lat3, 'Latitude', lev, du3, 50, 
        'Forecast Dust Mixing Ratio [ppbm] @ 60$^\circ$ W (Greenland to North Pole)',
	'greenland.forecast_DU.%s_%sz.png'%(yyyymmdd,hh),
	yyyymmdd, hh, yyyymmddt)

# Svalbard to NP
   plot(lat4, 'Latitude', lev, du4, 50, 
	'Forecast Dust Mixing Ratio [ppbm] @ 15.5$^\circ$ E (Svalbard to North Pole)',
	'svalbard.forecast_DU.%s_%sz.png'%(yyyymmdd,hh),
	yyyymmdd, hh, yyyymmddt)

# Alert to Svalbard
   plot(lon5, 'Longitude', lev, du5, 50, 
	'Forecast Dust Mixing Ratio [ppbm] @ 80$^\circ$ N (Alert to Svalbard)',
	'alert_svalbard.forecast_DU.%s_%sz.png'%(yyyymmdd,hh),
	yyyymmdd, hh, yyyymmddt)
