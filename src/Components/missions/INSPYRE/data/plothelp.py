from pyobs.icartt import ICARTT
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.path as mpath
from matplotlib import cm #cm is colormap
import matplotlib.dates as mdates
import matplotlib.colors as colors
import matplotlib.ticker as mticker
import matplotlib
matplotlib.use('agg')

def plot_timeseries(ictFile):

#   Get the ICARTT file describing the trajectory
    if ictFile.find('ER2') > 0:
        aircraft = 'ER2'
        i0 = ictFile.find('ER2')+4
    if ictFile.find('GV')  > 0:
        aircraft = 'GV'
        i0 = ictFile.find('GV')+3
    if ictFile.find('Learjet') > 0:
        aircraft = 'Learjet'
        i0 = ictFile.find('Learjet')+8
    m = ICARTT(ictFile)
    alt, lon, lat, tyme = m.Nav['Altitude'], m.Nav['Longitude'], m.Nav['Latitude'], m.Nav['Time']
    yyyymmdd = ictFile[i0:-4]
    dateout  = ictFile[i0:i0+4]+"-"+ictFile[i0+4:i0+6]+"-"+ictFile[i0+6:i0+8]
    print(ictFile, aircraft, yyyymmdd)

    fig, ax = plt.subplots(figsize=(16, 6))
    plt.subplots_adjust(left=0.07,bottom=0.1,right=0.99,top=0.9)
    if aircraft == "ER2":
        ax.set_ylim(0,22)
    else:
        ax.set_ylim(0,16)
    plt.ylabel('GPS Altitude [km]')
    dtFmt = mdates.DateFormatter('%H:%M') # define the formatting
    plt.gca().xaxis.set_major_formatter(dtFmt) # apply the format to the desired axis
    cf,  = ax.plot(tyme,alt/1000.,color='red',linewidth=4,label='Altitude',zorder=100)
    for item in ([ax.title, ax.xaxis.label, ax.yaxis.label] +
             ax.get_xticklabels() + ax.get_yticklabels()):
        item.set_fontsize(20)
    return fig, ax, aircraft, yyyymmdd, cf

