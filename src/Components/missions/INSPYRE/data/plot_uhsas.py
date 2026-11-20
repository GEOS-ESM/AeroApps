#!/usr/bin/env python3
import sys
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
import os
from datetime import timedelta
from optparse import OptionParser   # Command-line args
from plothelp import plot_timeseries

uhsas_nbin    = 94

uhsas_diam_nm = [70.15, 72.17, 74.25, 76.39, 78.6, 80.86, 83.19, 85.59, 88.06,
                 90.6, 93.21, 95.89, 98.66, 101.5, 104.4, 107.4, 110.5, 113.7,
                 117, 120.4, 123.8, 127.4, 131.1, 134.9, 138.8, 142.8, 146.9,
                 151.1, 155.5, 159.9, 164.5, 169.3, 174.2, 179.2, 184.4, 189.7,
                 195.1, 200.8, 206.5, 212.5, 218.6, 224.9, 231.4, 238.1, 244.9,
                 252, 259.3, 266.7, 274.4, 282.3, 290.5, 298.9, 307.5, 316.3,
                 325.5, 334.8, 344.5, 354.4, 364.6, 375.1, 386, 397.1, 408.5,
                 420.3, 432.4, 444.9, 457.7, 470.9, 484.5, 498.5, 512.8, 527.6,
                 542.8, 558.5, 574.6, 591.1, 608.2, 625.7, 643.7, 662.3, 681.4,
                 701, 721.2, 742, 763.4, 785.4, 808, 831.3, 855.3, 880, 905.3,
                 931.4, 958.3, 985.9]

uhsas_dlogD   = 0.01234

#Get the HALO RGB colors
import csv
from matplotlib.colors import LinearSegmentedColormap
with open('/home/pcolarco/lib/halo_colorbar.csv', newline='') as csvfile:
    reader = csv.DictReader(csvfile, delimiter=',',fieldnames=['r','g','b'])
    i = 0
    rgb = []
    for row in reader:
        rgb.append([int(row['r'])/255.,int(row['g'])/255.,int(row['b'])/255.])
cm = LinearSegmentedColormap.from_list(
        'my_map', rgb)

def uhsas_psd(x):
#   Assemble UHSAS PSD into a simple 2d array
    nt  = len(x.tyme)
    psd = np.zeros((uhsas_nbin,nt))
    for i in range(0,uhsas_nbin):
        varn = f"UHSAS_Bin{i+1:02}"
        varv = x.__dict__[varn][:]
        psd[i,:] = varv
    return psd

def plotext(ictFile):

    fig, ax, aircraft, yyyymmdd, cf  = plot_timeseries(ictFile)

#   See if there is a CAPS file to read
    if aircraft == "GV":
        fn = f"../data/INSPYRE-PUTLS-UHSAS_GV_{yyyymmdd}.ict"        
        if os.path.exists(fn):
            ax3 = ax.twinx()
            ax3.spines.right.set_position(("axes",1))
            ax3.set(ylabel="Diameter [nm]")
            ax3.set_yscale("log")
            for item in ([ax3.title, ax3.xaxis.label, ax3.yaxis.label] +
                ax3.get_xticklabels() + ax3.get_yticklabels()):
                item.set_fontsize(16)
            x = ICARTT(fn)
            psd = uhsas_psd(x)
            clevs = np.arange(2,4.05,.05)
            cv = ax3.contourf(x.tyme,uhsas_diam_nm,np.log10(psd),clevs,
                              cmap=cm,zorder=50,extend="both")
            cbar = plt.colorbar(cv,ticks=[2,3,4],pad=0.1,extend="max",
                    format=mticker.FixedFormatter(['100','1000','10000']))
            cbar.ax.tick_params(labelsize=12)
            cbar.set_label(label="# cm-3",
                    size=16,rotation=270.,labelpad=25)
        else:
            ax.legend(handles=[cf])
    else:
        ax.legend(handles=[cf])

    yvals = [70,100,150,200,300,400,500,600,700,800,1000]
    ax3.set_yticks(yvals)
    ax3.set_yticklabels(yvals)
    dateout  = yyyymmdd[0:4]+"-"+yyyymmdd[4:6]+"-"+yyyymmdd[6:8]
    plt.title('UHSAS PSD, %s track: '%(aircraft)+dateout, size=20)
    ofname = f"./INSPYRE-{aircraft}_UHSAS_{yyyymmdd}.png"
    plt.savefig(ofname)
    plt.close(fig)

if __name__ == "__main__":

    parser = OptionParser(usage="Usage: %prog [options] modelname date0",
                          version='xxx' )
    (options, args) = parser.parse_args()
 
#  GET ICT FILE FROM INPUT ARGUMENT LIST
    if len(args) == 1:
        ict        = args[0]
    else:
        parser.error("must have 0 argument: icartt")

    plotext(ict)
