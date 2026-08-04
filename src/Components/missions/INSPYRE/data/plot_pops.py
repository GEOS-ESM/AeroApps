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

pops_nbin    = 73

pops_diam_nm = [204.7, 209.1, 214, 219.3, 225.1, 231.4, 238.4, 245.9, 254.2,
                263.3, 273.2, 284.1, 295.9, 308.8, 323, 338.5, 355.4, 373.9,
                394.1, 416.2, 440.4, 466.8, 495.7, 527.2, 559.1, 583.5,
                601.7, 620.4, 639.7, 659.6, 680.2, 701.3, 723.2, 745.7, 768.9,
                792.8, 817.5, 842.9, 869.1, 896.2, 924.1, 952.8, 982.5, 1013,
                1045, 1077, 1111, 1145, 1181, 1218, 1255, 1295, 1335,
                1376, 1419, 1463, 1509, 1556, 1604, 1654, 1706, 1759, 1814,
                1870, 1928, 1988, 2050, 2114, 2180, 2248, 2317, 2390, 2464]

pops_dlogD   = [0.009014, 0.009645, 0.01031, 0.01099, 0.01171, 0.01245, 0.01321,
                0.014, 0.0148, 0.01563, 0.01646, 0.01731, 0.01817, 0.01903,
                0.01989, 0.02075, 0.02161, 0.02245, 0.02329, 0.02411, 0.02491,
                0.02569, 0.02645, 0.02718, 0.02382, 0.01331, 0.01331, 0.01331,
                0.01331, 0.01331, 0.01331, 0.01331, 0.01331, 0.01331, 0.01331,
                0.01331, 0.01331, 0.01331, 0.01331, 0.01331, 0.01331, 0.01331,
                0.01331, 0.01331, 0.01331, 0.01331, 0.01331, 0.01331, 0.01331,
                0.01331, 0.01331, 0.01331, 0.01331, 0.01331, 0.01331, 0.01331,
                0.01331, 0.01331, 0.01331, 0.01331, 0.01331, 0.01331, 0.01331,
                0.01331, 0.01331, 0.01331, 0.01331, 0.01331, 0.01331, 0.01331,
                0.01331, 0.01331, 0.01331]

def pops_psd(x):
#   Assemble POPS PSD into a simple 2d array
    nt  = len(x.tyme)
    psd = np.zeros((pops_nbin,nt))
    for i in range(0,pops_nbin):
        varn = f"POPS_Bin{i+1:02}"
        varv = x.__dict__[varn][:]
        psd[i,:] = varv
    return psd

def plotext(ictFile):

    fig, ax, aircraft, yyyymmdd, cf  = plot_timeseries(ictFile)

    ax3 = ax.twinx()
    ax3.spines.right.set_position(("axes",1))
    ax3.set(ylabel="# cm-3 dNdlogD")
    ax3.set_yscale("log")
    for item in ([ax3.title, ax3.xaxis.label, ax3.yaxis.label] +
             ax3.get_xticklabels() + ax3.get_yticklabels()):
        item.set_fontsize(16)



#   See if there is a CAPS file to read
    if aircraft == "GV":
        fn = f"../data/INSPYRE-PUTLS-POPS_GV_{yyyymmdd}.ict"        
        if os.path.exists(fn):
            x = ICARTT(fn)
            psd = pops_psd(x)
            clevs = [1,2,5,10,20,50,100,200,500,1000]
            cv = ax3.contourf(x.tyme,pops_diam_nm,psd,clevs,zorder=50)
        else:
            ax.legend(handles=[cf])
    else:
        ax.legend(handles=[cf])

    dateout  = yyyymmdd[0:4]+"-"+yyyymmdd[4:6]+"-"+yyyymmdd[6:8]
    plt.title('POPS PSD, %s track: '%(aircraft)+dateout, size=20)
    ofname = f"./INSPYRE-{aircraft}_POPS_{yyyymmdd}.png"
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
