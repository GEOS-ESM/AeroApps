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

def plotext(ictFile):

    fig, ax, aircraft, yyyymmdd, cf  = plot_timeseries(ictFile)

    ax3 = ax.twinx()
    ax3.spines.right.set_position(("axes", 1.075))
    ax3.set(ylabel="Black Carbon [ng m-3]")
    ax3.set_yscale("log")
    for item in ([ax3.title, ax3.xaxis.label, ax3.yaxis.label] +
             ax3.get_xticklabels() + ax3.get_yticklabels()):
        item.set_fontsize(16)



#   See if there is SP2 file to read
    if aircraft == "GV":
        fn = f"../data/INSPYRE-SP2-BC-1Hz_GV_{yyyymmdd}.ict"
        if os.path.exists(fn):
            x = ICARTT(fn)
#           Hack for time offset in first files
            p4, = ax3.plot(x.tyme+timedelta(hours=6),x.BC_mass_90_550_nm_amb,color="dimgray",label="SP2 BC Mass",zorder=50)
#            p4, = ax3.plot(x.tyme,x.BC_mass_90_550_nm_amb,color="black",label="SP2 BC Mass",zorder=50)
            ax.legend(handles=[cf,p4])
        else:
                ax.legend(handles=[cf])
    else:
            ax.legend(handles=[cf])

    dateout  = yyyymmdd[0:4]+"-"+yyyymmdd[4:6]+"-"+yyyymmdd[6:8]
    plt.title('Aerosol Concentration, %s track: '%(aircraft)+dateout, size=20)
    ofname = f"./INSPYRE-{aircraft}_SP2_{yyyymmdd}.png"
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
