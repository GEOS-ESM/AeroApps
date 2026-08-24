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
    ax3.spines.right.set_position(("axes",1))
    ax3.set(ylabel="Extinction/Scattering [Mm-1]")
    ax3.set_yscale("log")
    for item in ([ax3.title, ax3.xaxis.label, ax3.yaxis.label] +
             ax3.get_xticklabels() + ax3.get_yticklabels()):
        item.set_fontsize(16)



#   See if there is a CAPS file to read
    if aircraft == "GV":
        fn = f"../data/INSPYRE-PUTLS-CAPS_GV_{yyyymmdd}.ict"
        if os.path.exists(fn):
            x = ICARTT(fn)
#           Hack for time offset in first files
            p4, = ax3.plot(x.tyme,x.Ext_530nm_submicron/x.stdPT,color="dimgray",label="CAPS Extinction",zorder=50)
            p5, = ax3.plot(x.tyme,x.Scat_530nm_submicron/x.stdPT,color="cornflowerblue",label="CAPS Scattering",zorder=50)
            ax.legend(handles=[cf,p4,p5])
        else:
                ax.legend(handles=[cf])
    else:
            ax.legend(handles=[cf])

    dateout  = yyyymmdd[0:4]+"-"+yyyymmdd[4:6]+"-"+yyyymmdd[6:8]
    plt.title('Aerosol Scattering & Extinction, %s track: '%(aircraft)+dateout, size=20)
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
