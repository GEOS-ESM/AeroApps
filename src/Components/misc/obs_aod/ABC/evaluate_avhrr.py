"""
Reduces AVHRR according to avering kernel and speciation.
"""

import os
import sys

from numpy             import arange
from matplotlib.pyplot import savefig

from datetime  import date, datetime, timedelta

from sknet     import loadnet

from abc_avhrr import ABC_AVHRR

from pyobs     import kde


Bins = dict ( all = (-3, 0),
              fdu = (-3, 0),
              fss = (-3, -1.2),
              fcc = (-3, 0),
              fsu = (-3, -1.2),
            )

Species = dict ( all = 'All',
                 fdu = 'Dust',
                 fss = 'Sea-salt',
                 fcc = 'Carbonaceous',
                 fsu = 'Sulfate',
            )


if __name__ == "__main__":

    ident = 'avhrr'
    expid = 'nnr_001e'

    RootDir = '/nobackup/AVHRR/Level2/NPZ'

    net = loadnet('Net/%s.avhrr_Tau.net'%expid)

    for year in (2005,2006,2007,2008,2009):

        a = ABC_AVHRR(Path='/nobackup/AVHRR/Level2/NPZ/%d/*.npz'%year,
                      N_bal=400*1000,Verb=False)

        a.net = net
        a.Input = net.InputNames
        a.Target = net.TargetNames

        for s in list(Species.keys()):

            figFile = '%s.%s.%d.%s.png'%(ident,expid,year,s)
            statFile = '%s.%s.%d.%s.txt'%(ident,expid,year,s)

            if os.path.exists(figFile):
                print("Skipping <%s> ..."%Species[s])
                continue
            else:
                print("Working on <%s> for %s ..."%(Species[s],expid))

            if s == 'all':
                I = a.iValid 
            else:
                I = a.iValid & (a.__dict__[s]>0.5)

            out, reg = a.test(fname=statFile)
            results = a.eval(I)
            targets = a.getTargets(I)

            x_bins, y_bins, P = kde.calc_kde2d(targets,results,
                                           x_range=Bins[s],y_range=Bins[s],
                                           Nx=128, Ny=128)

            kde.plot_kde2d(x_bins, y_bins, P,
                       Title  = 'AVHHR %s Aerorol Optical Depth (%d)'\
                                %(Species[s],year),
                       xLabel = r'MODIS $\tau_{550}$',
                       yLabel = r'NNR $\tau_{550}$',
                       formatter=kde.aodFormat(),
                       )
            
            savefig(figFile)
