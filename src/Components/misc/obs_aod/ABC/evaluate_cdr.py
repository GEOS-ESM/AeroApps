"""
Reduces AVHRR according to avering kernel and speciation.
"""

import os
import sys

from numpy             import arange, log, ones
from matplotlib.pyplot import savefig

from datetime  import date, datetime, timedelta

from sknet     import loadnet

from abc_avhrr import ABC_AVHRR

from pyobs     import kde, NPZ

from grads     import GrADS

Bins = dict ( all = (-3, 0),
              fdu = (-3, 0),
              fss = (-3, 0), # (-3, -1.2),
              fcc = (-3, 0),
              fsu = (-3, 0), # (-3, -1.2),
            )

Species = dict ( all = 'All',
                 fdu = 'Dust',
                 fss = 'Sea-salt',
                 fcc = 'Carbonaceous',
                 fsu = 'Sulfate',
            )

MISSING = -1e20

#--------------------------------------------------------------------------------
def aodInterpAngs(lambda_,tau1,tau2,lambda1,lambda2):
    """
       Log-interpolated AOD.
    """
    I = (tau1>0) & (tau2>0)
    angstrom = -log(tau1[I]/tau2[I])/log(lambda1/lambda2)
    tau = MISSING * ones(len(tau1))
    tau[I] = tau2[I] * (lambda2/lambda_)**angstrom
    return tau

#--------------------------------------------------------------------------------
if __name__ == "__main__":

    ident = 'avhrr'
    expid = 'cdr'

    RootDir = '/nobackup/AVHRR/Level2/NPZ'

    ext_Nc = '/nobackup/MERRAero/ext_Nc.ddf'

    ga = GrADS(Window=False,Echo=False)
    fh = ga.open(ext_Nc)
    
    #for year in (2005,2006,2007,2008,2009):
    for year in (2008,):

        a = ABC_AVHRR(Path='/nobackup/AVHRR/Level2/NPZ/%d/*.npz'%year,
                      N_bal=400*1000,Verb=False)

        ga('set lev 675e-9')
        # a.tau_675 = ga.sampleXYT('taod',a.lon,a.lat,a.tyme,Verbose=True)
        n = NPZ('mTau_675.npz')
        a.tau_675 = n.tau_675

        mtau_630 = aodInterpAngs(630.,a.tau_550,a.tau_675,550.,675.)

        a.iValid = a.iValid & (mtau_630>0)
        
        ##for s in Species.keys():
        for s in ('fss', 'fsu'):

            figFile = '%s.%s.%d.%s.png'%(ident,expid,year,s)
            statFile = '%s.%s.%d.%s.txt'%(ident,expid,year,s)

            if os.path.exists(figFile):
                print "Skipping <%s> ..."%Species[s]
                continue
            else:
                print "Working on <%s> for %s ..."%(Species[s],expid)

            if s == 'all':
                I = a.iValid 
            else:
                I = a.iValid & (a.__dict__[s]>0.5)

            targets = log( mtau_630[I]+0.01)  # from MERRAero
            results = log(a.tau_630[I]+0.01)  # from CDR
     
            x_bins, y_bins, P = kde.calc_kde2d(targets,results,
                                               x_range=Bins[s],y_range=Bins[s],
                                               Nx=128, Ny=128)

            kde.plot_kde2d(x_bins, y_bins, P,
                       Title  = 'AVHHR %s Aerorol Optical Depth (%d)'\
                                %(Species[s],year),
                       xLabel = r'MERRAero $\tau_{630}$',
                       yLabel = r'CDR $\tau_{630}$',
                       formatter=kde.aodFormat(),
                       )
            
            savefig(figFile)
