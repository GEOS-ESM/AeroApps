#!/usr/bin/env python

"""

   Simple scripts to collect simple innovation stats for SQC.

"""

from pylab import *
from pyods import *

import scipy.stats as stats

if __name__ == "__main__":

#    basen = 'misr.200806'
    basen = 'modo550.200806'
    ch = 550
    log_transf = True
    eps = 0.01
    
    kx_names = {}
    kx_names['301'] = "MOD04 (Ocean)"
    kx_names['302'] = "MOD04 (Land)"
    kx_names['311'] = "MYD04 (Ocean)"
    kx_names['312'] = "MYD04 (Land)"
    kx_names['313'] = "MISR"
    kx_names['314'] = "OMI"
    kx_names['316'] = "PARASOL (Land)"
    kx_names['317'] = "PARASOL (Ocean)"
    kx_names['320'] = "MYD04 (Deep Blue)"
    kx_names['324'] = "AVHRR NNR (PATMOSX)"

#   Concatenate obs/o-f for given period of time
#   --------------------------------------------
    OBS = ()
    OMF = ()
    for nymd in range(20080601,20080607):
        for nhms in range(0,240000,30000):
            ods = ODS(basen+".ods",nymd,nhms,only_good=True)
            try:
                ods = ods.select(lev=ch)
                nobs = ods.nobs
                OBS = OBS + (ods.obs,)
                OMF = OMF + (ods.omf,)
            except:
                nobs = 0
            print nymd, nhms, nobs

#   Concatenate and undo log-transform
#   ----------------------------------
    Lo = concatenate(OBS)
    Ld = concatenate(OMF)
    Lf = Lo - Ld
    o = exp(Lo) - eps
    f = exp(Lf) - eps
    d = o - f

#   Stats
#   -----
    Go = stats.norm(loc=o.mean(),scale=o.std())
    Gd = stats.norm(loc=d.mean(),scale=d.std())
    GLo = stats.norm(loc=Lo.mean(),scale=Lo.std())
    GLd = stats.norm(loc=Ld.mean(),scale=Ld.std())

#   Plot
#   ----
    clf()

#   AOD
#   ---
    subplot(221)
    x = linspace(-0.5, 1., 100) 
    hist(o,100,range=(-0.5,1),normed=True)
    plot(x,Go.pdf(x),'r',linewidth=2)
    title(r'$\tau^o_{%d}$'%ch,fontsize=12)

#   LAOD
#   -----
    subplot(222)
    x = linspace(-5., 1., 100)
    hist(Lo,100,range=(-5.,1),normed=True)
    plot(x,GLo.pdf(x),'r',linewidth=2)
    title(r'Log($(0.01 + \tau^o_{%d}$)'%ch,fontsize=12)

#   AOD Inc
#   -------
    subplot(223)
    x = linspace(-0.5, 0.5, 100)
    hist(d,100,range=(-0.5,0.5),normed=True)
    plot(x,Gd.pdf(x),'r',linewidth=2)
    title(r'$\Delta \tau^o_{%d}$'%ch,fontsize=12)

#   LAOD Inc
#   --------
    subplot(224)
    x = linspace(-2., 2., 100) 
    hist(Ld,100,range=(-2.,2),normed=True)
    plot(x,GLd.pdf(x),'r',linewidth=2)
    title(r'$\Delta$Log($(0.01 + \tau^o_{%d}$)'%ch,fontsize=12)




