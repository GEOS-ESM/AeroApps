#!/usr/bin/env python

"""

   Simple scripts to collect simple innovation stats for SQC.

"""

from pylab import *
from pyods import *

import scipy.stats as stats

if __name__ == "__main__":

    basen = 'logtau'
    ch = 550
    log_transf = True

    ods = ODS(basen+".ods",20080629,120000,only_good=True)

    ods = ods.select(lev=ch)
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

    clf()

    print " kx   mean   stdv  |    0    25%    50%    75%    100%"
    print "---- ------ ------ | ------ ------ ------ ------ ------ "
    x = linspace(-5,2,100)
    i = 0
    for kx in ( 301, 302, 313,  311, 312, 320 ):

        i = i + 1

        ods_ = ods.select(kx=kx)
        
        d_f = ods_.omf

        m = d_f.mean()
        s = d_f.std()

        if d_f.size > 0.:

            g = stats.norm(loc=m,scale=s)

            subplot(2,3,i)
            if log_transf:
                hist(d_f,100,range=(-5,2),normed=True) 
            else:
                hist(d_f,100,range=(-1,1),normed=True) 

            plot(x,g.pdf(x),'r',linewidth=2)
            title("%s [%dnm]"%(kx_names[str(kx)],ch),fontsize=12)
            args = (kx,d_f.mean(),d_f.std())+tuple(prctile(d_f))

            print " %d %6.3f %6.3f | %6.3f %6.3f %6.3f %6.3f %6.3f "%args

        savefig(basen+"-%dnm.png"%ch)
