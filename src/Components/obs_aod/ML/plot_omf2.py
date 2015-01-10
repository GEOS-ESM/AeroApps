#!/usr/bin/env python

from grads import *
from pylab import *

if __name__ == "__main__":

    channel = 470
    yy = 2008
    mm = 06
    yymm = 2008*100 + mm

    
    qname = {}
    qname['aod_modo'] = 'MODIS/TERRA Ocean'
    qname['aod_modl'] = 'MODIS/TERRA Land'
    qname['aod_mydo'] = 'MODIS/AQUA  Ocean'
    qname['aod_mydl'] = 'MODIS/AQUA  Land'
    qname['aod_deep'] = 'MODIS/AQUA  Deep-Blue'
    qname['aod_misr'] = 'MISR'
    qname['aod_anet'] = 'AERONET'
    qname['aod_parl'] = 'PARASOL Land'
    qname['aod_paro'] = 'PARASOL Ocean'
    qname['aod_omi']  = 'OMI'

    ga = GrADS(Window=False,Echo=False)

    ga("sdfopen a0005.gritas_aod.obs.%s.nc"%yymm)
    ga("sdfopen a0005.gritas_aod.omf.%s.nc"%yymm)
    ga("sdfopen a0005.gritas_aod.oma.%s.nc"%yymm)

    ga('set lev %d'%channel)
    ga('set gxout grfill')

    for q in ( 'aod_modo', 'aod_modl', 'aod_mydo', 'aod_mydl',
               'aod_misr', 'aod_deep', 'aod_omi' ):

        print "Plotting "+q

        ga('clear')
        ga("d %s.1"%q)
        ga('cbarn')
        ga('draw title %dnm OBS Log(eps+AOD) - %s [%d-%d]'%(channel,qname[q],yy,mm))
        ga('gxyat obs.%s_%dnm.%d.png'%(q,channel,yymm))
        
        ga('clear')
        ga("xydiff %s.2"%q)
        ga('draw title %dnm O-F Log(eps+AOD) - %s [%d-%d]'%(channel,qname[q],yy,mm))
        ga('gxyat omf.%s_%dnm.%d.png'%(q,channel,yymm))
        ga('print omf.%s_%dnm.%d.eps'%(q,channel,yymm))
        
        ga('clear')
        ga("xydiff %s.3"%q)
        ga('draw title %dnm O-A Log(eps+AOD) - %s [%d-%d]'%(channel,qname[q],yy,mm))
        ga('gxyat oma.%s_%dnm.%d.png'%(q,channel,yymm))



     
    
