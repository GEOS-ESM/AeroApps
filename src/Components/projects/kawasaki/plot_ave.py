import os
import sys

from numpy    import linspace, array, savez
from glob     import glob
from datetime import datetime, timedelta
from string   import Template

from matplotlib.pyplot import *

from grads import GrADS

from npz import NPZ

class myGrADS(GrADS):

    def Cmd(self,tmpl,d):
        """
        Wrapper for cmd() 
        """
        #print Template(tmpl).substitute(d)
        self.__call__(Template(tmpl).substitute(d))

def plot_clm(trange,clm):

    c = NPZ('NPZ/kawasaki_japan.climatology.2000-2010.npz')
    month = range(12)
    figure(dpi=120)

    # Surface Mass
    for q in ('dusmass','so4smass','bcsmass','ocsmass','sssmass','so2smass'):
        if q=='sssmass':
            f = 0.5e9
        else:
            f = 1e9
        plot(month,f*clm[q],'-o',linewidth=1.5,
             label=q.replace('smass','').upper())
    xticks(month,['Jan','Feb','Mar','Apr','May','Jun','Jul',
                  'Aug','Sep','Oct','Nov','Dec'])
    ylabel(r'$\mu$g/m$^3$')
    xlabel('Climatological Month')
    title( trange+' Surface Concentration over Japan')
    #legend(loc='upper left')
    #legend(loc='upper right')
    grid()
    plot(month,100.*c.all/c.all.sum(),'k-o',linewidth=3,label='% KD Cases')
    axis([-0.5,11.5,0,18])
    legend(loc='upper left',prop={'size':10})
    #ax2 = twinx()
    
    savefig('clm.sfc.'+trange+'.png',dpi=180)

    return

    # AOD
    clf()
    for q in ('duexttau', 'suexttau', 'ocexttau','bcexttau','ssexttau'):
        plot(month,clm[q],'-o',linewidth=1.5,label=q[0:2].upper())
    xticks(month,['Jan','Feb','Mar','Apr','May','Jun','Jul',
                  'Aug','Sep','Oct','Nov','Dec'])
    title( trange+' Climatological AOD over Japan')
    legend()
    grid()
    axis([-0.5,11.5,0,0.25])
    savefig('clm.aod.'+trange+'.png',dpi=180)

#....................................................................

if __name__ == "__main__":

    # Instantiate GrADS
    # -----------------
    
    y1, y2 = (2000,2015)

    d = dict( y1=y1, y2=y2,
              fvIn='/discover/nobackup/projects/gmao/share/dasilva/fvInput/g5chem/sfc' )
              

    clm = dict()
    ano = dict()

    for q in ('ducmass','so4cmass','bccmass','occmass','so2cmass','sscmass',
              'dusmass','so4smass','bcsmass','ocsmass','so2smass','sssmass',
              'duexttau','suexttau','bcexttau','ocexttau','ssexttau' ):
 
        print '[] working on <%s>'%q
        
        d['q'] = q

        # Climatology
        # -----------
        ga = myGrADS(Window=False,Echo=False)
        ga.Cmd("""
           sdfopen $fvIn/countries.x576_y361.nc4
           sdfopen Anom/MERRA2_clm.$q.$y1-$y2.nc4
           set dfile 2 

           set x 1
           set y 1
           set t 1 12

           define $q = tloop(aave(${q}clm*if(countries.1(t=1),==,200,1,-u),lon=128,lon=148,lat=30,lat=46))

          """,d)

        clm[q] = ga.expr(q)

    # Plot climatology
    plot_clm('%d-%d'%(y1,y2),clm)
