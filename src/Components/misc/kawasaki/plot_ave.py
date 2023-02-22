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
    month = list(range(12))

    # 16:9 Figure
    # -----------
    figure(dpi=120,figsize=(11,11*9./16))

    # Surface Mass
    for q in ('dusmass25','so4smass','bcsmass','ocsmass',
              'sssmass25','so2smass'):
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
    axis([-0.5,11.5,0,13])
    legend(loc='upper center',prop={'size':10})
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

def plot_clm2(trange,clm,showKD=True):

    if showKD:
        c = NPZ('NPZ/kawasaki_japan.climatology.%s.npz'%trange)

    month = list(range(12))

    # 16:9 Figure
    # -----------
    figure(dpi=120,figsize=(11,11*9./16))

    so4 = clm['so4smass']
    sss = clm['sssmass25']
    dus = clm['dusmass25']
    ocs = clm['ocsmass']
    bcs = clm['bcsmass']
    
    f = 1e9
    plot(month,f*so4,'m-',linewidth=2,label=r'SO$_4$')
    plot(month,f*sss,'b-',linewidth=2,label='Sea Salt')
    plot(month,f*dus,'-',color='orange',linewidth=2,label='Dust')
    plot(month,f*ocs,'-',color='brown',linewidth=2,label='OC')
    plot(month,f*bcs,'-',color='grey',linewidth=2,label='BC')

    ylabel(r'$\mu$g/m$^3$')
    title( trange+' Surface PM2.5 Concentration over Japan')
    grid()
    axis([-0.5,11.5,0,7])
    legend(loc='upper right',prop={'size':10})

    if showKD:
        ax2 = twinx()
        p = 100*c.all.astype('float')/c.all.sum()
        ax2.plot(month,p,'k-o',linewidth=3,label='% KD Cases')
        ylabel('% KD Cases per Month')
        axis([-0.5,11.5,6,12])
        legend(loc='upper left',prop={'size':10})

    xticks(month,['Jan','Feb','Mar','Apr','May','Jun','Jul',
                  'Aug','Sep','Oct','Nov','Dec'])

    savefig('clm2.sfc.'+trange+'.png',dpi=120)

    # AOD
    # ---
    clf()
    so4 = clm['suexttau']
    sss = clm['ssexttau']
    dus = clm['duexttau']
    ocs = clm['ocexttau']
    bcs = clm['bcexttau']

    plot(month,so4,'m-',linewidth=2,label=r'SO$_4$')
    plot(month,sss,'b-',linewidth=2,label='Sea Salt')
    plot(month,dus,'-',linewidth=2,color='orange',label='Dust')
    plot(month,ocs,'-',linewidth=2,color='brown',label='OC')
    plot(month,bcs,'-',linewidth=2,color='grey',label='BC')

    title( trange+' Aerosol Optical Depth over Japan')
    grid()
    #axis([0,11,0,0.13])
    legend(loc='upper right',prop={'size':10})

    if showKD:
        ax2 = twinx()
        p = 100*c.all.astype('float')/c.all.sum()
        ax2.plot(month,p,'k-o',linewidth=3,label='% KD Cases')
        ylabel('% KD Cases per Month')
        axis([0,11,6,12])
        legend(loc='upper left',prop={'size':10})

    xticks(month,['Jan','Feb','Mar','Apr','May','Jun','Jul',
                  'Aug','Sep','Oct','Nov','Dec'])

    savefig('clm2.aod.'+trange+'.png',dpi=120)

def plot_clm_stack(trange,clm,showKD=True):

    if showKD:
        c = NPZ('NPZ/kawasaki_japan.climatology.%s.npz'%trange)

    month = list(range(12))

    # 16:9 Figure
    # -----------
    figure(dpi=120,figsize=(11,11*9./16))

    so4 = clm['dusmass25']+clm['so4smass']+clm['bcsmass']\
        + clm['ocsmass']+clm['sssmass25']
    dus = clm['dusmass25']+clm['bcsmass']\
        + clm['ocsmass']+clm['sssmass25']
    sss = clm['bcsmass']\
        + clm['ocsmass']+clm['sssmass25']
    ocs = clm['bcsmass']+clm['ocsmass']
    bcs = clm['bcsmass']
    
    f = 1e9
    fill_between(month,0,f*so4,color='magenta')
    fill_between(month,0,f*dus,color='orange')
    fill_between(month,0,f*sss,color='blue')
    fill_between(month,0,f*ocs,color='brown')
    fill_between(month,0,f*bcs,color='grey')

    plot(month,f*so4,'m-',linewidth=2,label=r'SO$_4$')
    plot(month,f*dus,'-',linewidth=2,color='orange',label='Dust')
    plot(month,f*sss,'b-',linewidth=2,label='Sea Salt')
    plot(month,f*ocs,'-',linewidth=2,color='brown',label='OC')
    plot(month,f*bcs,'-',linewidth=2,color='grey',label='BC')

    ylabel(r'$\mu$g/m$^3$')
    title( trange+' Surface PM2.5 Concentration over Japan')
    grid()
    axis([0,11,0,18])
    legend(loc='upper right',prop={'size':10})

    if showKD:
        ax2 = twinx()
        p = 100*c.all.astype('float')/c.all.sum()
        ax2.plot(month,p,'k-o',linewidth=2,label='% KD Cases')
        ylabel('% KD Cases per Month')
        axis([0,11,6,12])
        legend(loc='upper left',prop={'size':10})

    xticks(month,['Jan','Feb','Mar','Apr','May','Jun','Jul',
                  'Aug','Sep','Oct','Nov','Dec'])

    savefig('clm.sfc_s.'+trange+'.png',dpi=120)

    # AOD
    # ---
    clf()

    so4 = clm['duexttau']+clm['suexttau']+clm['bcexttau']\
        + clm['ocexttau']+clm['ssexttau']
    dus = clm['duexttau']+clm['bcexttau']\
        + clm['ocexttau']+clm['ssexttau']
    sss = clm['bcexttau']\
        + clm['ocexttau']+clm['ssexttau']
    ocs = clm['bcexttau']+clm['ocexttau']
    bcs = clm['bcexttau']

    
    fill_between(month,0,so4,color='magenta')
    fill_between(month,0,dus,color='orange')
    fill_between(month,0,sss,color='blue')
    fill_between(month,0,ocs,color='brown')
    fill_between(month,0,bcs,color='grey')

    plot(month,so4,'m-',linewidth=2,label=r'SO$_4$')
    plot(month,dus,'-',linewidth=2,color='orange',label='Dust')
    plot(month,sss,'b-',linewidth=2,label='Sea Salt')
    plot(month,ocs,'-',linewidth=2,color='brown',label='OC')
    plot(month,bcs,'-',linewidth=2,color='grey',label='BC')

    title( trange+' Aerosol Optical Depth over Japan')
    grid()
    axis([0,11,0,0.42])
    legend(loc='upper right',prop={'size':10})

    if showKD:
        ax2 = twinx()
        p = 100*c.all.astype('float')/c.all.sum()
        ax2.plot(month,p,'k-o',linewidth=2,label='% KD Cases')
        ylabel('% KD Cases per Month')
        axis([0,11,6,12])
        legend(loc='upper left',prop={'size':10})

    xticks(month,['Jan','Feb','Mar','Apr','May','Jun','Jul',
                  'Aug','Sep','Oct','Nov','Dec'])

    savefig('clm.aod_s.'+trange+'.png',dpi=120)

#....................................................................

if __name__ == "__main__":

    # Instantiate GrADS
    # -----------------
    
    y1, y2 = (2000,2015)
    showKD = False
    # y1, y2 = (1980,1999)
    #y1, y2 = (2011,2013)

    d = dict( y1=y1, y2=y2,
              fvIn='/discover/nobackup/projects/gmao/share/dasilva/fvInput/g5chem/sfc' )
              

    clm = dict()
    ano = dict()

    for q in ('ducmass25','so4cmass','bccmass','occmass','so2cmass','sscmass25',
              'dusmass25','so4smass','bcsmass','ocsmass','so2smass','sssmass25',
              'duexttau','suexttau','bcexttau','ocexttau','ssexttau' ):
 
        print('[] working on <%s>'%q)
        
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
    plot_clm2('%d-%d'%(y1,y2),clm,showKD)
    plot_clm_stack('%d-%d'%(y1,y2),clm,showKD)
