import os
import sys

from numpy    import load, linspace, array, savez
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

#...................................................

def _getTime(y1,y2):
    Tyme = []
    for year in range(y1,y2+1):
        for month in range(1,13):
            Tyme.append(year + float(month-1)/12.)
    return array(Tyme)
     
def _getTyme(y1,y2):
    Tyme = []
    for year in range(y1,y2+1):
        for month in range(1,13):
            Tyme.append(datetime(year,month,15))
    return array(Tyme)
     
def _toYearly(Year,q):
    Q = []
    for year in range(Year[0],Year[-1]+1):
        I = (Year==year)
        Q.append(q[I].mean())
    return array(Q)

def _getYearly(y1,y2):
    n = load('NPZ/m2_japan.monthly.%d-%d.npz'%(y1,y2))
    tyme = _getTyme(y1,y2)
    Year = array([t.year for t in tyme])
    clm = dict()
    for q in n:
        print("[] working on <%s>"%q)
        clm[q] = _toYearly(Year,n[q])
    return clm

def plot_Japan(trange,clm,showKD=True):

    if showKD:
        c = NPZ('NPZ/kawasaki_japan.climatology.%s.npz'%trange)        

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

    ylabel(r'$\mu$g/m$^3$')
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

    savefig('clm2.aod.'+trange+'.png',dpi=120,bbox_inches='tight')

def plot_Japan_stack(trange,clm,showKD=True,Yearly=False):

    y1, y2 = trange.split('-')
    y1, y2 = int(y1), int(y2)
    if Yearly:
        tyme = list(range(y1,y2+1))
        tag = 'yy'
    else:
        tag = 'mm'
        tyme = _getTime(y1,y2) # monthly

    if showKD:
        if Yearly:
            c = NPZ('NPZ/kawasaki_japan.yearly.01-12.npz')
            c.tyme = array(list(range(1969,2011)))
            I = (c.tyme>=y1)
        else:
            c = NPZ('NPZ/kawasaki_japan.monthly.1969-2010.npz')
            c.tyme = _getTime(1969,2010)
            I = (c.tyme>=tyme[0])

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
    fill_between(tyme,0,f*so4,color='magenta')
    fill_between(tyme,0,f*dus,color='orange')
    fill_between(tyme,0,f*sss,color='blue')
    fill_between(tyme,0,f*ocs,color='brown')
    fill_between(tyme,0,f*bcs,color='grey')

    plot(tyme,f*so4,'m-',linewidth=2,label=r'SO$_4$')
    plot(tyme,f*dus,'-',linewidth=2,color='orange',label='Dust')
    plot(tyme,f*sss,'b-',linewidth=2,label='Sea Salt')
    plot(tyme,f*ocs,'-',linewidth=2,color='brown',label='OC')
    plot(tyme,f*bcs,'-',linewidth=2,color='grey',label='BC')

    ylabel(r'$\mu$g/m$^3$')
    title('Surface PM2.5 Concentration over Japan')
    grid()
    if Yearly:
        axis([y1,y2,0,16])
        legend(loc='upper center',prop={'size':10})
    else:
        legend(loc='upper right',prop={'size':10})

    if showKD:
        ax2 = twinx()
        p = c.all.astype('float')/1000.
        ax2.plot(c.tyme[I],p[I],'k-',linewidth=1,label='KD Cases')
        ylabel('Average No. KD Cases [K]')
        axis([y1,y2,0,1.5])
        legend(loc='upper left',prop={'size':10})

    savefig('m2.%s.sfc_s.'%tag+trange+'.png',dpi=120,bbox_inches='tight')

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

    fill_between(tyme,0,so4,color='magenta')
    fill_between(tyme,0,dus,color='orange')
    fill_between(tyme,0,sss,color='blue')
    fill_between(tyme,0,ocs,color='brown')
    fill_between(tyme,0,bcs,color='grey')

    plot(tyme,so4,'m-',linewidth=2,label=r'SO$_4$')
    plot(tyme,dus,'-',linewidth=2,color='orange',label='Dust')
    plot(tyme,sss,'b-',linewidth=2,label='Sea Salt')
    plot(tyme,ocs,'-',linewidth=2,color='brown',label='OC')
    plot(tyme,bcs,'-',linewidth=2,color='grey',label='BC')

    title('Aerosol Optical Depth over Japan')
    grid()
    #axis([0,11,0,0.42])
    legend(loc='upper right',prop={'size':10})

    if showKD:
        ax2 = twinx()
        p = c.all.astype('float')/1000.
        ax2.plot(c.tyme[I],p[I],'k-',linewidth=1,label='KD Cases')
        ylabel('Average No. KD Cases [K]')
        axis([y1,y2,0,1.5])
        legend(loc='upper left',prop={'size':10})

    savefig('m2.%s.aod_s.'%tag+trange+'.png',dpi=120,bbox_inches='tight')

#....................................................................

if __name__ == "__main__":

    y1, y2 = (1980,2015)
    trange = '%d-%d'%(y1,y2)
    
    clm = _getYearly(y1,y2)
    plot_Japan_stack(trange,clm,showKD=True,Yearly=True)

else:

    # Instantiate GrADS
    # -----------------
    
    y1, y2 = (2000,2015)
    y1, y2 = (1980,2015)
    trange = '%d-%d'%(y1,y2)

    showKD = False

    d = dict( y1=y1, y2=y2,
              fvIn='/discover/nobackup/projects/gmao/share/dasilva/fvInput/g5chem/sfc',
              ctlPath = '/home/adasilva/opendap/m2/opendap' )

    clm = dict()

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
           open $ctlPath/tavgM_2d_aer_Nx
           set dfile 2 

           set x 1
           set y 1
           set time jan$y1 dec$y2

           define $q = tloop(aave(${q}*if(countries.1(t=1),==,200,1,-u),lon=128,lon=148,lat=30,lat=46))

          """,d)

        clm[q] = ga.expr(q).data

    # Save it for later
    # -----------------
    savez('NPZ/m2_japan.monthly.%s.npz'%trange,**clm)    

    # Plot climatology
    # plot_clm2('%d-%d'%(y1,y2),clm,showKD)
    # plot_clm_stack('%d-%d'%(y1,y2),clm,showKD)
