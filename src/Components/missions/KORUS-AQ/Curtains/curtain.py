"""

Plot extinction curtains.

"""

from pyobs      import ICARTT

from numpy             import tile, linspace, array, ma
from matplotlib.pyplot import contourf, xlabel, ylabel, title, grid, plot, \
                              figure, gca, clf, cm, savefig, axes, colorbar, legend
from types import *
import sys
import unicodedata

from grads import GrADS

from datetime import datetime, timedelta


UNDEF = 1.0E+10

class Curtain(object):
    """
    Produce Curtain plots from ASCII Flight Plans.
    """

    def __init__(self,meteo,chem,ict,aircraft='Aircraft'):
        """
        Load 
        """

        # Opengrads
        # ---------
        ga = GrADS(Echo=False,Window=False)
        self.ga = ga # save for later

        # Load Vertical coordinates
        # -------------------------
        fh_met = ga.open(meteo)
        fh_chm = ga.open(chem)
        ga('set t 1 %d'%fh_met.nt)
        ga('set z 1')
        self.pblz = ga.expr('(pblh + phis/9.81)/1000.') # in km, ASL 
        ga('set z 1 %d'%fh_met.nz)   # met has 26 levs
        self.h = ga.expr('h')
        ga('set dfile 2')
        ga('set z 1 %d'%fh_chm.nz)

 
        # Flight path
        # -----------
        f = ICARTT(ict)
        self.Longitude, self.Latitude, self.tyme = f.Lon, f.Lat, f.tyme
        try:
            self.Altitude = f.AltP   # real flight, in km (from 60s nav)
        except:
            self.Altitude = 1000*f.Alt_km # Flight plan
        t0 = self.tyme[-1]
        t0 = datetime(t0.year,t0.month,t0.day)
        self.Hour = array([(t-t0).total_seconds()/(60.*60.) for t in self.tyme])

        self.takeoff = self.tyme[0]
        self.landing = self.tyme[-1]
        self.aircraft = aircraft.upper()

#---
    def loadExt(self):

        ga = self.ga
        self.duext = ga.expr('duext')
        self.ssext = ga.expr('ssext')
        self.bcext = ga.expr('bcext')
        self.ocext = ga.expr('ocext')
        self.suext = ga.expr('suext')

        self.ccext = self.bcext+self.ocext
        self.ext  = self.duext+self.ssext+self.bcext+self.ocext+self.suext

#---
    def loadConc(self):

        ga = self.ga
        self.du = ga.expr('airdens*du')
        self.ss = ga.expr('airdens*ss')
        self.bc = ga.expr('airdens*bc')
        self.oc = ga.expr('airdens*oc')
        self.su = ga.expr('airdens*so4')
        self.so2 = ga.expr('airdens*so2')

        self.cc = self.bc+self.oc

        self.co2 = ga.expr('airdens*co2')
        self.co = ga.expr('airdens*co')
        self.coffas = ga.expr('airdens*conbas')
        self.cobbae = ga.expr('airdens*cobbae')
        self.cobbot = ga.expr('airdens*(cobbgl-cobbae)')
        
#---
    def contourf(self,q,Title=None,Alt=False,N=None,figFile=None,hmax=9,**kwopts):
        """
        Plots a curtain contour plot of time vs. height.
        Input data assumed to be in pressure coordinates
        """

        # UNDEF safeguard
        I = q.data>0.01*UNDEF # hackish
        q.data[I] = UNDEF
        q.mask[I] = True
        
        nt, nz = q.shape
        nh = self.h.shape[1]
        
        clf()
        ax = gca()
        ax.set_axis_bgcolor('black') 
        
        h = self.h.mean(axis=0)/1000.
        h[0] = 0.

        # Hack to cope with the fact that asm_Np has 42 levs while
        # ext_Np has 26.
        # --------------------------------------------------------
        if nh > nz:
            h = _fixLev(h)
            I = h<=hmax
        
        if N is None:
            contourf(self.Hour,h[I],q[:,I].T,**kwopts)
        else:
            contourf(self.Hour,h[I],q[:,I].T,N,**kwopts)

        _colorbar()
        
        plot(self.Hour,self.Altitude,'m',linewidth=2,label=self.aircraft+' Altitude')
        if 'pblz' in self.__dict__.keys():
            plot(self.Hour,self.pblz,'k-',linewidth=2,label='PBL Height')
        legend(loc='upper right')
            
        grid()
        xlabel('UTC Hour on %s'%self.landing.date().ctime().replace('00:00:00 ',''))
        ylabel(' Altitude (km)')
        if Title is not None:
            title(Title)

        if figFile is not None:
            savefig(figFile,dpi=180)
    

#---------------------------------------------------------------------------------
def _colorbar():                                                                
    """ Draw a colorbar """                                                     
                                                                                
    ax = gca()                                                                  
    bbox = ax.get_position()                                                    
    l,b,w,h = bbox.bounds # mpl >= 0.98                                         
    cax = axes([l+w+0.01, b, 0.04, h]) # setup colorbar axes.
    cbar = colorbar(cax=cax)
    labels = [unicodedata.normalize('NFKD', item.get_text()).encode('ascii','ignore') for item in cbar.ax.get_yticklabels()]
    #labels = ['%3.3f' % float(item) for item in labels]
    labels = ['%g' % float(item) for item in labels]
    
    cbar.ax.set_yticklabels(labels)        
    axes(ax)  # make the original axes current again

    #---
def _fixLev(h):
    
    lev26 = [1000, 975, 950, 925, 900,      850,      800,      750,      700, 650,
              600, 550, 500, 450, 400, 350, 300, 250, 200, 150, 100, 70, 50,     30, 20, 10 ]

    lev42 = [1000, 975, 950, 925, 900, 875, 850, 825, 800, 775, 750, 725, 700, 650, 
              600, 550, 500, 450, 400, 350, 300, 250, 200, 150, 100, 70, 50, 40, 30, 20, 10,
              7, 5, 4, 3, 2, 1, 0.7, 0.5, 0.4, 0.3, 0.1 ]

    h_ = []

    for i in range(len(h)):
        if lev42[i] in lev26:
            h_ += [h[i],]

    h_ = array(h_)

    return h_
            
#---------------------------------------------------------------------------------------

if __name__ == "__main__":

    sdir = '/home/adasilva/iesa/kaq/sampled/plan'
    meteo =  sdir + '/KORUSAQ-GEOS5-METEO-DC8_PLAN_20160503_R0.nc'
    chem  =  sdir + '/KORUSAQ-GEOS5-CHEM-DC8_PLAN_20160503_R0.nc'
    ict = '../Plans/fltplan_dc8_20160503.ict'

    c = Curtain(meteo,chem,ict)
