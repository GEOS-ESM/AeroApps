"""

Plot extinction curtains.

"""

from flightplan import FlightPlan
from nav        import DC8, ER2

from numpy             import tile, linspace, array, ma
from matplotlib.pyplot import contourf, xlabel, ylabel, title, grid, plot, \
                              figure, gca, clf, cm, savefig, axes, colorbar, legend
from types import *
import sys
import unicodedata

UNDEF = 1.0E+10

class Curtain(object):

    def addVar(self,url,Vars=None,Verbose=True,Levels=None):
        """
        Sample variable along flight path.
        """
        from grads import GrADS
        ga = GrADS(Window=False,Echo=False)
        fh = ga.open(url)
        if Levels is not None:
            ga('set lev %s'%Levels)

        if Vars is None:
            Vars = ga.query('file').vars
        elif type(Vars) is StringType:
            Vars = [Vars,]
        for var in Vars:
            if Verbose:
                print ' Working on <%s>'%var
            q = ga.sampleXYT(var,self.Longitude,self.Latitude,self.Tyme,Verbose=Verbose).data
            self.__dict__[var] = ma.MaskedArray(q,mask=q>=UNDEF)
               
    def sampleExt(self,asm_Np,asm_Nv,ext_Np,flx_Nx):
        """
        Sample GEOS-5 Variables at flight path.
        """
        # Retrieve the data & interpolate to flight path
        # -----------------
        self.addVar(flx_Nx,('pblh',))
        self.addVar(asm_Nv,('phis',))
        self.addVar(asm_Np,('h',),Levels='1000 150')
        self.addVar(ext_Np,('duext','ssext','bcext','ocext','suext'),Levels='1000 150')

        # Derived quantities
        # ------------------
        self.pblz = (self.phis/9.8 + self.pblh)/1000.
        self.ext  = self.duext+self.ssext+self.bcext+self.ocext+self.suext
    
    def sampleAer(self,asm_Nv,aer_Nv,flx_Nx):
        """
        Sample GEOS-5 Variables at flight path.
        """
        # Retrieve the data & interpolate to flight path
        # -----------------
        self.addVar(flx_Nx,('pblh',))
        self.addVar(asm_Nv,('phis',))
        self.addVar(asm_Nv,('h',),Levels='72 40')
        self.addVar(aer_Nv,Levels='72 40')

        # Derived quantities
        # ------------------
        self.pblz = (self.phis/9.8 + self.pblh)/1000.

        ####self.ext  = self.duext+self.ssext+self.bcext+self.ocext+self.suext
    

class CurtainFP(FlightPlan,Curtain):
    """
    Produce Curtain plots from ASCII Flight Plans.
    """
    def contourf(self,q,Title=None,Alt=False,N=None,clevs=None,vmaxmin=None,figFile=None,**kwopts):
        """
        Plots a curtain contour plot of time vs. height.
        Input data assumed to be in pressure coordinates
        If N is given, plot N contours.
        If clevs is given, plot contours for the contour levels in the list.
        If vmaxmin = (vmin, vmax) is given, plot linearly spaced levels between
        vmin and vmax.
        """

        nt, nz = q.shape
        nh = self.h.shape[1]
        
        if self.aircraft == 'DC8':
            Alt = True # Always plot altitude for DC8
        
        clf()
        ax = gca()
        ax.set_axis_bgcolor('black') 
        
        h = self.h.mean(axis=0)/1000

        # Hack to cope with the fact that asm_Np has 42 levs while
        # ext_Np has 26.
        # --------------------------------------------------------
        if nh > nz:
            h = _fixLev(h)
                
        if N is None and clevs is None and vmaxmin is None:
            ## Python picks the levels
            contourf(self.Hour,h,q.T,**kwopts)
        elif N != None:
            ## Python picks N levels
            contourf(self.Hour,h,q.T,N,**kwopts)
        elif clevs != None:
            ## clevs gives the levels
            contourf(self.Hour,h,q.T,clevs,**kwopts)
        elif vmaxmin != None:
            ## Linearly spaced levels between vmaxmin[0] and vmaxmin[1]
            if len(vmaxmin) < 2:
                print 'vmaxmin should be a 2-element list'
                sys.exit()
            clevs = linspace(vmaxmin[0],vmaxmin[1],32)
            contourf(self.Hour,h,q.T,**kwopts)

        _colorbar()
        if Alt:
            plot(self.Hour,self.Altitude,'w',linewidth=2,label=self.aircraft+' Altitude')
        if 'pblz' in self.__dict__.keys():
            plot(self.Hour,self.pblz,'w--',linewidth=2,label='PBL Height')
        legend(loc='upper right')
            
        grid()
        xlabel('UTC Hour on %s'%self.takeoff.date().ctime().replace('00:00:00 ',''))
        ylabel(' Altitude (km)')
        if Title is not None:
            title(Title)

        if figFile is not None:
            savefig(figFile,dpi=180)
    
class CurtainDC8(DC8,Curtain):
    """
    Produce Curtain plots from ASCII Flight Plans.
    """
    def contourf(self,q,Title=None,Alt=False,N=None,figFile=None,**kwopts):
        """
        Plots a curtain contour plot of time vs. height.
        Input data assumed to be in pressure coordinates
        """

        nt, nz = q.shape
        nh = self.h.shape[1]
        
        if self.aircraft == 'DC8':
            Alt = True # Always plot altitude for DC8
        
        clf()
        ax = gca()
        ax.set_axis_bgcolor('black') 
        
        h = self.h.mean(axis=0)/1000

        # Hack to cope with the fact that asm_Np has 42 levs while
        # ext_Np has 26.
        # --------------------------------------------------------
        if nh > nz:
            h = _fixLev(h)
        
        if N is None:
            contourf(self.Hour,h,q.T,**kwopts)
        else:
            contourf(self.Hour,h,q.T,N,**kwopts)

        _colorbar()
        if Alt:
            plot(self.Hour,self.Altitude,'k',linewidth=2,label=self.aircraft+' Altitude')
        if 'pblz' in f.__dict__.keys():
            plot(self.Hour,self.pblz,'k--',linewidth=2,label='PBL Height')
        legend(loc='upper right')
            
        grid()
        xlabel('UTC Hour on %s'%self.takeoff.date().ctime().replace('00:00:00 ',''))
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
    labels = ['%3.3f' % float(item) for item in labels]
    
    cbar.ax.set_yticklabels(labels)        
    axes(ax)  # make the original axes current again

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




