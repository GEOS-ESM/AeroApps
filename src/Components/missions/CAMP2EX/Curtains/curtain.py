"""

Plot extinction curtains.

"""

from pyobs.icartt import ICARTT

from numpy import tile, linspace, array, ma, ones, zeros, interp, \
                        isnan, NaN, savez
from matplotlib.pyplot import contourf, xlabel, ylabel, title, grid, plot, \
                              figure, gca, clf, cm, savefig, axes, colorbar, legend
from types import *
import sys
import unicodedata
from csv import DictReader

from grads import GrADS

from datetime import datetime, timedelta
from dateutil.parser import parse as isoparser

UNDEF = 1.0E+10

def _subsetLev(h42):
    
    lev26 = [1000, 975, 950, 925, 900,      850,      800,      750,      700, 650,
              600, 550, 500, 450, 400, 350, 300, 250, 200, 150, 100, 70, 50,     30, 20, 10 ]

    lev42 = [1000, 975, 950, 925, 900, 875, 850, 825, 800, 775, 750, 725, 700, 650, 
              600, 550, 500, 450, 400, 350, 300, 250, 200, 150, 100, 70, 50, 40, 30, 20, 10,
              7, 5, 4, 3, 2, 1, 0.7, 0.5, 0.4, 0.3, 0.1 ]

    nt = h42.shape[0]
    h26 = ma.ones((nt,26),fill_value=UNDEF)
    
    j = 0
    for i in range(42):
        if lev42[i] in lev26:
            h26[:,j] = h42[:,i]
            j += 1

    return h26

def _zInterp(z,h,v):
        nz = len(z)
        nt, nh = h.shape
        x = zeros((nt,nz))+UNDEF
        x = ma.masked_array(x,isnan(x),fill_value=UNDEF)
        for t in range(nt):
            hmin = h[t,:].min()
            #I = z>hmin
            J = v.mask[t,:]==False
            #x.data[t,I] = interp(z[I],h[t,J],v[t,J],left=UNDEF)
            x.data[t,:] = interp(z,h[t,J],v[t,J],left=UNDEF)
            #print t, hmin, x[t,:].min(), v[t,:].min()

        x.mask[x.data==UNDEF] = True

        return x.T
            
def _getTrackCSV(csvFile):
    """
    Get trajectory from a CSV with (lon,lat,time) coordinates.
    Notice that *time* must be specified in ISO format, e.g.,

                  2014-02-05T12:30:45
    """
    CSV = DictReader(open(csvFile))
    lon, lat, tyme = [], [], []
    for row in CSV:
        lon  += [float(row['lon']),]
        lat  += [float(row['lat']),]
        tyme  += [isoparser(row['time']),]
        
    return ( array(lon), array(lat), array(tyme) )


class Curtain(object):
    """
    Produce Curtain plots from ASCII Flight Plans.
    """

    def __init__(self,meteo,chem,ext,coords,aircraft='Aircraft',zmax=8,nz=160,prs=True):
        """
        Load 
        """

        # Opengrads
        # ---------
        ga = GrADS(Echo=False,Window=False)
        self.ga = ga # save for later

        # Open files
        # ----------
        self.fh_met = ga.open(meteo)
        self.fh_chm = ga.open(chem)
        self.fh_ext = ga.open(ext)
        
        # Load Vertical coordinates
        # -------------------------
        ga('set t 1 %d'%self.fh_met.nt)
        ga('set z 1')
        self.pblz = ga.expr('(pblh + phis/9.81)/1000.') # in km, ASL 
        self.hs = ga.expr('(phis/9.81)/1000.') # in km, ASL 
        ga('set z 1 %d'%self.fh_met.nz)   # met has 26 levs
        self.H_met = ga.expr('h/1000.')
        ga('set dfile 2')
        ga('set z 1 %d'%self.fh_chm.nz)

        self.nt, self.nh_met = self.fh_met.nt, self.fh_met.nz
        self.nh_chm = self.fh_chm.nz

        # Created reduced H for chem file
        # -------------------------------
        if prs:
            self.H_chm = _subsetLev(self.H_met)
        else:
            self.H_chm = self.H_met # native vertical coords
            

        # Constant height grid
        # --------------------
        self.z = linspace(0,zmax,nz)
        self.nz = nz
 
        # Flight coordinates
        # ------------------
        if coords[-4:] == '.ict' or coords[-7:] == '.ict.gz':
            f = ICARTT(coords)
            self.Longitude, self.Latitude, self.tyme = f.Nav['Longitude'][::60], f.Nav['Latitude'][::60], f.tyme[::60]
            self.Altitude = f.Nav['Altitude'][::60]/1000. # km
        else:
            raise ValueError('invalid coords file %s'%coords)

        #print 'Altitude',self.Altitude

        # UTC hour
        # --------
        t0 = self.tyme[-1]
        t0 = datetime(t0.year,t0.month,t0.day)
        self.Hour = array([(t-t0).total_seconds()/(60.*60.) for t in self.tyme])

        self.takeoff = self.tyme[0]
        self.landing = self.tyme[-1]
        self.aircraft = aircraft.upper()

#---
    def loadMet(self):

        ga = self.ga
        ga('set dfile %d'%self.fh_met.fid)
        ga('set z 1 %d'%self.fh_met.nz)
        self.T  = ga.expr('T')
        self.RH = ga.expr('100*RH')
        try:
            self.CLOUD = ga.expr('cloud')
        except:
            pass

    def load(self,fh,Vars,factor=None):
        """
        Load a list of variables:
        """
        ga = self.ga
        ga('set dfile %d'%fh.fid)
        ga('set z 1 %d'%fh.nz)
        for v in Vars:
            if factor is not None:
                self.__dict__[v] = ga.expr(factor+'*'+v)
            else:
                self.__dict__[v] = ga.expr(v)
                        
#--
    def zInterp(self,v5):
        """
        Vertically interpolates the GEOS-5 variable v5
        to fixed heights.
        """
        nh = v5.shape[1]
        if nh == self.nh_chm:
            return _zInterp(self.z,self.H_chm,v5)
        elif nh == self.nh_met:
            return _zInterp(self.z,self.H_met,v5)
        else:
            raise ValueError("invalid vertical dimension, nh=%d"%nh)
            
#---
    def contourf(self,q,Title=None,Alt=False,N=None,figFile=None,pblc='m',**kwopts):
        """
        Plots a curtain contour plot of time vs. height.
        """

        clf()
        ax = gca()
        ax.set_axis_bgcolor('black') 
        

        if N is None:
            contourf(self.Hour,self.z,self.zInterp(q),**kwopts)
        else:
            contourf(self.Hour,self.z,self.zInterp(q),N,**kwopts)

        _colorbar()

        if self.aircraft.upper() != 'ER2':
            plot(self.Hour,self.Altitude,'b',linewidth=2,label=self.aircraft+' Altitude')
        if 'pblz' in list(self.__dict__.keys()):
            plot(self.Hour,self.pblz,pblc+'-',linewidth=2,label='PBL Height')
        legend(loc='upper right')
            
        #grid(color='w')
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
    cbar.ax.set_axis_bgcolor('white') 
    axes(ax)  # make the original axes current again

    #---

    #---------------------------------------------------------------------------------------

if __name__ == "__main__":

    sdir = '/home/adasilva/iesa/kaq/sampled/plan'
    meteo =  sdir + '/KORUSAQ-GEOS5-METEO-DC8_PLAN_20160503_R0.nc'
    chem  =  sdir + '/KORUSAQ-GEOS5-CHEM-DC8_PLAN_20160503_R0.nc'
    ict = '../Plans/fltplan_dc8_20160503.ict'

    c = Curtain(meteo,chem,ict)
