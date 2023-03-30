"""
Reads Level 2 MOD14/MYD14 granules for a single day and returns a
single object with the relevant data.

This software is hereby placed in the public domain.
Arlindo.daSilva@nasa.gov
"""

import os
from types import *
from numpy import zeros, ones, sqrt, std, mean, unique, concatenate, where, array, load, savez

from pylab import prctile, linspace, plot, xlabel, ylabel, title

from pyhdf.SD import *

from scipy.stats import kde

from grads import GrADS

from datetime  import date, timedelta, datetime

from glob import glob

DAY = timedelta(seconds=60*60*24)

#....................................................................
class Sample(object):
    def __init__(self):
        return

class MxD14_L2(object):
    """
    This class implements the MODIS Level 2 fire products, usually
    referred to as MOD14 (TERRA satellite) and MYD14 (AQUA satellite).
    """

    def __init__ (self,Path,RoundHour=False,Verb=0,qc_thresh=0.0,dt_thresh=-999.):
       """
       Reads individual granules or a full day of Level 2 MOD14/MYD14 files
       present on a given *Path* and returns a single object with
       all data concatenated. On input, 

         Path -- can be a single file, a single directory, of a list
                 of files and directories.  Directories are
                 transversed recursively. If a non MOD14/MYD14 Level 2
                 file is encountered, it is simply ignored. 
         RoundHour -- GMT and local times should be rounded to whole hours
         Verb -- Verbose level:
                 0 - really quiet (default)
                 1 - Warns if invalid file is found
                 2 - Prints out non-zero number of fires in each file.

       """

#      Initially are lists of numpy arrays for each granule
#      ------------------------------------------------
       self.verb = Verb
       self.algo = None # fire property algorithm; set by DOZIER
       self.sat  = None # Satellite name
       self.col  = None # collection, e.g., 005
       self.lon  = [] # longitude (degrees)
       self.lat  = [] # latitude (degrees)
       self.pixar= [] # pixel area (km2)
       self.pow  = [] # fire radiative power (MW)
       self.R2   = [] # reflectance of channel 2 (0.858 microns)
       self.T21  = [] # brightness temperature for MODIS channel 21 (4 microns)
       self.T31  = [] # brightness temperature for MODIS channel 31 (11 microns)
       self.DT   = [] # mean backround bright temp difference
       self.Tb21 = [] # background brightness temperature channel 21 (4 microns)
       self.Tb31 = [] # background brightness temperature channel 31 (11 microns)
       self.qc   = [] # detection confidence (%)
       self.yyyy = [] # Year
       self.jjj  = [] # Julian day (within year)
       self.hh  =  [] # hour of granule
       self.nn  =  [] # minute of granule
       self.tgmt = [] # GMT time (deribed from hh & nn)
       self.tloc = [] # Local time
       self.tga  = [] # GrADS time string

       self.met = None    # Met Fields
       self.sample = None # will replace .met 
       self.h_F = None # flaming heat flux
       self.r_F = None # fraction of flaming energy

#      Stop here if empty object derired
#      ---------------------------------
       if Path is None:
           return
       
#      Read each granule, appending them to the list
#      ---------------------------------------------
       if type(Path) is not ListType:
           Path = sorted(glob(Path)) # expand pattern if need be
       self._readList(Path)

#      Make each attribute a single numpy array, screening for qc=0
#      ------------------------------------------------------------
       self.qc   = concatenate(self.qc)
       self.T21 = concatenate(self.T21)
       self.T31 = concatenate(self.T31)
       self.DT = concatenate(self.DT)
       self.Tb21 = concatenate(self.Tb21)
       self.Tb31 = concatenate(self.Tb31)
       m = (self.qc > qc_thresh) & (self.Tb21 > 0.0) & (self.Tb31 > 0.0) & \
           ((self.T31-self.Tb31) > dt_thresh)
       nraw = self.qc.size
       self.qc  = self.qc[m]
       nrej = nraw - self.qc.size
       self.Tb21 = self.Tb21[m]
       self.Tb31 = self.Tb31[m]
       self.T21 = self.T21[m]
       self.T31 = self.T31[m]
       self.DT = self.DT[m]
       self.lon = concatenate(self.lon)[m]
       self.lat = concatenate(self.lat)[m]
       self.pixar = concatenate(self.pixar)[m] # in km2
       self.pow = concatenate(self.pow)[m]
       self.R2  = concatenate(self.R2)[m]
       self.yyyy = concatenate(self.yyyy)[m]
       self.jjj  = concatenate(self.jjj)[m]
       self.hh   = concatenate(self.hh)[m]
       self.nn   = concatenate(self.nn)[m]
       self.tgmt = concatenate(self.tgmt)[m]
       self.tloc = concatenate(self.tloc)[m]
       self.tga = concatenate(self.tga)[m]

#      Ensure type of some attributes
#      ------------------------------
       self.yyyy = self.yyyy.astype(int)
       self.jjj  = self.jjj.astype(int)
       self.hh   = self.hh.astype(int)
       self.nn   = self.nn.astype(int)

       # Build datetime for each fire
       # ----------------------------
       self.tyme = array([gat2dt(gat) for gat in self.tga])
       
#      Compute pixel area from sample number
#      -------------------------------------
       self.pixar = _pixar(self.pixar)

       self.m = (self.qc>0) # dozier will redefine it
       
       if self.verb>0:
           print_stats('__header__','Summary for %s with %d Fires, %d Rejected (%4.2f%%)'\
                       %(Path[0],nraw,nrej,100.*nrej/nraw) )
           print_stats('PixArea',self.pixar)
           print_stats('R2  ',self.R2)
           print_stats('pow ',self.pow)
           print_stats('T21 ',self.T21)
           print_stats('Tb21',self.Tb21)
           print_stats('DT21',self.T21-self.Tb21)
           print_stats('T31 ',self.T31)
           print_stats('Tb31',self.Tb31)
           print_stats('DT31',self.T31-self.Tb31)
           print_stats('DT',self.DT)
           print_stats('Q/C ',self.qc)
           print_stats('__footer__')

#---
    def _readList(self,List,RoundHour=False):
        """
        Recursively, look for files in list; list items can
        be files or directories.
        """
        for item in List:
            if os.path.isdir(item):      self._readDir(item,RoundHour)
            elif os.path.isfile(item):   self._readGranule(item,RoundHour)
            else:
                print("%s is not a valid file or directory, ignoring it"%item)
#---
    def _readDir(self,dir,RoundHour=False):
        """Recursively, look for files in directory."""
        for item in os.listdir(dir):
            path = dir + os.sep + item
            if os.path.isdir(path):      self._readDir(path,RoundHour)
            elif os.path.isfile(path):   self._readGranule(path,RoundHour)
            else:
                print("%s is not a valid file or directory, ignoring it"%item)

#---
    def _readGranule(self,filename,RoundHour=False):
        """Reads one MOD14/MYD14 granule with Level 2 fire data."""

#       Don't fuss if the file cannot be opened
#       ---------------------------------------
        try:
            hfile = SD(filename)
        except HDF4Error:
            if self.verb > 2:
                print("- %s: not recognized as an HDF file"%filename)
            return 

#       No fires, nothing to do
#       ----------------------
        if hfile.select('FP_longitude').checkempty():
            if self.verb > 2:
                print("- %s:  no  fires"%filename)
            return

#       Read select variables
#       ---------------------
        self.lon.append(hfile.select('FP_longitude').get())
        self.lat.append(hfile.select('FP_latitude').get())
        self.pixar.append(hfile.select('FP_sample').get())
        self.pow.append(hfile.select('FP_power').get())
        self.R2.append(hfile.select('FP_R2').get())
        self.T21.append(hfile.select('FP_T21').get())
        self.T31.append(hfile.select('FP_T31').get())
        self.DT.append(hfile.select('FP_MeanDT').get())
        self.Tb21.append(hfile.select('FP_MeanT21').get())
        self.Tb31.append(hfile.select('FP_MeanT31').get())
        self.qc.append(hfile.select('FP_confidence').get())

#       Satellite name
#       --------------
        str = hfile.attributes()['MOD03 input file']
        if self.sat is None:
            sat = str.split('.')[0].split('/')[-1][0:3]
            if sat == "MYD":
                sat = 'Aqua'
            elif sat == "MOD":
                sat = 'Terra'
            else:
                raise ValueError('Unknown MOD03 file type: '+sat)

            self.sat = sat

#       Collection
#       ----------
        if self.col is None:
            i = str.index('D03.A')-2
            base = str[i:i+23]   # e.g., MOD03.A2003001.0000.005
            self.col = base.split(".")[3]

#       Granule time
#       ------------
        temp1=str.split('.')[1].split('-')[0][1:12]
        temp2=str.split('.')[2].split('.')[0][0:4]
        timestamp = temp1+temp2

        yyyy = int(timestamp[0:4])
        jjj = int(timestamp[4:7])
        hh = int(timestamp[7:9])
        nn = int(timestamp[9:11])

        tgmt = _local_time(hh,nn,zeros(1),RoundHour) # Local time
        tga = [_gatime(yyyy,jjj,hh,nn),] # GrADS time string
        tloc = _local_time(hh,nn,self.lon[-1],RoundHour) # Local time

        n = self.lon[-1].size

        self.yyyy.append(yyyy * ones(n))
        self.jjj.append(jjj * ones(n))
        self.hh.append(hh * ones(n))
        self.nn.append(nn * ones(n))

        self.tgmt.append(tgmt * ones(n))
        self.tga.append(n*tga)
        self.tloc.append(tloc)
    
        if self.verb > 1:
            print("- %s: %4d fires"%(filename,n))

#---

    def attach(self,filename,Vars=None,Cache=False):
        """
        Attach variables from gridded meteorological files, interpolating
        from the gridded values to the (lat,lon) of the fire. If the variables
        are 3D, the whole curtain is included. If *Vars* is not specified,
        all non-coordinate variables on file are included. 

        When *Cache* is true, the variables are saved locally to directory
        "__cache__". When the input file name is "__cache__", the variables are read
        from cache.

        Note: this method is being deprecated in favor o sampleFile().
        
        """

#       Expects GrADS v2
#       ----------------
        ga = GrADS(Bin='grads',Window=False,Echo=False)

#       Open the file
#       -------------
        fh = ga.open(filename)
        ga('set lon -180 180') # good for interpolation

#       Either all varibles on file or user subset
#       ------------------------------------------
        if Vars == None:
            vinfo = fh.var_info
        else:
            vinfo = []
            for v,k,l in fh.var_info:
                if v in Vars:
                    vinfo.append((v,k,l))
            if len(vinfo)==0:
                print("IndexError: requested variables - ", Vars)
                raise IndexError("cannot find any matchig variable in file %f"\
                    %filename)

#       For each observation, find the correspondng time on file
#       --------------------------------------------------------
        utimes = unique(self.tga)   # unique obs times in grads format
        self.tgaf = self.tga.copy() # will hold times on file for each ob
        for tga in utimes:
            ga('set time %s'%tga,Quiet=True)
            qh = ga('query time',Quiet=True)
            self.tgaf[self.tga==tga] = ga.rword(1,3)

#       Loop over each desired variable and interpolate it to
#        fire location 
#       -----------------------------------------------------
        self.met = {}
        n = self.lon.size
        levs = 1000. * zeros(fh.nz)
        for v,nlevs,l in vinfo:
            if nlevs==0: 
                nlevs=1
                y_f = zeros(n)
            else:
                y_f = zeros((n,nlevs))
            for tgaf in unique(self.tgaf): 
                ga('set time %s'%tgaf,Quiet=True)
                ga('set z %d %d'%(1,nlevs),Quiet=True)
                m = (self.tgaf == tgaf) # gather obs for this time
                lon_, lat_ = (self.lon[m], self.lat[m])
                print("- Interpolating %5d %s obs at %s"%(lon_.size,v,tgaf))
                y_f[m], levs = ga.interp(v,lon_,lat_) # interp & scatter
            self.met[v] = y_f
        self.met['levs'] = levs # record vertical levels

#---
    def sampleFP(self, npzFile=None, dir='/home/adasilva/iesa/aerosol/experiments/seac4rs_01',I=None):
        """
        Get Met fields from GEOS-5 Forward Processing. Variable selection customized for
        PR modeling.
        """
        from grads import GrADS

        ga = GrADS(Window=False,Echo=False)
    
        ga.open(dir+'/tavg1_2d_flx_Nx.ddf')
        self.sampleFile(ga,onlyVars=('ustar', 'bstar', 'pblh', 'hflux', 'eflux', 'rhoa', 'tsh'),I=I)
        ga('close 1')

        ga.open(dir+'/inst3_3d_asm_Nv.ddf')
        ga('set z 1 72')
        self.sampleFile(ga,onlyVars=('phis','u','v','t','qv','o3','delp'),
                        npzFile=npzFile, I=I)
        ga('close 1')

#---
    def sampleFile(self, ga, npzFile=None, onlyVars=None, I=None, Verbose=True):
        """
        Interpolates all variables of default file in ga and optionally
        save them to file *npzFile*
        """
        
        if self.sample == None:
            self.sample = Sample() # first time
            
        if onlyVars is None:
            fh = ga.query('file',Quiet=True)
            onlyVars = fh.vars

        if I is None:
            I = list(range(len(self.lon)))

        self.sample.I = I  # save for consistency check
        
        tymes = self.tyme[I]
        lons = self.lon[I]
        lats = self.lat[I]

        # Loop over variables on file
        # ---------------------------
        for v in onlyVars:
            if Verbose:
                print("<> Sampling ", v)
            var = ga.sampleXYT(v,lons,lats,tymes,Verbose=Verbose)
            self.sample.__dict__[v] = var.data

        if npzFile is not None:
            savez(npzFile,lons=lons,lats=lats,tymes=tymes,**self.sample.__dict__)            
       
#---
    def sampleLoadz(self,npzFile):
        """
        Loads sample from npz file.
        """
        self.sample = Sample()
        npz = load(npzFile)
        for v in list(npz.keys()):
            self.sample.__dict__[v] = npz[v]

#----
    def checkpoint(self,npzFile):
        """
        Save attributes to disk file.
        """
        savez(npzFile,**self.__dict__)

    def restart(self,npzFile):
        """
        Retrieve attributes from disk file.
        """
        f = load(npzFile)
        for v in list(f.keys()):
            self.__dict__[v] = f[v]

#......................................................................................

def granules(t1,t2,product='MOD14',
             rootdir='/home/adasilva/iesa/aerosol/data/MODIS/Level2'):
    """
    Return the granules, given a range of dates.
    """
    oneday = timedelta(hours=24)
    d, d2 = t1.date(), t2.date()
    Files = []
    while d<=d2:
        year = d.year
        doy =  (d - date(year-1,12,31)).days
        g = rootdir+'/%s/%d/%03d/%s.*.*.*.hdf'%(product,year,doy,product)
        Files += sorted(glob(g))
        d += oneday
    return Files

__Months__ = ['JAN','FEB','MAR','APR','MAY','JUN','JUL','AUG','SEP','OCT','NOV','DEC']

def gat2dt(gat):
    """
    Convert grads time to datetime.
    """
    time, date = gat.upper().split('Z')
    if time.count(':') > 0:
        h, m = time.split(":")
    else:
        h = time
        m = '0'
    mmm = date[-7:-4]
    dd, yy = date.split(mmm)
    mm = __Months__.index(mmm) + 1
    dt = datetime(int(yy),int(mm),int(dd),int(h),int(m))
    return dt

def dt2gat(t):
    """
    Convert datetime to grads time.
    """
    gat = "%d:%dZ%d%s%d"%(t.hour,t.minute,t.day,__Months__[t.month-1],t.year)
    return gat

def _local_time(hh,nn,lon,RoundHour=False):
    tlocal =  float(hh) + (float(nn)/60.) + (12./180.) * lon
    if RoundHour:
        tlocal = floor(tlocal+0.5) # round up to nearest hour
    tlocal = where(tlocal<  0, tlocal+24, tlocal)
    tlocal = where(tlocal>=24, tlocal-24, tlocal)
    return tlocal

def _gatime(yyyy,jjj,hh=None,nn=None):
    """Forms a GrADS time string from year, julian day, hour and minute"""
    t = date(yyyy,1,1) + (jjj - 1)*DAY
    if hh is None:
        return '%02d%s%4d'%(t.day,t.ctime()[4:7],yyyy)
    else:
        return '%02d:%02dZ%02d%s%4d'%(hh,nn,t.day,t.ctime()[4:7],yyyy)

def _nint(str):
    """Nearest integer, internal use."""
    x = float(str)
    if x >= 0: return int(x+0.5)
    else:      return int(x-0.5)

def _pixar(x):
    """
    Compute pixel area (in km2) given the sample number *x*. Polynomial approximation
    given in p. 34 of MODIS Fire Product User's Guide.
    """

#   Polynomial coeffs
#   -----------------
    c = [9.7421684,-0.091159223,0.00051138175,-1.7683231e-6,3.8048273e-9,
        -5.0660609e-12, 4.0471196e-15, -1.7739490e-18, 3.2795410e-22]

#   Horner's rule
#   -------------
    area = c[8] * ones(x.size)
    for i in (7,6,5,4,3,2,1,0):
        area = area * x + c[i]

    return area

def do_kde(X,range=None,N=256):
    if range is None:
        prc = prctile(X.ravel())
        a = prc[0]
        b = prc[4]
    else:
        a, b = range
    bins = linspace(a,b,N)
    kernel = kde.gaussian_kde(X.ravel())
    return bins, kernel(bins)

def plot_kde(X,a=None,b=None,N=256,Title=None,Label=None):
    if a==None:
        prc = prctile(X.ravel())
        a = prc[0]
        b = prc[4]
    if Title is None: 
        Title = 'Kernel Density Function'
    bins = linspace(a,b,N)
    kernel = kde.gaussian_kde(X.ravel())
    plot(bins,kernel(bins))
    ylabel('PDF')
    title(Title)
           
def print_stats(name,x=None):
    "Prints simple stats"
    if type(name) is not StringType:
        x = name
        name = 'mean,stdv,rms,min,25%,median,75%,max: '
    if name == '__header__':
        print('')
        n = (80 - len(x))/2
        print(n * ' ' + x)
        print(n * ' ' + len(x) * '-')
        print('')
        print('   Name       mean      stdv      rms      min     25%    median     75%      max')
        print(' ---------  -------  -------  -------  -------  -------  -------  -------  -------')
    elif name == '__sep__':
        print(' ---------  -------  -------  -------  -------  -------  -------  -------  -------')
    elif name == '__footer__':
        print(' ---------  -------  -------  -------  -------  -------  -------  -------  -------')
        print('')
    else:
        ave = x.mean()
        std = x.std()
        rms = sqrt(ave*ave+std*std)
        prc = prctile(x)
        print('%10s  %7.2f  %7.2f  %7.2f  %7.2f  %7.2f  %7.2f  %7.2f  %7.2f  '%\
            (name,ave,std,rms,prc[0],prc[1],prc[2],prc[3],prc[4]))


#............................................................................

if __name__ == "__main__":

    t1 = datetime(2013,8,2)
    t2 = datetime(2013,8,31)

    for p in ( 'MOD14', 'MYD14'):
        
        Files = granules(t1,t2,product=p)
    
        f = MxD14_L2(Files,Verb=1)

        I_na = (f.lon>-170)&(f.lon<-50)&(f.lat>15)&(f.lat<80)
        
        f.sampleFP(npzFile='%s.sample_nam.2013-08.npz'%p,I=I_na)
