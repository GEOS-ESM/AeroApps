"""
Implements interface to Pete's NetCDF Level 2m files for the LIDAR crowd. 

"""

from netCDF4 import Dataset 

from datetime import datetime, timedelta
from matplotlib.dates import num2date, date2num
from numpy import array, savez, load

SDS = dict (
      COORD = ('longitude', 'latitude', 'h', 'h0', 'delp', 'time' ),
      METEO = ('u', 'v', 't', 'rh','ps' ),
      AERO  = ('du','ss','bc','oc','SO4','SO2','extinction', 'ssa', 'backscat', 'tau', 'aback_sfc', 'aback_toa', 'ext2back'),
    )

ALIAS = dict (
    longitude = 'lon',
    latitude = 'lat',
    )

class LIDAR_L2(object):
    """
    Simple interface to LIDAR Level 2m files.
    """

    def __init__(self,filename,Verbose=False):
        """
        Ingests may attributes from filename.
        """

       
        # Initially are lists of numpy arrays for each granule
        # ------------------------------------------------
        self.SDS = SDS['COORD'] + SDS['METEO'] + SDS['AERO']

        # Read select variables
        # ---------------------
        f = Dataset(filename)
        for sds in self.SDS:
            v = f.variables[sds]
            if sds in ALIAS:
                sds = ALIAS[sds]
            self.__dict__[sds] = v 

        self.nt, self.nz =self.h.shape
        

        # handle time
        # -----------
        self.tyme = []
        
        for t in self.time[:]:
            nymd = int(t)
            yy = int(nymd/10000)
            mm = int((nymd - 10000*yy)/100)
            dd = nymd - (yy*10000 + mm * 100) 
            fday = t - nymd
            self.tyme.append(datetime(yy,mm,dd)+timedelta(seconds=int(fday*86400)))
            
        self.tyme = array(self.tyme)
        self.tnum = date2num(self.tyme)
       
    def sampleFile(self, inFile, npzFile=None, onlyVars=None, Verbose=True):
        """
        Interpolates all variables of inFile and optionally
        save them to file *npzFile*
        """
        from gfio import GFIOctl, GFIOHandle

        # Instiate grads and open file
        # ----------------------------
        fh = GFIOctl(inFile)
        self.sample = GFIOHandle(inFile)
        if onlyVars is None:
            onlyVars = fh.vname

        nt = self.lon.shape
        nz = self.nz
        tymes = self.tyme
        lons = self.lon
        lats = self.lat

        # Loop over variables on file
        # ---------------------------
        for v in onlyVars:
            if Verbose:
                print("<> Sampling ", v)
            var = fh.sample(v,lons,lats,tymes,Verbose=Verbose)
            if len(var.shape) == 1:
                self.sample.__dict__[v] = var
            elif len(var.shape) == 2:
                var = var.T # shape should be (nobs,nz)
                self.sample.__dict__[v] = var
            else:
                raise IndexError('variable <%s> has rannk = %d'%len(var.shape))

        if npzFile is not None:
            savez(npzFile,**self.sample.__dict__)            

       
    def sampleLoadz(self,npzFile):
        """
        Loads sample from npz file.
        """
        from grads.gahandle import GaHandle
        self.sample = GaHandle(npzFile)
        npz = load(npzFile)
        for v in list(npz.keys()):
            self.sample.__dict__[v] = npz[v]
                
    def addVar(self,ga,expr='albedo',vname=None,clmYear=None,tight=True):
        """
        Given a grads object *ga* having the correct MERRA file as default,
        interpolates *var* to obs location and saves it as an attribute
        named *vname*.

        If *tight* is True, domain will be restricted conserve memory. This feature
        has proven somewhat unstable for reasons yet TBD.
        """

        U = MISSING * ones(self.nobs)
        if vname == None:
            vname = expr

        # nearest time
        # ------------
        t = _gatime(self.nymd,self.nhms)
        if clmYear != None:
            t = t[:-4] + str(clmYear) # replace year
        ga('set time '+t,Quiet=True)

        # To conserve memory, restrict domain with 1 gridpoint halo
        # ---------------------------------------------------------
        if tight:
            fh = ga.query("file")
            x1, x2  = self.lon.min(),self.lon.max()
            y1, y2  = self.lat.min(),self.lat.max()
            ga('set lon %f %f'%(x1,x2),Quiet=True)
            ga('set lat %f %f'%(y1,y2),Quiet=True)
            qh = ga.query("dims")
            x1, x2 = (qh.xi[0]-1,qh.xi[1]+1)
            y1, y2 = (max(1,qh.yi[0]-1),min(fh.ny,qh.yi[1]+1)) # in [1,ny]
            ga('set x %d %d'%(x1,x2),Quiet=True) # make sure x range is int
            ga('set y %d %d'%(y1,y2),Quiet=True) # make sure y range is int
            expr_ = ga.exp(expr)
        else:
            expr_ = ga.expr(expr)
        u, levs = ga.interp(expr_, self.lon, self.lat )
        U = u.data
        if len(shape(U)) == 0:
             U = U * ones(1) # so that we can slice it later

        self.__dict__[vname] = U

#..............................................................
if __name__ == "__main__":

     import os
    
#     inDir = '/nobackup/LIDAR/dR_Fortuna-2-4-b4'
#     inFile = inDir+'/dR_Fortuna-2-4-b4.calipso_532nm.20090715.nc'
#     aer_v = '/home/adasilva/GAAS/LIDAR/aer_Nv.ddf'
#     npzFile = os.path.basename(inFile.replace('532nm','aer').replace('.nc','.npz'))

#    f = LIDAR_L2(inFile)
#    f.sampleFile(aer_v,npzFile=npzFile)
    
     inDir = '/nobackup/2/vbuchard/LIDAR/dR_Fortuna-2-4-b4'
     inFile = inDir+'/dR_Fortuna-2-4-b4.calipso_532nm.20090715.nc'
     aer_v = '/nobackup/2/vbuchard/LIDAR/dR_Fortuna-2-4-b4/aer_v.ddf'
     npzFile = os.path.basename(inFile.replace('532nm','aer').replace('.nc','.npz'))

     f = LIDAR_L2(inFile)
     f.sampleFile(aer_v,npzFile=npzFile)
