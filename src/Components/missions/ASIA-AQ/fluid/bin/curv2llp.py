#
# Regridding curvilinear to regular lat-lon-pressure grid.
#

from netCDF4 import Dataset
from scipy import interpolate
import numpy as npy
from sys import stdout

UNDEF   = -999.0      # undefined value (outside of domain, or from file)
MISSING = -888.0     # will be filed in by interpolation

def getMinMax(Mask,I,J,cX,axis):
    rX = Mask.copy()
    rX[J,I] = cX
    rX = npy.ma.masked_array(rX,rX==MISSING)
    xmin, xmax = rX.min(axis=axis), rX.max(axis=axis)
    return (xmin,xmax)

def _varList(f):
    from grads import GrADS
    if isinstance(f,Dataset):
        return f.variables.keys()
    elif isinstance(f,GrADS):
        fh = f.query('file')
        return fh.vars
    else:
        raise NotImplementedError, 'Only NetCDF or GrADS supported for now'

def _varGet(f,name):
    from grads import GrADS, numtypes
    if isinstance(f,Dataset):
        return f.variables[name]
    elif isinstance(f,GrADS):
        fh = f.query('file')
        qh = f.query('ctlinfo')
        nz = fh.var_levs[fh.vars.index(name)]
        if nz==0:
            f('set z 1')
            shp = (fh.nt,fh.ny,fh.nx)
        else:
            shp = (fh.nt,nz,fh.ny,fh.nx)
        var = numtypes.GaField(npy.zeros(shp))
        var.name = name
        var.long_name = fh.Vars[name].long_name
        var.missing_value = qh.undef
        return var
    else:
        raise NotImplementedError, 'Only NetCDF or GrADS supported for now'

def _varGetData(f,n,name):
    from grads import GrADS
    if isinstance(f,Dataset):
        return f.variables[name][n]
    elif isinstance(f,GrADS):
        fh = f.query('file')
        qh = f.query('ctlinfo')
        f('set t %d'%(n+1))
        f('set x 1 %d'%fh.nx)
        f('set y 1 %d'%fh.ny)
        nz = fh.var_levs[fh.vars.index(name)]
        if nz==0:
            f('set z 1')
        else:
            f('set z 1 %d'%nz)
        f.flush()
        var = f.expr(name)
        var[var.mask] = UNDEF
        # print '\n ', name, ' --- ', var.shape, var.min(), var.max()
        return var
    else:
        raise NotImplementedError, 'Only NetCDF or GrADS supported for now'
    
class Curv2LLP(object):

    def __init__(self,xLon,xLat,rP=None):
        """
        This class implements a somewhat generic regridding capability from
        a generic lat-lon curvilinear grid on model vertical levels to a regular
        lat-lon-pressure grid.
        
        On input:
           xLon --- 2D array of longitudes on curvlinear grid (deg E), OR
                    1D array of longitudes on regular grid (deg E) 
           xLat --- 2D array of latitudes on curvlinear grid (deg N), OR
                    1D array of latitudes on regular grid (deg N)
           cP   --- 3D array of pressures (hPa)

        Arrays are assumed to have shape (y,x) and (z,y,x)
        
        """

        # Curvilinear grid
        # ----------------
        if len(xLon.shape) == 2:

            cLon, cLat = xLon, xLat            
            ny, nx = cLon.shape
            self.curv = True # curvilinear grid
            
            # Regular grid coordinates
            # -------------------------
            self.lon = npy.linspace(cLon.min(),cLon.max(),nx)
            self.lat = npy.linspace(cLat.min(),cLat.max(),ny)
        
            # Indices of curvilinear grid on lat-lon grid
            # -------------------------------------------
            self.I = (0.5+(cLon-self.lon[0])/(self.lon[1]-self.lon[0])).astype('int')
            self.J = (0.5+(cLat-self.lat[0])/(self.lat[1]-self.lat[0])).astype('int')
            
            # Mask region outside of domain
            # -----------------------------
            Mask = MISSING * npy.ones((ny,nx))
            xmin, xmax = getMinMax(Mask,self.I,self.J,cLon,axis=1)
            ymin, ymax = getMinMax(Mask,self.I,self.J,cLat,axis=0)
            for j in range(ny):
                M = (self.lon[:]<xmin[j])|(self.lon[:]>xmax[j])
                Mask[j,M] = UNDEF
            for i in range(nx):
                M = (self.lat[:]<ymin[i])|(self.lat[:]>ymax[i])
                Mask[M,i] = UNDEF
            self.Mask = Mask

        # Regular lat/lon grid
        # --------------------
        elif len(xLon.shape) == 1:
            ny, nx = len(xLat), len(xLon)
            self.lon, self.lat = xLon, xLat
            self.curv = False # Not a curvilinear grid

        else:
            raise ValueError, 'Invalid dimensionality of xLon: %d'%len(xLon.shape)

        # Horizontal grid: destination size same as source grid
        # -----------------------------------------------------
        self.ny, self.nx = ny, nx
            
        # Vertical aspects
        # ----------------
        self.k1, self.alpha = [], []
        self.rP = rP
        if rP is not None:
            self.np = len(rP)
        else:
            self.np = 0

    def get_cP(self,n,ncin):
        """
        Given time level "n", return 3D pressure for this time.
        This is an abstract method, it needs to be re-defined
        when subclassing this class.
        """
        raise NotImplementedError, 'Method get_cP not implemented in Curv2LLP' 

    def set_pInterp(self,n,ncin):
        """
        Retrieve vertical pressure and compute
        vertical interpolation weights for time step "n".
        """
        cP = self.get_cP(n,ncin)
        self.k1, self.alpha = [], []
        if self.rP is not None:
            bot_top = (cP[0,0,0]>cP[1,0,0]) # bottom-top or top-bottom?
            self.nz = cP.shape[0]
            for p in self.rP:
                k1 = UNDEF*npy.ones((self.ny,self.nx)).astype('int')
                a = UNDEF*npy.ones((self.ny,self.nx))
                if bot_top:
                    self.kdel = 1
                    for k in range(self.nz-1):
                        k1 = npy.where( (cP[k+1]<=p)&(p<cP[k,:,:]), k, k1 )
                        a  = npy.where(k==k1, (p-cP[k])/(cP[k+1]-cP[k]), a)
                else:
                    self.kdel = -1
                    for k in range(1,self.nz):
                        k1 = npy.where( (cP[k-1]<=p)&(p<cP[k,:,:]), k, k1 )
                        a  = npy.where(k==k1, (p-cP[k])/(cP[k-1]-cP[k]), a)
                self.k1 += [k1,]
                self.alpha += [ a, ]

    def regrid(self,cV):
        """
        Regrid a 2D or 3D variable from a curvlinear to a regular
        lat-lon[-pressure] grid. 
        """
        if len(cV.shape)==2:
            return self.regrid2d(cV[:,:])
        elif len(cV.shape)==3:
            return self.regrid3d(cV[:,:,:])
        else:
            raise ValueError, 'Invalid dimmensionality of field'

    def regrid2d(self,cV):
        """
        Regrid a 2D variable from a curvlinear to a regular lat-lon grid.
        """

        # not a curvilinear grid, nothing to do
        # -------------------------------------
        if not self.curv:
            return cV      
        
        ny, nx = cV.shape

        # Put Var on regular lat/lon grid (with holes)
        # --------------------------------------------
        rV = self.Mask.copy()
        rV[self.J,self.I] = cV

        # Fill missing values in middle of domain
        # ---------------------------------------
        for j in range(ny):  # 1D interpolation in lon
            z = rV[j,:]
            B = (z==MISSING)
            if any(B):
                G = (z!=UNDEF)&(z!=MISSING)
                if sum(G)>=2:
                    f = interpolate.interp1d(self.lon[G],z[G],kind='linear',
                                         bounds_error=False,fill_value=UNDEF)
                    rV[j,B] = f(self.lon[B])

        return rV

    def regrid3d(self,cV):
        """
        Regrid a 3D variable from a curvlinear to a regular lat-lon grid.
        """

        rV = UNDEF * npy.ones((self.np,self.ny,self.nx))
        
        # Interpolate to every desired vertical pressure level
        # ----------------------------------------------------
        for k in range(self.np):
            k1, a = self.k1[k], self.alpha[k]
            I = k1!=UNDEF
            kmin, kmax = int(k1[I].min()), int(k1[I].max())
            cU = UNDEF * npy.ones((self.ny,self.nx))
            for k_ in range(kmin,kmax+1):
                cU = npy.where(k1==k_,
                               (1-a)*cV[k_,:,:]+a*cV[k_+self.kdel,:,:], cU)
            rV[k,:,:] = self.regrid(cU)

        return rV
        
    def writeNC ( self, tyme, options, ncin, zlib=False ):
        """
        Write a NetCDF file with sampled GEOS-5 variables at the station locations
        described by (lon,lat).
        
        On input:

           tyme  ---  array of time levels to write to file  (usually 1)
           options -  Command line options:
                      outFile   ---  output filenane
                      format    ---  file format
                      dryrun    ---  writes out file with zeros
                      verbose
        
        """

        # Open output NC file
        # -------------------
        nc = Dataset(options.outFile,'w',format=options.format)

        # Set global attributes, borrow from input file
        # ---------------------------------------------
        if isinstance(ncin,Dataset):
            for att in ncin.__dict__.keys():
                nc.__setattr__(att,ncin.__dict__[att]) 
        nc.History = 'Converted to pressure-lat-lon using curv2llp'
        nc.Conventions = 'CF'

        # Create dimensions
        # -----------------
        nx = nc.createDimension('lon', self.nx )
        ny = nc.createDimension('lat', self.ny )
        if self.np>0:
            np = nc.createDimension('lev', self.np )
        nt = nc.createDimension('time', len(tyme)) 
     
        # Coordinate variables
        # --------------------
        lon = nc.createVariable('lon','f4',('lon',),zlib=zlib)
        lon.long_name = 'Longitude'
        lon.units = 'degrees_east'
        lon[:] = self.lon[:]
    
        lat = nc.createVariable('lat','f4',('lat',),zlib=zlib)
        lat.long_name = 'Latitude'
        lat.units = 'degrees_north'
        lat[:] = self.lat[:]

        if self.np > 0:
            lev = nc.createVariable('lev','f4',('lev',),zlib=zlib)
            lev.long_name = 'Pressure'
            lev.units = 'hPa'
            lev.positive = 'down'
            lev.axis = 'z'
            lev[:] = self.rP

        time = nc.createVariable('time','i4',('time',),zlib=zlib)
        time.long_name = 'Time'
        t0 = tyme[0]
        time.units = 'seconds since %s'%t0.isoformat(' ')
        time[:] = npy.array([(t-t0).total_seconds() for t in tyme])

        # Create variables on NetCDF file.
        # -------------------------------
        ncVar = dict()
        # Vars = options.Vars or ncin.variables.keys()
        Vars = options.Vars or _varList(ncin)
        for v in Vars:
            if v in options.ignore: continue
            if options.verbose: stdout.write("[] Creating <%s> ... "%v)
            # var = ncin.variables[v]
            var = _varGet(ncin,v)
            if len(var.shape) == 3: # 2D variable
                dim = ('time','lat','lon')
                if options.verbose: stdout.write("2d ... ")
            elif len(var.shape) == 4: # 3D variable
                if self.np == 0:
                    if options.verbose: stdout.write("ignored\n")
                    continue
                dim = ('time', 'lev','lat','lon')
                if options.verbose: stdout.write("3d ... ")
            else:
                if options.verbose:
                    stdout.write("ignored\n")
                continue
            this = nc.createVariable(var.name,'f4',dim,zlib=zlib)
            try:
                this.long_name = var.__dict__[u'long_name']
            except:
                this.long_name = v
            this.missing_value = UNDEF # out of domain
            try:
                this.units = var.__dict__[u'units']
            except:
                this.units = 'unknown'
            ncVar[v] = this
            if options.verbose: stdout.write("ok\n")

        # Write variables for each time
        # -----------------------------
        for n in range(len(tyme)):
            if options.verbose:
                stdout.write("[] Writing on %s ... "%tyme[0].isoformat())
            self.set_pInterp(n,ncin) # set interpolation weights
            for v in Vars:
                if v in options.ignore: continue
                var = _varGetData(ncin,n,v)
                if len(var.shape) == 2: # 2D variable
                    shp = ( self.ny, self.nx )
                elif len(var.shape) == 3: # 3D variable
                    if self.np == 0: continue
                    shp = ( self.np, self.ny, self.nx )
                else: continue
                this = ncVar[v]
                if options .dryrun:
                    this[n,:] = npy.zeros(shp)[:]
                else:
                    this[n,:] = self.regrid(var)[:]
                if options.verbose: stdout.write("%s ... "%v)
            if options.verbose: stdout.write("done\n")
            
        # Close the file
        # --------------
        nc.close()

        if options.verbose:
            print("<> wrote %s file %s"%(options.format,options.outFile))
  
