#--------------------------------------------------------------------------
#
#    Copyright (C) 2006-2008 by Arlindo da Silva <dasilva@opengrads.org>
#    All Rights Reserved.
#
#    This program is free software; you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation# using version 2 of the License.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY# without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with this program# if not, please consult  
#              
#              http://www.gnu.org/licenses/licenses.html
#
#    or write to the Free Software Foundation, Inc., 59 Temple Place,
#    Suite 330, Boston, MA 02111-1307 USA
#
#------------------------------------------------------------------------

""" 
This module extends the GrADS client class by providing methods for
exchanging n-dimensional NumPy array data between Python and
GrADS.
"""

__version__ = '1.1.3'


from gacore       import *
from numtypes     import *
from simplekml    import SimpleKML

from numpy        import zeros, ones, average, newaxis, sqrt, pi, cos, inner, \
                         arange, fromfile, float32, ma, reshape, ndarray, \
                         abs, size, meshgrid, shape, tile

from numpy.linalg import svd, lstsq

class GaNum(GaCore):
    """
    This class extends the GrADS client class by providing methods
    for data exchange and numerical computations:

    _Methods provided:
       exp  -  exports a GrADS expression into a NumPy array, with metada
       imp  -  imports NumPy array (+metadata) into GrADS
       eof  -  compute Empirical Orthogonal Functions (EOFS) from expressions 
       lsq  -  least square parameter estimation from expressions 

    """

#........................................................................

    def exp (self, expr):
        """
        Exports GrADS expression *expr*, returning a GrADS Field.

            F = self.exp(expr)

        where

            F  ---  GrADS field
            
        Generalized Expressions
        =======================

        For convenience, *expr* can also be a GrADS Field.  In such
        cases, the input Field is just returned back. This
        *overloading* is useful for writing high level functions that
        work equaly on a GrADS expression to be exported or on GrADS
        fields already known to python.

        Limitation
        ==========

        This function does not handle varying ensemble dimensions in
        GrADS v2.

        """

#       If IPC extension is not available, then try expr() instead
#       ----------------------------------------------------------
        if not self.HAS_IPC:
            return self.expr(expr)

#       For convenience, allows calls where expr is not a string, in which
#        case it returns back the input field or raise an exception
#       -------------------------------------------------------------------
        if type(expr) in StringTypes:
            pass # OK, will proceed to export it from GrADS
        elif isinstance(expr,GaField):
            return expr # just return input
        elif isinstance(expr,ndarray):
            return expr # this is handy for 'lsq'
        else:
            raise GrADSError, "input <expr> has invalid type: %s"%type(expr)

#       Retrieve dimension environment
#       ------------------------------
        dh = self.query("dims", Quiet=True) 
        t1, t2 = dh.t
        z1, z2 = dh.z 
        nx, ny, nz, nt = (dh.nx, dh.ny, dh.nz, dh.nt)
   

#       Shortcut for 2D slices (any 2 dimensions)
#       -----------------------------------------
        if dh.rank ==2:
            return self._exp2d(expr)

#       Initial implementation: require x,y to vary for rank>2
#       Note: remove this restriction is not very hard, but requires
#             special handling the different dimension permutations separately
#             given the way GrADS invokes functions for XZ, YZ, ZT, etc
#       ----------------------------------------------------------------------
        if nx==1: raise GrADSError, 'lon must be varying but got nx=1'
        if ny==1: raise GrADSError, 'lat must be varying but got ny=1'

#       Loop over time/z, get a GrADS 2D slice at a time/z
#       --------------------------------------------------
        l = rc = 0  
        Data = None # defer allocations until we know the size
        grid = GaGrid(expr)
        grid.meta = zeros((nt,nz,20),dtype=float32)
        grid.denv = dh 
        grid.time = [] 
        grid.lev = zeros(nz,dtype=float32)         
        try:
            for t in range(t1,t2+1):
                self.cmd("set t %d"%t,Quiet=True) 
                self.cmd("q time",Quiet=True)
                grid.time.append(self.rword(1,3))
                k = 0
                for z in range(z1,z2+1):
                    self.cmd("set z %d"%z,Quiet=True)
                    field = self._exp2d(expr)
                    if Data is None:
                        ny_, nx_ = field.shape # may differ from dh.nx/dh.ny
                        Data = zeros(shape=(nt,nz,ny_,nx_), dtype=float32)
                    Data[l,k,:,:] = field.data
                    grid.lev[k] = field.grid.lev[0]
                    grid.meta[l,k,:] = field.grid.meta
                    k = k + 1
                l = l + 1

#           Record lat/lon
#           --------------
            grid.lat = field.grid.lat
            grid.lon = field.grid.lon
            amiss = grid.meta[0]

#           Remove dimensions with size 1
#           -----------------------------
            if nz==1: 
                Data = Data.reshape(nt,ny_,nx_)
                grid.dims = [ 'time', 'lat', 'lon' ]
                grid.meta = grid.meta.reshape(nt,20)
            elif nt==1: 
                Data = Data.reshape(nz,ny_,nx_)
                grid.dims = [ 'lev', 'lat', 'lon' ]
                grid.meta = grid.meta.reshape(nz,20)
            else:
                grid.dims = [ 'time', 'lev', 'lat', 'lon' ]

        except:
            self.setdim(dh)
            raise GrADSError, 'could not export <%s>'%expr


        grid.tyme = array([gat2dt(t) for t in grid.time])

#       Restore dimension environment
#       -----------------------------
        self.setdim(dh)
        return GaField(Data, name=expr, grid=grid, 
                       mask=(Data==amiss), dtype=float32)
    
#........................................................................

    def _exp2d ( self, expr, dh=None ):
        """ 
        Exports GrADS expression *expr* as a GrADS Field.
        The stdio pipes are used for data exchange.
        This is an internal version handling 2D xy slices.
        In here, *expr* must be a string.
        """

        if dh==None:
            dh = self.query("dims",Quiet=True)

#       Check environmnet
#       -----------------
        nx, ny = (dh.nx, dh.ny)
        if dh.rank !=2:
            raise GrADSError, 'expecting rank=2 but got rank=%d'%dh.rank

#       Create output handle, fill in some metadata
#       -------------------------------------------
        grid = GaGrid(expr)
        grid.denv = dh

#       Issue GrADS command, will check rc later
#       -----------------------------------------
        if self.Version[1] is '1':
            cmd = 'ipc_define void = ipc_save('+expr+',-)\n'
        else:
            cmd = 'define void = ipc_save('+expr+',-)\n'

        self.Writer.write(cmd)

#       Position stream pointer after <EXP> marker
#       ------------------------------------------
        got = ''
        while got[:5] != '<EXP>' :
            got = self.Reader.readline()

#       Read header
#       -----------
        grid.meta = fromfile(self.Reader,count=20,dtype=float32)

        amiss = grid.meta[0]
        id = int(grid.meta[1])
        jd = int(grid.meta[2])
        nx_ = int(grid.meta[3])
        ny_ = int(grid.meta[4])

#        if id!=0 or jd!=1:
        if id<0 or id>3 or jd<0 or jd>3 or id==jd:
            self.flush()
            raise GrADSError, \
                  'invalid exchange metadata (idim,jdim)=(%d,%d) - make sure <%s> is valid and that lon/lat is varying.'%(id,jd,expr)

#       Read data and coords
#       --------------------
        try:
            array_ = fromfile(self.Reader,count=nx_*ny_,dtype=float32)
            grid.lon = fromfile(self.Reader,count=nx_,dtype=float32)
            grid.lat = fromfile(self.Reader,count=ny_,dtype=float32)
        except:
            self.flush()
            raise GrADSError, 'problems exporting <'+expr+'>, fromfile() failed'

#       Annotate grid - assumes lon, lat
#       --------------------------------
        dims = ( 'lon', 'lat', 'lev', 'time' )
        grid.dims = [dims[jd],dims[id]]
        grid.time = [ dh.time[0] ]
        grid.lev = ones(1,dtype=float32) * float(dh.lev[0])

#       Check rc from asynchronous ipc_save
#       -----------------------------------
        rc = self._parseReader(Quiet=True)
        if rc:
            self.flush()
            raise GrADSError, 'problems exporting <'+expr+'>, ipc_save() failed'

        grid.tyme = array([gat2dt(t) for t in grid.time])

#       Create the GaField object
#       -------------------------
        data = array_.reshape(ny_,nx_)
        return GaField(data, name=expr, grid=grid, mask=(data==amiss) )
    
#........................................................................

    def imp ( self, name, Field ):
        """
        Sends a GrADS Field containing a NumPy array and associated 
        grid information to GrADS, defining it in GrADS as *name*.
        Notice that *Field* can be an instance of the GaField
        class or a tuple with the (Array,Grid) components.

        Limitation
        ==========

        This function does not handle varying ensemble dimensions in
        GrADS v2.

        """

#       If IPC extension is not available, barf
#       ---------------------------------------
        if not self.HAS_IPC:
            raise GrADSError, "IPC extension not available - cannot import!"

#       Resolve Field
#       -------------
        if isinstance(Field,GaField):
            grid = Field.grid
        else:
            raise GrADSError, "Field has invalid type"
                
#       Retrieve dimension environment
#       ------------------------------
        dh = self.query("dims", Quiet=True) 
        t1, t2 = dh.t
        z1, z2 = dh.z 
        nx, ny, nz, nt = (dh.nx, dh.ny, dh.nz, dh.nt)
        nxy = nx * ny

#       Initial implementation: require x,y to vary
#       Note: remove this restriction is not very hard, but requires
#             special handling the different dimension permutations separately
#             given the way GrADS invokes functions for XZ, YZ, ZT, etc
#       ----------------------------------------------------------------------
        if nx==1: raise GrADSError, 'lon must be varying but got nx=1'
        if ny==1: raise GrADSError, 'lat must be varying but got ny=1'

#       Determine actual load command
#       -----------------------------
        if name == '<display>':
            cmd = 'display ipc_load()\n'
            if nz>1 and nt>1:
                raise GrADSError, \
                      'for <display> only one of z/t can vary'+\
                      ' but got (nz,nt)=(%d,%d)'%(nz,nt) 
        else:
            if self.Version[1] is '1':
                cmd = 'ipc_define %s = ipc_load()\n'%name
            else:
                cmd = 'define %s = ipc_load()\n'%name

#       Tell GrADS to start looking for data in transfer stream
#       -------------------------------------------------------
        try:
            self.cmd("ipc_open - r")
        except GrADSError:
            raise GrADSError, '<ipc_open - r> failed; is IPC installad?'
        self.Writer.write(cmd) # asynchronous transfer

#       Reshape and get original t/z offset
#       -----------------------------------
        t1_, z1_ = (grid.denv.t[0], grid.denv.z[0])
        nt_, nz_ = (grid.denv.nt,grid.denv.nz)
        nx_ = len(grid.lon)
        ny_ = len(grid.lat)
        nxy_ = nx_ * ny_
        data = Field.data.reshape(nt_,nz_,ny_,nx_)
        meta = grid.meta.reshape(nt_,nz_,20)

#       Write the data to transfer stream
#       ----------------------------------
        try:
            for t in range(t1,t2+1):
                l = t - t1_
                for z in range(z1,z2+1):
                    k = z - z1_
                    mx = int(meta[l,k,3])
                    my = int(meta[l,k,4])
                    if mx!=nx_ or my!=ny_:
                        self.flush()
                        raise GrADSError, \
                             'nx/ny mismatch; got (%d,%d), expected (%d,%d)'%\
                             (mx,my,nx_,ny_)
                    meta[l,k,:].tofile(self.Writer)
                    data[l,k,:,:].tofile(self.Writer)
                    grid.lon.tofile(self.Writer)
                    grid.lat.tofile(self.Writer)
                    self.Writer.flush()
        except:
            self.flush()
            self.setdim(dh)
            raise GrADSError, \
                  'could not import <%s>, tofile() may have failed'%name


#       Check rc from asynchronous ipc_save
#       -----------------------------------
        rc = self._parseReader(Quiet=True)
        self.flush()

#       Restore dimension environment
#       -----------------------------
        self.setdim(dh)
        self.cmd("ipc_close")
        if rc:
            raise GrADSError, 'problems importing <'+name+'>, ipc_load() failed'

#........................................................................

    def expr (self, expr):
        """
        Evaluates a GrADS expression returning a GrADS Field. This is similar
        to the exp() method except that the resulting GaField cannot be
        imported back into GrADS. It relies on *gacore* methods eval()
        and coords() to retrieve the data and coordinate information. 
        """

#       For convenience, allows calls where expr is not a string, in which
#        case it returns back the input field or raise an exception
#       -------------------------------------------------------------------
        if type(expr) in StringTypes:
            pass # OK, will proceed to retrieve it from GrADS
        elif isinstance(expr,GaField):
            return expr # just return input
        elif isinstance(expr,ndarray):
            return expr # this is handy for 'lsq'
        else:
            raise GrADSError, "input <expr> has invalid type: %s"%type(expr)

        d = self.eval(expr)
        c = self.coords()
        g = GaGrid(expr,coords=c)

        Data = reshape(d,c.shape)
        F = GaField(Data,mask=(Data==c.undef),name=expr,grid=g)

        return F

#........................................................................

    def eof ( self, expr, transf='anomaly', metric='area', keep=None):
        """ 
        Given a GrADS generalized expression *expr*, calculates Empirical 
        Orthogonal Functions (EOFS) using Singular Value Decomposition (SVD). 
        
            V, d, c = self.eof(expr)
        
        where

            V  ---  A specialized GrADS Field holding eigenvectors
                    with *add* offset and *scale* factors to aid
                    subsequent decomposition in terms of V
            d  ---  NumPy array with eigenvalues
            c  ---  NumPy array with principal components

        The other optional parameters are:

        transf     
            Type of pre-processing transform to be applied:
            None     ---  time series as is
            anomaly  ---  remove time mean
            z-score  ---  remove time mean and divide 
                          by standard deviation
        metric    
            Determines whether to scale the timeseries prior
            to calculation of the EOFs; this is equivalent to 
            choosing a particular norm. Acceptable values are:
            None   ---  do not scale
            'area' ---  multiply by cos(lat), the default

        keep
            How many eigenvectors to keep:
            None  ---  in this case keep as many vectors as 
                       there are timesteps (nt) 
            n     ---  keep "n" eigenvectors

        Notice that *expr* on input is a *generalized expression* in the
        sense that it can contain a string with a GrADS expression to be
        evaluated or a valid GrADS field. See method *exp* for additional
        information.

        IMPORTANT: the input (masked) array mask (undef values) must be
                   the same for all times.
        
        """

#       At least 2 time steps
#       ---------------------
        dh = self.query("dims",Quiet=True)
        if dh.nt < 2:
            raise GrADSError, \
                  'need at least 2 time steps for EOFS but got nt=%d'%dh.nt
        nt, nz, ny, nx = (dh.nt, dh.nz, dh.ny, dh.nx)

#       Export N-dimensional array
#       --------------------------
        u = self.exp(expr)
        g = u.grid

#       Reshape as 4D
#       -------------
        nx = len(g.lon) # may differ from dh.nx
        ny = len(g.lat) # may differ from dh.ny
        u = u.reshape(nt,nz,ny,nx)

#       Remove time mean if necessary
#       -----------------------------
        offset = ma.zeros((nz,ny,nx),dtype=float32) # place holder
        if transf==None:
            pass
        elif transf=='anomaly' or transf=='z-score':
            offset = average(u,axis=0)
            u = u - offset[newaxis,:,:,:]
        else:
            raise GrADSError, 'Unknown transf <%s>'%transf
    
#       Scale by stdv if z-scores required
#       ----------------------------------
        scale = ma.ones((nz,ny,nx),dtype=float32) # place holder
        if transf=='z-score':
            scale = sqrt(average(u*u,axis=0))
            u = u / scale[newaxis,:,:,:]

#       Apply metric if so desired
#       Note: may need delp for cases with nz>1
#       ---------------------------------------
        if metric=='area':
            factor = sqrt(cos(pi*g.lat/180.))
            u = u * factor[newaxis,newaxis,:,newaxis]
            scale = scale * factor[newaxis,newaxis,:,newaxis]

#       Singular value decomposition, reuse u
#       -------------------------------------
        I = u.mask[0,:,:,:]==False  # un-masked values
        fill_value = u.fill_value   # save it for later
        pc, d, u = svd(u.data[:,I],full_matrices=0)

#       Trim eigenvectors
#       -----------------
        if keep==None:
            nv = nt
        else:
            nv = min(keep,nt)
            u = u[0:nv,:]
            d = d[0:nv]
            pc = pc[:,0:nv]

#       Adjust grid properties
#       ----------------------
        g.dims[0] = 'eof'
        g.time = arange(nv)
        g.eof = arange(nv)

#       Eigenvalues/coefficients
#       ------------------------
        d = d * d / (nt - 1)
        pc = (nt -1) * pc.transpose()

#       Normalize eigenvectors
#       ----------------------
        for i in range(nv):
            vnorm = _norm(u[i,:])
            u[i,:] = u[i,:] / vnorm
            pc[i,:] = pc[i,:] * vnorm
    
        # Scatter eigenvectors
        # --------------------
        g.meta = g.meta[0:nv] 
        u = _scatter(u,I,(nv,nz,ny,nx),fill_value)

#       Let's make sure "u" is a bonafide GaGield
#       -----------------------------------------
        u = GaField(u.data, name=expr, grid=g, mask=u.mask)
        u.offset = offset.squeeze()
        u.scale = scale.squeeze()

#       Note: since GrADS v1 does not know about e-dimensions yet,
#       we let it think that the EOF dimension is the time dimension

#       All done
#       --------
        return (u, d, pc)

#.....................................................................

    def lsq (self, y_expr, x_exprs, Bias=False, Mask=None):
        """
        Given a target GrADS expression *y_expr* and a tuple of predictor
        GrADS expressions *x_exprs*, returns a NumPy array with linear 
        regression coefficients 

            c, info = self.lsq(y_expr, x_exprs)

        where *info* contains information about the minimization:

            info.residuals  ---  sum square of residuals
            info.rank       ---  rank of the predictor matrix
            info.s          ---  singular values of predictor matrix

        When Bias=False (default) the residual

            y - c[0] * x[:,0] + c[1] * x[:,1] + ... + c[n] * x[:,N-1]

        is minimized in the last square sense, where *x* and *y* are
        NumPy arrays associated with the following GrADS Fields:

            Y = self.exp(y_expr)
            X[:,n] = self.exp(x_exprs[n]),  n = 0, ..., N-1

       When Bias=True, the additional predictor array 

            x[N] = ones(y)

       is included in the calculation, resulting in an output array of
       size N+1. This is equivalent to allowing for an *intercept* in
       the regression equation.

       The optional Mask is a 1D logical array with the location of the
       data to be included (see compress() method).

       On input, all expressions are *generalized expressions* in the
       sense that they can contain a string with a GrADS expression to be
       evaluated or a valid GrADS field. See method *expr* for additional
       information.

       """

        N = len(x_exprs)
        if N<1:
            raise GrADSError, \
                'expecting at least one predictor but got %d'%N
        if Bias: N = N + 1
        
#       Retrieve target
#       ---------------
        f = self.exp(y_expr)
        y = f.ravel()
        if Mask!=None:  
            y = y.compress(Mask)
        M = y.size

#       Retrieve predictors
#       -------------------
        X = ones((M,N),dtype=float32)
        for n in range(len(x_exprs)):
            f = self.exp(x_exprs[n])
            x = f.ravel()
            if Mask!=None: 
                x = x.compress(Mask)
            X[:,n] = x

#       Perform LS minimization
#       -----------------------
        info = GaHandle('lsq')
        (c, info.residuals, info.rank, info.s) = lstsq(X,y)

#       All done
#       --------
        return (c, info)

#   ..................................................................

    def _interpXY ( self, expr, lons, lats, levs=None, dh=None, **kwopts):
        """
        Evaluates GrADS expression (or GaField) and interpolates it to
        the the (longitude,latitude) points given the input arrays
        (lons,lats) on input. Both x/y dimensions must be
        varying. When the z-dimenson is varying as well a curtain
        slice is returned. For now, the time dimension must be fixed.
        Example:

          tau, levs = ga.interp('duexttau',lons,lats)

        where *levs* is an 1D array with the versical levels. The optional
        **kwopts arguments are passwd to the interpolate() function.

        Note: the basemap interpolation routine requires longitudes in
        the range [-180,180]. When *expr* is a string the longitudes are
        set to the correct range. However, when *expr* is a GaField
        the user must make sure this is the case or undefs may result.
        """

#       Check dim environment
#       ---------------------
        if dh==None:
            dh = self.query("dims", Quiet=True)
        if dh.nx==1 or dh.ny==1:
            raise GrADSError, \
            "expecting varying x/y dimensions but got (nx,ny) = (%d,%d)"\
            %(dh.nx,dh.ny)
        if dh.nt>1:
            raise GrADSError, \
            "sorry, cannot interpolate with varying time dimension"

#       Evaluate GrADS expression
#       -------------------------
        if dh.lon[0]>180. or dh.lon[1]>180:
            self.cmd('set lon -180 180',Quiet=True) # assume global grid
        Z = self.exp(expr)
        g = Z.grid
        self.cmd('set x %s %s'%dh.x,Quiet=True)

#       Loop over vertical levels
#       -------------------------
        n = size(lons)
        lon_, lat_ = reshape(lons,(n,1)), reshape(lats,(n,1))
        if len(Z.shape)==2:
            y = interpolate(Z, g.lon, g.lat, lon_, lat_, 
                       masked=True, **kwopts)
        else: 
            y = ma.masked_array(zeros((n,1,dh.nz),dtype=float32))
            for z in range(dh.nz): # because Interp needs 2D arrays on input
                y[:,:,z] = interpolate(Z[z], g.lon, g.lat, lon_, lat_,
                                       masked=True, **kwopts)
                
#       Return array with same shape as the input lons/lats, 
#        with possibly an additional dimension in case the
#        z-dimension is varying
#       ----------------------------------------------------
        y = ma.masked_array(y,dtype=float32) # for consistency of types
        return (y.squeeze(),g.lev)

#   ..................................................................

    def sampleXY ( self, expr, lons, lats, levs=None, dh=None, **kwopts):
        """
        Evaluates GrADS expression (or GaField) and interpolates it to
        the the (longitude,latitude) points given the input arrays
        (lons,lats) on input. Both x/y dimensions must be
        varying. When the z-dimenson is varying as well a curtain
        slice is returned. The output is a special case of a GaField,
        where the first axis contains the "observational dimension". 
        
        Example:

          tau = ga.sampleXY('duexttau',lons,lats)

        The trajectory coordinates are (lons,lats) returned in the
        grid atrributes, e.g., tau.grid.lon, tau.grid.lat.

        The optional **kwopts arguments are passed to the
        interpolate() function.

        Note: the basemap interpolation routine requires longitudes in
        the range [-180,180]. When *expr* is a string the longitudes are
        set to the correct range. However, when *expr* is a GaField
        the user must make sure this is the case or undefs may result.
        """

        # Inputs must be 1D arrays
        # ------------------------
        if len(lons.shape)!=1 or len(lats.shape)!=1:
            raise GrADSError, "lons, lats, time must be 1D arrays"
        
        
        # Retrieve dimension environment
        # ------------------------------
        dh = self.query("dims", Quiet=True) 

        # Loop over time, performing interpolation
        # ----------------------------------------
        g = GaGrid("sampleXY")
        g.time = []
        V  = ma.masked_array(zeros((len(lons),dh.nt,dh.nz)),dtype=float32)
        for t in dh.ti:
            n = t - dh.ti[0]
            self.cmd('set t %d'%t, Quiet=True)
            v, g.lev = self._interpXY ( expr, lons, lats, 
                                           levs=levs, 
                                           **kwopts)
            if len(v.shape)==1:
                V[:,n,0] = v
            else:
                V[:,n,:] = v
            qh =self.query("time",Quiet=True)
            g.time.append(qh.t1)

        g.dims = ['obs',]
        if dh.nt>1:
            g.dims.append('time')
        if dh.nz>1:
            g.dims.append('lev')
        
        g.lon, g.lat = (lons, lats) # "obs" coordinates
        g.tyme = array([gat2dt(t) for t in g.time])

#       Restore dimension environment
#       -----------------------------
        self.setdim(dh)
        V = V.squeeze()
        return GaField(V.data, name=expr, grid=g, 
                       mask=V.mask, dtype=float32)
#   ..................................................................

    def sampleXYT ( self, expr, lons, lats, tyme=None,
                    levs=None, dh=None, Verbose=False, **kwopts):
        """
        Evaluates GrADS expression (or GaField) and interpolates it to
        the the (longitude,latitude,time) points given the input arrays
        (lons,lats,tyme) on input. Both x/y dimensions must be
        varying. When the z-dimenson is varying as well a curtain
        slice is returned. If *tyme* is not specified it reverts to
        method sampleXY(). The output is a special case of a GaField,
        where the first axis contains the "observational dimension". 
        
        Example:

          tau = ga.sampleXYT('duexttau',lons,lats,

        The trajectory coordinates (lons,lats,tyme) are returned in
        the grid atrributes, e.g., tau.grid.lon, tau.grid.lat,
        tau.grid.tyme.
    
        The optional **kwopts arguments are passed to the
        interpolate() function.  Notice that *tyme* is an array of
        datetime objects.

        Note: the basemap interpolation routine requires longitudes in
        the range [-180,180]. When *expr* is a string the longitudes are
        set to the correct range. However, when *expr* is a GaField
        the user must make sure this is the case or undefs may result.
        """

        # Revert back to InterpXY if no time is specified
        # -----------------------------------------------
        if tyme is None:
            return self.sampleXY ( expr, lons, lats, 
                                   levs=levs, dh=dh, 
                                   **kwopts)

        # Inputs must be 1D arrays
        # ------------------------
        if len(lons.shape)!=1 or len(lats.shape)!=1 or len(tyme.shape)!=1:
            raise GrADSError, "lons, lats, tyme must be 1D arrays"
        
        # Retrieve dimension environment
        # ------------------------------
        dh = self.query("dims", Quiet=True) 

        # Find GrADS times bracketing the input time array
        # ------------------------------------------------
        self.cmd('set time %s'%dt2gat(tyme[0]),Quiet=True)
        qh = self.query("dims",Quiet=True)
        tbeg = int(qh.t[0])
        if tyme[0] < gat2dt(qh.time[0]):
            tbeg = tbeg - 1
        self.cmd('set time %s'%dt2gat(tyme[-1]),Quiet=True)
        qh = self.query("dims",Quiet=True)
        tend = int(qh.t[0])
        if tyme[-1] > gat2dt(qh.time[0]):
            tend = tend + 1

        # Check if (tbeg,tend) is in range of default file
        # ------------------------------------------------
        fh = self.query("file",Quiet=True)
        if tbeg<1 or tbeg>fh.nt:
            raise GrADSError, "(tbeg,tend) outside of range (1,%d)"%fh.nt

        # Find time step
        # --------------
        dt = self._getDatetime(tbeg+1) - self._getDatetime(tbeg)
        dt_secs = dt.total_seconds()
        
        # Loop over time, producing XY interpolation at each time
        # -------------------------------------------------------
        V, I = [], []
        for t in range(tbeg,tend+1):
            now = self._getDatetime(t) # grads time is set to t
            if Verbose: print " [] XY Interpolating at ", now
            i = (tyme>=now-dt) & (tyme<=now+dt)
            if any(i):
                self._tightDomain(lons[i],lats[i]) # minimize I/O
                v, levs = self._interpXY(expr, lons[i], lats[i], levs=levs, **kwopts)
            else:
                v = None
            V.append(v)
            I.append(i)
            
        # Now perform the time interpolation
        # ----------------------------------
        N = len(lons) 
        if len(levs)>1:
            shp = [N,len(levs)]
        else:
            shp = [N,]
        v  = ma.masked_array(zeros(shp),dtype=float32)
        v1, v2 = v.copy(), v.copy() # scratch space
        n = 0
        for t in range(tbeg,tend):
            now = self._getDatetime(t) 
            v1[I[n]], v2[I[n+1]] = V[n], V[n+1]
            j = (tyme>=now) & (tyme<=now+dt)
            if any(j): 
                a = array([r.total_seconds()/dt_secs for r in tyme[j]-now],dtype=float32) 
                if len(shp)==2: # has vertical levels
                    a = tile(a,(shp[1],1)).T # replicate array
                v[j] = (1-a) * v1[j] + a * v2[j]
            n += 1

        # Grid
        # ----
        g = GaGrid("sampleXYT")
        g.lev = levs
        g.lon, g.lat, g.tyme = (lons, lats, tyme)
        g.time = array([dt2gat(t) for t in g.tyme])
        g.dims = ['obs',]
        if dh.nz>1:
            g.dims.append('lev')
        
        # Restore dimension environment
        # -----------------------------
        self.setdim(dh)
        
        return GaField(v.data, name=expr, grid=g, 
                       mask=v.mask, dtype=float32)
    
#   ..................................................................

    def sampleKML ( self, expr, kml_filename,
                    speed=90.,t0=None,metric=True,noSinglePoint=True,
                    **kwopts):
        """
        Evaluates GrADS expression (or GaField) and interpolates it to
        the the (longitude,latitude) points given in the input KML
        file (from Google Maps).

        On input, *speed* is the average speed, in km/h if *metric* is
        True, or miles/h otherwise. The initial time *t0* (a datetime
        object in UTC, not local time) defaults to now if not specified.
        
        Both x/y dimensions must be varying. When the z-dimenson is
        varying as well a curtain slice is returned.
        
        The output is a special case of a GaField, where the first
        axis contains the"observational dimension".

        Example:

          tau = ga.sampleKML('duexttau','directions_to_acadia.kml')

        The route coordinates are returned in the grid atrributes,
        e.g., tau.grid.lon, tau.grid.lat, tau.grid.tyme, tau.grid.dst,
        where *dst* is the distance from the first point in the route.
        
         The optional **kwopts arguments are passed to the sampleXYT()
        method.

        """

        kml = SimpleKML(kml_filename)
        lon, lat, dst, tyme = kml.getCoords(speed=speed,t0=t0,metric=metric,
                                            noSinglePoint=noSinglePoint)
        var = self.sampleXYT(expr,lon,lat,tyme,**kwopts)
        var.grid.dst = dst # distance from start
        
        return var
    
#   ..................................................................
    interp = _interpXY  # deprecated, use sample functions instead
#   ..................................................................

    def _getDatetime(self,t):
        """
        Return datetime given grads time index "t" or "time"
        Side effect: the grads time is set to "t".
        """
        if type(t) == type(1):
            self.cmd('set t %d'%t,Quiet=True)
        elif type(t) == type("abc"):
            self.cmd('set time %s'%t,Quiet=True)
        qh = self.query("dims",Quiet=True)
        return gat2dt(qh.time[0])

#   ..................................................................

    def _tightDomain(self,lons,lats):
        """
        Reduce the (lat,lon_ domain as to bracket the coordinates
        on input.
        Side effect: dimension environment is modified.
        """
        self.cmd('set lon %f %f'%(lons.min(),lons.max()),Quiet=True)
        self.cmd('set lat %f %f'%(lats.min(),lats.max()),Quiet=True)
        fh = self.query("file",Quiet=True)
        qh = self.query("dims",Quiet=True)
        x1, x2 = (qh.xi[0]-1,qh.xi[1]+1)
        y1, y2 = (max(1,qh.yi[0]-1),min(fh.ny,qh.yi[1]+1)) # in [1,ny]
        self.cmd('set x %d %d'%(x1,x2),Quiet=True)
        self.cmd('set y %d %d'%(y1,y2),Quiet=True)
   
#.....................................................................

def _norm(x):
    """
    L-2 norm, internal use
    """
    return sqrt(inner(x,x))

def _scatter(u,I,shp,fill_value):
    """
    Scatter input array according to index mask I.
    """
    v = ma.masked_array(data=zeros(shp,dtype='float32')+fill_value,
                        mask=ones(shp,dtype='bool'),
                        fill_value=fill_value)
    v = v.squeeze()
    v.data[:,I.squeeze()] = u[:,:]
    v.mask[:,I.squeeze()] = False
    
    return v

def interpolate(datain,xin,yin,xout,yout,checkbounds=False,masked=False,order=1):
    """
    Note: This function borrowed from basemap. Reproduced here to remove basemap
          dependency.
    
    Interpolate data (``datain``) on a rectilinear grid (with x = ``xin``
    y = ``yin``) to a grid with x = ``xout``, y= ``yout``.

    .. tabularcolumns:: |l|L|

    ==============   ====================================================
    Arguments        Description
    ==============   ====================================================
    datain           a rank-2 array with 1st dimension corresponding to
                     y, 2nd dimension x.
    xin, yin         rank-1 arrays containing x and y of
                     datain grid in increasing order.
    xout, yout       rank-2 arrays containing x and y of desired output grid.
    ==============   ====================================================

    .. tabularcolumns:: |l|L|

    ==============   ====================================================
    Keywords         Description
    ==============   ====================================================
    checkbounds      If True, values of xout and yout are checked to see
                     that they lie within the range specified by xin
                     and xin.
                     If False, and xout,yout are outside xin,yin,
                     interpolated values will be clipped to values on
                     boundary of input grid (xin,yin)
                     Default is False.
    masked           If True, points outside the range of xin and yin
                     are masked (in a masked array).
                     If masked is set to a number, then
                     points outside the range of xin and yin will be
                     set to that number. Default False.
    order            0 for nearest-neighbor interpolation, 1 for
                     bilinear interpolation, 3 for cublic spline
                     (default 1). order=3 requires scipy.ndimage.
    ==============   ====================================================

    .. note::
     If datain is a masked array and order=1 (bilinear interpolation) is
     used, elements of dataout will be masked if any of the four surrounding
     points in datain are masked.  To avoid this, do the interpolation in two
     passes, first with order=1 (producing dataout1), then with order=0
     (producing dataout2).  Then replace all the masked values in dataout1
     with the corresponding elements in dataout2 (using numpy.where).
     This effectively uses nearest neighbor interpolation if any of the
     four surrounding points in datain are masked, and bilinear interpolation
     otherwise.

    Returns ``dataout``, the interpolated data on the grid ``xout, yout``.
    """
    import numpy as np
    import numpy.ma as ma
    # xin and yin must be monotonically increasing.
    if xin[-1]-xin[0] < 0 or yin[-1]-yin[0] < 0:
        raise ValueError('xin and yin must be increasing!')
    if xout.shape != yout.shape:
        raise ValueError('xout and yout must have same shape!')
    # check that xout,yout are
    # within region defined by xin,yin.
    if checkbounds:
        if xout.min() < xin.min() or \
           xout.max() > xin.max() or \
           yout.min() < yin.min() or \
           yout.max() > yin.max():
            raise ValueError('yout or xout outside range of yin or xin')
    # compute grid coordinates of output grid.
    delx = xin[1:]-xin[0:-1]
    dely = yin[1:]-yin[0:-1]
    if max(delx)-min(delx) < 1.e-4 and max(dely)-min(dely) < 1.e-4:
        # regular input grid.
        xcoords = (len(xin)-1)*(xout-xin[0])/(xin[-1]-xin[0])
        ycoords = (len(yin)-1)*(yout-yin[0])/(yin[-1]-yin[0])
    else:
        # irregular (but still rectilinear) input grid.
        xoutflat = xout.flatten(); youtflat = yout.flatten()
        ix = (np.searchsorted(xin,xoutflat)-1).tolist()
        iy = (np.searchsorted(yin,youtflat)-1).tolist()
        xoutflat = xoutflat.tolist(); xin = xin.tolist()
        youtflat = youtflat.tolist(); yin = yin.tolist()
        xcoords = []; ycoords = []
        for n,i in enumerate(ix):
            if i < 0:
                xcoords.append(-1) # outside of range on xin (lower end)
            elif i >= len(xin)-1:
                xcoords.append(len(xin)) # outside range on upper end.
            else:
                xcoords.append(float(i)+(xoutflat[n]-xin[i])/(xin[i+1]-xin[i]))
        for m,j in enumerate(iy):
            if j < 0:
                ycoords.append(-1) # outside of range of yin (on lower end)
            elif j >= len(yin)-1:
                ycoords.append(len(yin)) # outside range on upper end
            else:
                ycoords.append(float(j)+(youtflat[m]-yin[j])/(yin[j+1]-yin[j]))
        xcoords = np.reshape(xcoords,xout.shape)
        ycoords = np.reshape(ycoords,yout.shape)
    # data outside range xin,yin will be clipped to
    # values on boundary.
    if masked:
        xmask = np.logical_or(np.less(xcoords,0),np.greater(xcoords,len(xin)-1))
        ymask = np.logical_or(np.less(ycoords,0),np.greater(ycoords,len(yin)-1))
        xymask = np.logical_or(xmask,ymask)
    xcoords = np.clip(xcoords,0,len(xin)-1)
    ycoords = np.clip(ycoords,0,len(yin)-1)
    # interpolate to output grid using bilinear interpolation.
    if order == 1:
        xi = xcoords.astype(np.int32)
        yi = ycoords.astype(np.int32)
        xip1 = xi+1
        yip1 = yi+1
        xip1 = np.clip(xip1,0,len(xin)-1)
        yip1 = np.clip(yip1,0,len(yin)-1)
        delx = xcoords-xi.astype(np.float32)
        dely = ycoords-yi.astype(np.float32)
        dataout = (1.-delx)*(1.-dely)*datain[yi,xi] + \
                  delx*dely*datain[yip1,xip1] + \
                  (1.-delx)*dely*datain[yip1,xi] + \
                  delx*(1.-dely)*datain[yi,xip1]
    elif order == 0:
        xcoordsi = np.around(xcoords).astype(np.int32)
        ycoordsi = np.around(ycoords).astype(np.int32)
        dataout = datain[ycoordsi,xcoordsi]
    elif order == 3:
        try:
            from scipy.ndimage import map_coordinates
        except ImportError:
            raise ValueError('scipy.ndimage must be installed if order=3')
        coords = [ycoords,xcoords]
        dataout = map_coordinates(datain,coords,order=3,mode='nearest')
    else:
        raise ValueError('order keyword must be 0, 1 or 3')
    if masked and isinstance(masked,bool):
        dataout = ma.masked_array(dataout)
        newmask = ma.mask_or(ma.getmask(dataout), xymask)
        dataout = ma.masked_array(dataout,mask=newmask)
    elif masked and is_scalar(masked):
        dataout = np.where(xymask,masked,dataout)
    return dataout

