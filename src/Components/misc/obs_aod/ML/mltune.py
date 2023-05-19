"""
   The *mltune* module implements an extension of the ODS class useful for Maximum-likelihhod
tuning of innovation statistics as in Dee and da Silva (1999).

"""
from numpy import *
from numpy.linalg.linalg import cholesky
from scipy.optimize    import fmin 

import pyods

a = 6376.e3 # earth's radius 

verbose = True

class MLTUNE(pyods.ODS):

    def rdist(self):
        """
    R = RDIST() returns the pairwise chordal distance 
    matrix for the earth surface locations defined by the 
    longitude-latitude attributes *lon*, *lat*.

    The lon-lat coordinates are assumed to be given in degrees; distances
    are returned in meters.

    See also static method *rdist*.
        """

        return rdist(self.lon,self.lat)

#                        --------------
#                        Static Methods
#                        --------------

def rdist(lon,lat,Chordal=False):
    """
    RDIST	Chordal/Great Circle distance.

    R = RDIST(LON,LAT) returns the pairwise chordal or great circle
        distance matrix for the earth surface locations defined by the
        longitude-latitude coordinate arrays LON,LAT

                   R(i,j) = dist(u(i),u(j)) 

        where u(i) has lon-lat coordinates LON(i),LAT(i), u(j) has
        lon-lat coordinates LON(j),LAT(j), and dist(.,.)  is chordal
        or great circle distance depending on the parameter *Chordal*.
		
       The lon-lat coordinates are assumed to be given in degrees; distances
       are returned in meters.

       """

    n = lon.size

#   Cartesian coords on unity sphere 
#   --------------------------------
    slon = sin((pi/180)*lon[:])
    clon = cos((pi/180)*lon[:])
    slat = sin((pi/180)*lat[:])
    clat = cos((pi/180)*lat[:])
    x, y, z = (clat*clon, clat*slon, slat) 

#   Compute distances on unit sphere
#   --------------------------------
    R = zeros((n,n))
    if Chordal:
        for i in range(n):
            xi, yi, zi = ( x[i], y[i], z[i] )
            r2 = (x-xi)**2 + (y-yi)**2 + (z-zi)**2
            r2[r2<0.0] = 0.0
            R[i] = sqrt(r2)
    else:
        for i in range(n):
            xi, yi, zi = ( x[i], y[i], z[i] )
            c = x*xi + y*yi + z*zi
            c[c>1.]  = 1.
            c[c<-1.] = -1.
            R[i] = arccos(c)

#   Multiply by Earth's radius
#   --------------------------
    R = a * R

    return R

#--
def powerlaw(r,L,extra=False):
    """
    POWERLAW Powerlaw function and derivatives.

       f = powerlaw(R,L) is the powerlaw function with length scale
       parameter L evaluated at R.

       The length scale L is defined here as L = sqrt(-1/f''(0)).

       f,df,ddf = powerlaw(R,L,extra=True) also returns df = (1/r) df/dr and
       ddf = r d/dr (1/r) df/dr = r d/dr df.

       Based on Matlab script of 14Nov97 by Dick Dee.
    """

    f = 1./(1. + 0.5*(r/L)**2)

    if extra:
        df  = (-1/(L**2))*(f**2);
        ddf = (-2*(r/L)**2)*(f*df);
        return (f,df,ddf)
    else:
        return f


#--
def cspline(r,L,extra=False):
    """
    CSPLINE Compactly supported spline function and derivatives.

       f = cspline(R,L) is the compactly supported spline function
       with length scale parameter L evaluated at R.

       The length scale L is defined here as L = sqrt(-1/f''(0)).

      f,df,ddf = cspline(R,L,extra=True) also returns df = (1/r) df/dr and
      ddf = r d/dr (1/r) df/dr = r d/dr df.

    Based on Matlab script of 14Nov97 by Dick Dee.
    """

    s = r.shape
    c = sqrt(10./3.)*L
    z = r.ravel()/c
    n = r.size

    if extra is False:

        f = zeros(n)

        i = (z<=1.)

        f[i] =  (- 1./4.)*z[i]**5 + ( 1./2.)*z[i]**4 \
             + ( 5./8.)*z[i]**3 + (- 5./3.)*z[i]**2 + 1

        i = (1<z)&(z<2)

        f[i] =  ( 1./12.)*z[i]**5 + (-1./2.)*z[i]**4  \
             + ( 5./8.)*z[i]**3 + (  5./3.)*z[i]**2 \
             + (-5.)*z[i]       + (4.) + (-2./3.)*z[i]**(-1)     

        return f.reshape(s)

    else:

        f   = zeros(n)
        df  = zeros(n)
        ddf = zeros(n)
        
        i = find(z<=1)

        f[i]   =  (- 1./4.)*z[i]**5 + ( 1./2.)*z[i]**4 \
               + ( 5./8.)*z[i]**3 + (- 5./3.)*z[i]**2 + (1) 
        df[i]  = ((- 5./4.)*z[i]**3 + ( 2  )*z[i]**2 \
               + (15./8.)*z[i]    + (-10./3.)               )/c^2
        ddf[i] = ((-15./4.)*z[i]**3 + ( 4  )*z[i]**2 \
               + (15./8.)*z[i]                            )/c^2
        
        i = find((1<z)&(z<2))

        f[i]   =  ( 1./12.)*z[i]**5 + (-1./2.)*z[i]**4 \
               + ( 5./8.)*z[i]**3 + (  5./3.)*z[i]**2 + (-5)*z[i]       \
               + (4) + (-2./3.)*z[i]**(-1)     
        df[i]  = (( 5./12.)*z[i]**3 + (-2  )*z[i]**2 + (15./8.)*z[i]    \
               + ( 10./3.)         + (-5)*z[i]**(-1) +       ( 2./3.)*z[i]**(-3))/c^2
        ddf[i] = (( 5./4. )*z[i]**3 + (-4  )*z[i]**2 + (15./8.)*z[i]    \
               +                   ( 5)*z[i]**(-1) +       (-6./3.)*z[i]**(-3))/c^2

        return (f.reshape(s),df.reshape(s),ddf.reshape(n))


#--
def winplaw(r,L,rstar=6.e6,extra=False):
    """
    WINPLAW Windowed powerlaw function.

       f = winplaw(R,L,RSTAR) is the windowed powerlaw function with
       length scale parameter L evaluated at R, and support parameter
       RSTAR. See DAO Office Note 97-??, Appendix A.

       The function is identically zero where R>RSTAR. If the third
       parameter is missing, then RSTAR = 6e6.

       The length scale L is defined here as L = sqrt(-1/f''(0)).

       f,df,ddf = powerlaw(R,L,extra=True) also returns df = (1/r) df/dr and
       ddf = r d/dr (1/r) df/dr = r d/dr df.

       Based on Matlab script of 14Nov97 by Dick Dee.
    """

    L2 = (rstar/2.) * sqrt(3./10.)
    L1 = L/sqrt(1. - (40./3.)*(L/rstar)**2)

    if extra is False:

        g = powerlaw(r,L1)
        h =  cspline(r,L2)
        f = g*h 

        return f

    else:

        g, dg, ddg = powerlaw(r,L1,extra=True)
        h, dh, ddh = cspline(r,L2,extra=True)
       
        f   = g*h
        df  = g*dh  + h*dg
        ddf = g*ddh + r*dg*r*dh + h*ddg

        return (f,df,ddf)

#--
def cov_generic(alpha,v,corrfun):
    """
     Covariance model for scalar residuals based on a generic
     corrrelation function. On input,

        alpha = (sigO, sigF, L)

        v --- ODS object for a single synoptic time

        corrfun --- correlation function

    """
    lon = v.lon     # list of longitudes
    lat = v.lat     # list of latitudes
    n = lon.size

    varO = alpha[0]**2        # observation error standard deviation
    varF = alpha[1]**2        # forecast error standard deviation
    L    = alpha[2]*1e6       # decorrelation length scale (L in mega meter)

    r = rdist(lon,lat)

    S = varO * eye(n) + varF * corrfun(r,L)
     
    return S

#--
def cov_plaw(alpha,v):
    """
    Power Law covariance model. On input, v is an ODS object object and 
      sigO = alpha[0]
      sigF = alpha[1]
      L    = alphap[2]
    """
    if len(alpha) != 3:
        raise ValueError("alpha must have size 3")
    return cov_generic(alpha,v,powerlaw)

#--
def cov_winplaw(alpha,v):
    """
    Compactly Suported Power Law covariance model. On input, v is an ODS 
    object object and 
      sigO = alpha[0]
      sigF = alpha[1]
      L    = alphap[2]
    """
    if len(alpha) != 3:
        raise ValueError("alpha must have size 3")
    return cov_generic(alpha,v,winplaw)

#--
def cor_plaw(alpha,v):
    """
    Power Law correlation model. On input, v is an object object and 
      sigO = cos(alpha[0])
      sigF = sin(alpha[0])
      L    = alphap[1]
    """
    if len(alpha) != 2:
        raise ValueError("alpha must have size 2")
    sigO = cos(alpha[0])
    sigF = sin(alpha[0])
    alpha_ = (sigO,sigF,alpha[1])
    return cov_generic(alpha_,v,powerlaw)

#--
def cor_winplaw(alpha,v):
    """
    Compactly Supported Power Law correlation model. On input, v is 
    an object object and 
      sigO = cos(alpha[0])
      sigF = sin(alpha[0])
      L    = alphap[1]
    """
    if len(alpha) != 2:
        raise ValueError("alpha must have size 2")
    sigO = cos(alpha[0])
    sigF = sin(alpha[0])
    alpha_ = (sigO,sigF,alpha[1])
    return cov_generic(alpha_,v,winplaw)

#--
def llfun(alpha, covmodel, V):
    """
    LLFUN Log-likelihood function for covariance parameter estimation:

                 F = llfun(alpha, covmodel, v) 

    evaluates the log-likelihood function for data from a white
    multivariate Gaussian time series.

       alpha      a vector of parameter values
       covmodel   function which evaluates the covariance
                  model as a function of alpha and v, e.g.,
                  cov_plaw.
       v          a list of ODS objects, one for each time

     Based on Matlab script of 05Dec97 by Dick Dee.
     """

    global verbose

    N = 0
    f = 0.0

    if verbose:
        k = 0
        print('Evaluating Log-likelihood at sigO=%6.3f, sigF=%6.3f, L=%5.1f'\
               %(alpha[0],alpha[1],1000.*alpha[2]))

#   Loop over time
#   --------------
    for v in V:

#        if verbose:
#            k = k + 1
#            print '  --> time step ', k

#       Handle case of no observations
#       ------------------------------
        if v.nobs < 1: continue

#       Evaluate the covariance model
#       -----------------------------
        Sk = covmodel(alpha, v)

        vk = v.omf       # data for this k (column vector)

#       The following four lines are equivalent to
#             fk = log(det(Sk)) + vk'*(Sk\vk)
#       but compute somewhat faster
#       ------------------------------------------
        Uk = cholesky(Sk)
        sk = solve(Uk,vk) 
        dk = diag(Uk)
        fk = sum(log(dk*dk)) + dot(sk,sk)

        f = f + fk            # accumulate the cost function
        N = N + sk.size

    if N>0: f = f / N         # normalize cost function 

    return f

#--
def llhess(alpha, covmodel, v):
    """
    LLHESS --- Hessian of Log-likelihood function estimation.

    Hessian = llhess(alpha, covmodel, v) returns a finite-difference
                 approximation of the Hessian.

       alpha     a vector of parameter values
       covmodel  function which evaluates the covariance
                 model as a function of A and V
       v         a list of ODS objects, one for each time

     Based on Matlab script of 05Dec97 by Dick Dee.
     """

#   Approximate the Hessian of f at alpha:
#   -------------------------------------
    na = len(alpha)
    da = alpha*1e-2
    hessf = zeros((na,na))
    for j in range(na):
        ap = alpha.copy();    ap[j] = alpha[j] + da[j]
        am = alpha.copy();    am[j] = alpha[j] - da[j]
        for i in range(j):
            app = ap.copy(); app[i] = ap[i] + da[i]
            apm = ap.copy(); apm[i] = ap[i] - da[i]
            amp = am.copy(); amp[i] = am[i] + da[i]
            amm = am.copy(); amm[i] = am[i] - da[i]
            fpp = llfun(app, covmodel, v)
            fpm = llfun(apm, covmodel, v)
            fmp = llfun(amp, covmodel, v)
            fmm = llfun(amm, covmodel, v)
            hessf[i,j] = (fpp - fpm - fmp + fmm)/(4*da[i]*da[j])
    for j in range(na):
        for i in range(j+1,na):
            hessf[i,j] = hessf[j,i]

    return hessf

#--
def estimate(alpha0,ods_ts,covmodel,**kwopts):
    """
    Maximum-likelihood estimation of covariance parameters.

        alpha, err = ml.etimate(alpha0,covmodel,ods_ts,**kwopts)

    This function produces maximum-likelihood estimates of covariance
    parameters for a multivariate timeseries. Required paramaters:

       alpha0    ---  Initial guess for parameters
       
       covmodel  ---  covariance function, e.g., cov_plaw 

       ods_ts    ---  1D array of ODS objects, each element 
                      corresponding to a time instant. Everytime
                      may have a different set of observation locations.

    The optional keyword arguments (**kwopts) are passed to the Nelder-Mead
    minimization function scipy.optimize.fmin().
       
    Based on Matlab script of 14Nov97 by Dick Dee.
    """

#   perform the optimization
#   ------------------------
    alpha = fmin(llfun,alpha0,args=(covmodel,ods_ts),**kwopts)

#   evaluate the log-likelihood function and the Hessian at the minimum:
#   -------------------------------------------------------------------
    f = llfun (alpha, covmodel, ods_ts)
    H = llhess(alpha, covmodel, ods_ts)

#   Approximate the error covariance of the parameter estimates:
#   -----------------------------------------------------------
    n = 0
    for k in range(len(ods_ts)):
        n = n + ods_ts[k].omf.size
    covalpha = inv(H/2)/n;
    sigalpha = sqrt(abs(diag(covalpha)))

#   display the results:
#   -------------------
    print("Minimum cost function value: %6.4f"%f)
    print("Parameter estimates and standard errors: ")
    alpha[2] = alpha[2] * 1000. 
    sigalpha[2] = sigalpha[2] * 1000. 
    alphan = ( 'sigO', 'sigF', '   L')
    for j in range(alpha.size):
        print("    %s = %6.2f +/- %4.3f"%(alphan[j], alpha[j], sigalpha[j]))
    alpha[2] = alpha[2] / 1000. 
    sigalpha[2] = sigalpha[2] / 1000. 

    return (alpha, sigalpha)

def getODSts(ods_tmpl,NYMD,NHMS,N=0,BBox=None,qch=None):
    """
    Given an ODS filename template, and a sequence of dates and times, returns a
    sequence of ODS objects. Optionally, the number of obs in each element
    of the sequence can be capped by a number N. The default is not 
    to cap. In addition, a bounding box BBox can be specified to select a
    particular region of the globe. BBox is a tuple of the form
                (lon_min,lat_min,lon_max,lat_max)
    """

    # ods = pyods.ODS(ods_tmpl) # make sure this is a bonafide ODS file

    if verbose:
        print("  Date    Time    NOBS     nobs      N")
        print("-------- ------ -------- -------- --------")

    V = ()
    for nymd in NYMD:
        cnymd = str(nymd)
        for nhms in NHMS:
            cnhms = str(nhms)
            chh = str(nhms/10000)
            fname = ods_tmpl.replace('%nymd',cnymd).\
                             replace('%nhms',cnhms).\
                             replace('%hh',chh)
            ods = pyods.ODS(fname,nymd,nhms,only_good=True)
            nobs = ods.nobs
            if ods.nobs < 1: continue
            if BBox is not None:
                I = (ods.lon>=BBox[0]) & (ods.lon<=BBox[2]) & \
                    (ods.lat>=BBox[1]) & (ods.lon<=BBox[3])
                if any(I):
                    ods = ods.__copy__(Indices=I)
                else:
                    continue
            if qch is not None:
                ods = ods.select(qch=qch)
            nobs2 = ods.nobs
            if N > 0:
                if ods.nobs < N: continue # must have at least N observations
                n = min(ods.nobs,N)
                ods = ods.shuffle().__copy__(Indices=list(range(n)))
            if ods.nobs > 10:
                if verbose:
                    print("%8d %6d %8d %8d %8d"%(nymd,nhms,nobs,nobs2,ods.nobs))
                V = V + (ods,)

    return V

#................................................................................

if __name__ == "__main__":

    r = 1.e6 * linspace(0.,8.,512)
    L = 1.5e6

    pl = powerlaw(r,L)
    wpl = winplaw(r,L)

    clf()
    plot(r/1e3,pl,label='Power Law')
    plot(r/1e3,wpl,label='Windowed Power Law')

    title('Correlation Functions')
    xlabel('Distance (km)')
    legend(loc="upper right")
