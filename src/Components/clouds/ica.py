"""
 Interface to ICA code.
"""

import scipy.stats as stats
import numpy as np
import scipy as sp
import time

from ICA_ import *
from eta  import getPe
from MAPL.constants import MAPL_GRAV, MAPL_RADIUS

# ---
# correlation functions used by clumpICA ...

def _exponential_CF(r):
  return np.exp(-r)

def _Gaussian_CF(r):
  print 'corr: %s seems to trigger a non-positive def error' % corr
  raise RuntimeError, 'corr: %s currently disallowed' % corr
  return np.exp(-0.5*r**2)

def _SOAR_CF(r):
  # Second-Order Auto-Regressive
  # Gaspari & Cohn (2.34)
  return (1+r)*np.exp(-r)

def _TOAR_CF(r):
  # Third-Order Auto-Regressive
  # Gaspari & Cohn (4.8)
  return (1+r*(1+r/3.))*np.exp(-r)

def _TriSCnvR3_CF(r):
  # Self-convolution of Triangle function over R^3 (compactly supported on [0,2])
  # Gaspari & Cohn (4.10)
  # TO DO: make more efficient by hieracy of parentheses
  i1 = (r <= 1.); i3 = (r >= 2.)
  i2 = np.logical_not(np.logical_or(i1, i3))
  c23 = 2/3.; c53 = 5/3.; r12 = 1/12.
  C[i1] = -0.25*r[i1]**5 + 0.5*r[i1]**4 + 0.625*r[i1]**3 - c53*r[i1]**2            + 1.
  C[i2] =   r12*r[i2]**5 - 0.5*r[i2]**4 + 0.625*r[i2]**3 + c53*r[i2]**2 - 5.*r[i2] + 4. - c23/r[i2]
  C[i3] = 0.
  return C

# ---

def genICA(ncols,DELP,T,QV,QL,QI,CLOUD,PTOP,mode,
           Longitude=None, Latitude=None,
           plim=100.,Lp=100.,Ls=0.1):
    """
    Generate *ncols* independent columns.
    Optionally add spatial coherence (see below).
    
    QV_, QL_, QI_, rc = genICA(ncols,DELP,T,QV,QL,QI,CLOUD,PTOP,mode)

    On input,

           ncols   --- number of columns to generate
           DELP    --- pressure thickness of each layer [Pa]
           T       --- temperature [K]
           QV      --- specific humidity [kg/kg]
           QL      --- specific cloud liquid water [kg/kg]
           QI      --- specific cloud ice water [kg/kg]
           CLOUD   --- cloud fraction [0-1]
           PTOP    --- top pressure [Pa]
           mode    --- mode of ICA generation:
             'HOMOCLD-MAXRAN':     homogeneous clouds, maximum-random overlap
             'TOTWPDF-GCOP-SKEWT': total water PDF, Gaussian copula overlap, Skewed-Triangle PDF

    Optionally,

         Regular TOTWPDF-GCOP-type ICA:
           plim    --- variability from surface to plim [hPa]
           Lp      --- vertical decorrelation length [hPa]
           Ls      --- state dependent decorrelation length
                         (Riishojgaard's formulation)

         For an ICA plus spatial coherence, add these:
           Longitude --- longitude & latitude of ncols columns.
           Latitude        If either None, just do regular ICA.
                           Otherwise, must have shape (ncols,).

    On output,

     QV_, QL_, QI_ --- water vapor and cloud condensate for each subcolum.
                       (ncols, nlayers)

     rc            --- return code, 0 if success
    """

    # default success
    rc = 0

    # ICA generation
    # --------------
    if mode == 'TOTWPDF-GCOP-SKEWT':
      plim_,Lp_ = (100.*plim,100.*Lp) # hPa->Pa
      QV_,QL_,QI_,rc = genica(ncols,PTOP,DELP,T,QV,QL,QI,CLOUD,plim_,Lp_,Ls)
      if rc: print 'Error on return from genICA, rc %d'%rc
    elif mode == 'HOMOCLD-MAXRAN':
      pref = getPe(len(DELP))
      QV_,QL_,QI_,rc = genica_geos5like(ncols,PTOP,DELP,T,QV,QL,QI,CLOUD,pref)
      if rc: print 'Error on return from genICA_GEOS5like, rc %d'%rc
    elif mode == 'HOMOCLD-COSP':
      QV_,QL_,QI_,rc = genica_cosp(ncols,PTOP,DELP,T,QV,QL,QI,CLOUD)
      if rc: print 'Error on return from genICA_COSP, rc %d'%rc
    else:
      raise ValueError, 'unknown ICA mode %s'%mode

    # handle errors more gently
    if rc:
      print '  ... dropping this generation'
      return (None, None, None, rc)

    # Optional spatial clumping
    # -------------------------
    if (Longitude is not None) and (Latitude is not None):
      QV_,QL_,QI_ = clumpICA(ncols,QV_,QL_,QI_,Longitude,Latitude,
        np.sum((QL_+QI_)*np.tile(DELP/MAPL_GRAV,(ncols,1)),axis=1)) # CWP

    return (QV_, QL_, QI_, rc)

#---

def clumpICA(ncols,QV_,QL_,QI_,
             Longitude,Latitude,svar,
             corrfn=_SOAR_CF, 
             Lh=5.):
    """
    Generate clumped (spatially coherent) columns from regular ICA columns.
    
        QV_, QL_, QI_ = clumpICA(ncols,QV_,QL_,QI_,Longitude,Latitude,svar)

    On input,

        ncols     --- number of columns                     scalar
        QV_       --- specific humidity                     (ncols,nlayers)
        QL_       --- specific cloud liquid water           (ncols,nlayers)
        QI_       --- specific cloud ice water              (ncols,nlayers)
        Longitude --- longitude [deg]                       (ncols,)
        Latitude  --- latitude [deg]                        (ncols,)
        svar      --- variable to use for spatial clumping  (ncols,)

    Optionally,

        corrfn --- correlation function from _CF() above    function
        Lh     --- horiz correlation length scale [km]      scalar

    On output,

        QV_, QL_, QI_ --- clumped subcolumns                (ncols,nlayers)

    """

    # sanity checks
    if not (QV_.shape == QL_.shape == QI_.shape):
      raise ValueError, 'spatial: QV_,QL_,QI_ not same shape!'
    if QV_.shape[0] != ncols or QV_.ndim != 2:
      raise ValueError, 'spatial: QV_,QL_,QI_ need shape (ncols,nlayers)'
    if not (Longitude.shape == Latitude.shape == svar.shape == (ncols,)):
      raise ValueError, 'spatial: Lon, Lat and svar need shape (ncols,)'
    if not (Lh > 0.): raise ValueError, 'spatial: require Lh > 0.'

    # distances between input locations [m]
    R = _rdist(Longitude,Latitude,
      Chordal=True) # FASTER

    # GCOP correlation matrix
    C = corrfn(R/(Lh*1.e3))

    # Cholesky decomposition
    # >> numpy version is MUCH slower
    # >> overwrite_a=True may be a little faster,
    #      but destroys C (which is OK here)
    H = sp.linalg.cholesky(C,lower=True,overwrite_a=True); del C

    # GCOP rank [0,1] generation
    rank = stats.norm.cdf(np.dot(H,np.random.randn(ncols))); del H

    # convert rank to (0,1,...,ncols-1)
    rank = np.clip(np.floor(rank*ncols).astype(int), 0, ncols-1)

    # find the index to sort variable svar
    ix = np.argsort(svar)

    # now sample from index according to GCOP ranks
    ix = ix[rank]

    # ix is now the index which spatially clumps svar
    # we apply it to sample all outputs accordingly
    QV_ = QV_[ix]; QL_ = QL_[ix]; QI_ = QI_[ix]

    # note on how it all works:
    # two close by points will have a high correlation and therefore similar GCOP ranks.
    # they will therefore choose nearby values in the sorted list, and therefore give
    # similar values of var. 

    return (QV_, QL_, QI_)

#---

def getRe(DELP,T,U,QILS,QIAN,PTOP):
    """
    Return Re in meters.
    """

    pref = getPe(DELP.shape[-1])

    REL,REI = simre(PTOP,DELP,T,pref,U,QILS,QIAN)

    return (REL,REI)

# ---

def getTau(DELP,REL,REI,QL,QI):
    
    TAUL,TAUI = simtau(DELP,REL,REI,QL,QI)
                       
    return (TAUL,TAUI)

# ---

def getQsat(PTOP,DELP,T):
    
    return qsatmcs(PTOP,DELP,T)
                       
# ---

def _rdist(lon,lat,Chordal=False):
    """
    RDIST       Chordal/Great Circle distance.

    R = RDIST(LON,LAT) returns the pairwise chordal or great circle
        distance matrix for the earth surface locations defined by
        the longitude-latitude coordinate arrays LON, LAT (each 1D):

                   R(i,j) = dist(u(i),u(j))

        where u(i) has lon-lat coordinates (LON(i),LAT(i)), u(j) has
        lon-lat coordinates (LON(j),LAT(j)), and dist(.,.) is chordal
        or great circle distance depending on the parameter *Chordal*.

       The lon-lat coordinates are assumed to be given in degrees;
       distances are returned in meters.

       """

    if lon.ndim != 1 or lat.ndim != 1:
      raise ValueError, 'LON, LAT must each be 1D arrays'
    n = lon.size

#   Cartesian coords on unit sphere
#   -------------------------------
    slon = np.sin(np.radians(lon))
    clon = np.cos(np.radians(lon))
    slat = np.sin(np.radians(lat))
    clat = np.cos(np.radians(lat))
    x, y, z = (clat*clon, clat*slon, slat)

#   Compute distances on unit sphere
#   --------------------------------
    R = np.zeros((n,n))

    # Distances are symmetric with zero diagonal.
    # So do lower triangle only, diag remains zero.
    for i in range(1,n):

      # We start by calculating the chordal distance on
      # the unit sphere. Even the great circle distance
      # will be based on this.
      r = np.sqrt(np.clip(
        (x[0:i]-x[i])**2 + (y[0:i]-y[i])**2 + (z[0:i]-z[i])**2,
        0.0, 4.0))
      if Chordal:
        pass
      else:
        # great circle angle
        r = 2.*np.arcsin(r/2.)
      # multiply by Earth's radius
      r = MAPL_RADIUS * r
      # load lower triangle
      R[i][0:i] = r

    # symmetrize
    R = R + R.T

    return R

# ---
