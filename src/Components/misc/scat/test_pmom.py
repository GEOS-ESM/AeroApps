#!/usr/bin/env python3

"""
Simple script to test Pete's new phase funcion loolup tables.

 pmoments is the phase matrix parameters after  Liou, eq. 5.2.118 - 5.2.122

 The way the phase quantities are defined follows from Wiscombe's 
 documentation (MIEV.doc) and Liou, eq. 5.2.112:
 pmom1 -> S1 * conj(S1)          = i1 in Liou
 pmom2 -> S2 * conj(S2)          = i2 in Liou
 pmom3 -> Re( S1 * conj(S2) )    = (i3 + i4)/2. in Liou
 pmom4 -> -Im( S1 * conj(S2) )   = -(i4 - i3)/2. in Liou

 pmoments is dimensioned (nLam,nRH,nBin,nMom,nPol)

   where nPol = 0 corresponds to pmom11
   where nPol = 1 corresponds to pmom12
   where nPol = 2 corresponds to pmom33
   where nPol = 3 corresponds to pmom34

"""

import warnings
warnings.simplefilter('ignore',DeprecationWarning)

from numpy    import linspace, zeros, arccos, pi
from pylab    import plot, xlabel, ylabel, title, legend, savefig, clf, semilogy
from pyhdf.SD import *
from scipy.special.orthogonal import legendre

if __name__ == "__main__":

    h = SD('ExtData/g5chem/x/optics_SU.v3_pmom.nc')

#   Coordinate variables
#   --------------------
    ch = 1.e9 * h.select('lambda').get()
    rh =        h.select('rh').get()
    r  =        h.select('radius').get()

#   Read phase function
#   -------------------
    pmom = h.select('pmom').get()
    npol, nmom, nr, nrh, nch = pmom.shape

#   Evaluate phase function in physical space
#   -----------------------------------------
    mu = linspace(-1., 1, 1001)
    m = mu.shape
    ich = 4
    irh = 35
    ibin = 0
    p11, p12, p33, p34 = (zeros(m),zeros(m),zeros(m),zeros(m))
    for n in range(255):
        print("Adding moment ", n)
        P = legendre(n)
        p11 += pmom[0,n,ibin,irh,ich] * P(mu)
        p12 += pmom[1,n,ibin,irh,ich] * P(mu)
        p33 += pmom[2,n,ibin,irh,ich] * P(mu)
        p34 += pmom[3,n,ibin,irh,ich] * P(mu)

#   Plot phase function
#   -------------------
    clf()
    angle = (180./pi) * arccos(mu)
    semilogy(angle,p11,label='$P_{11}$')
#    plot(angle,p12,label='$P_{12}$')
#    plot(angle,p33,label='$P_{33}$')
#    plot(angle,p34,label='$P_{34}$')
    xlabel('Scattering Angle')
    title('Phase Function - Sulfates')
    legend()
    savefig('pmom_su_v3.png')
