#!/usr/bin/env python
"""
Read an optics table
evaluate GSF with weights from table
to recreate phase matrix
"""

import os
from netCDF4 import Dataset
import numpy as np
from math import sqrt, cos, factorial
import matplotlib.pyplot as plt
import sys
from scipy.special.orthogonal import legendre

def recur_d_mn(m,n,s,theta):
    smin = max(abs(m),abs(n))
    if s < smin:
        return 0
    elif s == smin:
        if n >= m:
            xi = 1
        else:
            xi = (-1)**(m-n)
        sqrt_term = factorial(2*smin)/(factorial(abs(m-n))*factorial(abs(m+n)))
        third_term  = (1-cos(theta))**(0.5*abs(m-n))
        fourth_term = (1+cos(theta))**(0.5*abs(m+n))

        return xi*2**(-1*smin)*sqrt(sqrt_term)*third_term*fourth_term

    else:
        first_term = s*sqrt((s+1)**2 - m**2)*sqrt((s+1)**2 - n**2)
        first_term = 1/first_term
        second_term = (2*s+1)*(s*(s+1)*cos(theta) - m*n)*d_mn(m,n,s-1,theta)
        third_term = (s+1)*sqrt(s**2 - m**2)*sqrt(s**2 - n**2)*d_mn(m,n,s-2,theta)

        return first_term*(second_term - third_term)


def d_mn(m,n,s,theta,dm1=None,dm2=None):
    smin = max(abs(m),abs(n))
    if s < smin:
        return 0, 0
    elif (s == smin):
        if n >= m:
            xi = 1
        else:
            xi = (-1)**(m-n)
        sqrt_term = factorial(2*smin)/(factorial(abs(m-n))*factorial(abs(m+n)))
        third_term  = (1-cos(theta))**(0.5*abs(m-n))
        fourth_term = (1+cos(theta))**(0.5*abs(m+n))

        d = xi*(2**(-1*smin))*sqrt(sqrt_term)*third_term*fourth_term

        if (s == smin):
            return d, 0
        else:
            return d, dm1

    else:
        s = s-1
        first_term = s*sqrt((s+1)**2 - m**2)*sqrt((s+1)**2 - n**2)
        first_term = 1./first_term
        second_term = (2*s+1)*(s*(s+1)*cos(theta) - m*n)*dm1
        third_term = (s+1)*sqrt(s**2 - m**2)*sqrt(s**2 - n**2)*dm2

        d = first_term*(second_term - third_term)

        return d, dm1




if __name__ == '__main__':

    inDir =  '/home/pcastell/workspace/GAAS/src/Components/missions/A-CCP/ExtDataOsku'
    inFile = 'optics_SS.v3_7.GSFun.nc'
    irh    = -1
    irad   = 4
    ch     = 800.


    # read file
    nc = Dataset('{}/{}'.format(inDir,inFile))
    angle = nc.variables['scattering_angle'][:]
    theta      = np.radians(angle)
    pmom       = nc.variables['pmom'][:]
    pmom2       = nc.variables['pmom2'][:]
    pmatrix    = nc.variables['phase_matrix'][:]
    wav        = nc.variables['lambda'][:]*1e9   # nm

    # get wavlength index
    iwav = np.argmin(np.abs(wav-ch))

    
    nPol, nMom, nrad, nrh, nlambda = pmom.shape
    ntheta = len(theta)
    nMomL = 300
    nMom = 30

    # P11
    # ordered over nPol as P11, P12, P33, P34, P22, P44
    # special case, legendre polynomials
    ipol = 0
    gsf_coef = pmom[ipol,:,irad,irh,iwav]
    leg_coef = pmom2[ipol,:,irad,irh,iwav]    

    p11 = np.zeros(ntheta)
    p11l = np.zeros(ntheta)
    mu    = np.cos(np.radians(angle))
    for s in range(nMomL):
        P = legendre(s)
        p11 += gsf_coef[s]* P(mu)
        p11l += leg_coef[s]* P(mu)


    plt.plot(angle,p11,label='GSF')
    plt.plot(angle,p11l,label='Legendre')
    plt.legend()
    plt.title('P11')
    plt.show()

    # P44
    # ordered over nPol as P11, P12, P33, P34, P22, P44
    # special case, legendre polynomials
    ipol = 5
    gsf_coef = pmom[ipol,:,irad,irh,iwav]
    leg_coef = pmom2[ipol,:,irad,irh,iwav]    

    p44 = np.zeros(ntheta)
    p44l = np.zeros(ntheta)
    mu    = np.cos(np.radians(angle))
    for s in range(nMomL):
        P = legendre(s)
        p44 += gsf_coef[s]* P(mu)
        p44l += leg_coef[s]* P(mu)


    plt.plot(angle,p44,label='GSF')
    plt.plot(angle,p44l,label='Legendre')
    plt.legend()
    plt.title('P44')
    plt.show()


    # P12
    # ordered over nPol as P11, P12, P33, P34, P22, P44    
    ipol = 1
    m = 0.0
    n = 2.0    
    gsf_coef = pmom[ipol,:,irad,irh,iwav]
    leg_coef = pmom2[ipol,:,irad,irh,iwav]

    p12 = np.zeros(ntheta)
    for i,t in enumerate(theta):
        d0, dneg1 = d_mn(m,n,0,t)
        pfunc = d0/(1j**(n-m)).real
        p12[i] = gsf_coef[0]*pfunc 
        d1, d0 = d_mn(m,n,1,t,dm1=d0,dm2=dneg1)
        pfunc = d1/(1j**(n-m)).real
        p12[i] = p12[i] + gsf_coef[1]*pfunc  

        dm2 = d0
        dm1 = d1
        for s in range(2,nMom):
            dfunc, dm2 = d_mn(m,n,s,t,dm1=dm1,dm2=dm2)
            pfunc = dfunc/(1j**(n-m)).real

            p12[i] = p12[i] + gsf_coef[s]*pfunc  

            dm1 = dfunc

    p12l = np.zeros(ntheta)
    mu    = np.cos(np.radians(angle))
    for s in range(nMomL):
        P = legendre(s)
        p12l += leg_coef[s]* P(mu)

    plt.plot(angle,p12,label='GSF')
    plt.plot(angle,p12l,label='Legendre')
    plt.legend()
    plt.title('P12')
    plt.show()

    # P34
    # ordered over nPol as P11, P12, P33, P34, P22, P44    
    ipol = 3
    m = 0.0
    n = 2.0    
    gsf_coef = pmom[ipol,:,irad,irh,iwav]
    leg_coef = pmom2[ipol,:,irad,irh,iwav]

    p34 = np.zeros(ntheta)
    for i,t in enumerate(theta):
        d0, dneg1 = d_mn(m,n,0,t)
        pfunc = d0/(1j**(n-m)).real
        p34[i] = gsf_coef[0]*pfunc 
        d1, d0 = d_mn(m,n,1,t,dm1=d0,dm2=dneg1)
        pfunc = d1/(1j**(n-m)).real
        p34[i] = p34[i] + gsf_coef[1]*pfunc  

        dm2 = d0
        dm1 = d1
        for s in range(2,nMom):
            dfunc, dm2 = d_mn(m,n,s,t,dm1=dm1,dm2=dm2)
            pfunc = dfunc/(1j**(n-m)).real

            p34[i] = p34[i] + gsf_coef[s]*pfunc  

            dm1 = dfunc

    p34l = np.zeros(ntheta)
    mu    = np.cos(np.radians(angle))
    for s in range(nMomL):
        P = legendre(s)
        p34l += leg_coef[s]* P(mu)

    plt.plot(angle,p34,label='GSF')
    plt.plot(angle,p34l,label='Legendre')
    plt.legend()
    plt.title('P34')
    plt.show()


    # P33 & P22
    # ordered over nPol as P11, P12, P33, P34, P22, P44    
    
    # first get a2 + a3    
    m = 2.0
    n = 2.0  
    ipol = 4  
    gsf_coef22 = pmom[ipol,:,irad,irh,iwav]
    leg_coef22 = pmom2[ipol,:,irad,irh,iwav]

    ipol = 2  
    gsf_coef33 = pmom[ipol,:,irad,irh,iwav]  
    leg_coef33 = pmom2[ipol,:,irad,irh,iwav]

    a2p3 = np.zeros(ntheta)
    for i,t in enumerate(theta):
        d0, dneg1 = d_mn(m,n,0,t)
        pfunc = d0/(1j**(n-m)).real
        a2p3[i] = (gsf_coef22[0]+gsf_coef33[0])*pfunc 
        d1, d0 = d_mn(m,n,1,t,dm1=d0,dm2=dneg1)
        pfunc = d1/(1j**(n-m)).real
        a2p3[i] = a2p3[i] + (gsf_coef22[1]+gsf_coef33[1])*pfunc  

        dm2 = d0
        dm1 = d1
        for s in range(2,nMom):
            dfunc, dm2 = d_mn(m,n,s,t,dm1=dm1,dm2=dm2)
            pfunc = dfunc/(1j**(n-m)).real

            a2p3[i] = a2p3[i] + (gsf_coef22[s]+gsf_coef33[s])*pfunc  

            dm1 = dfunc

    # next get a2 - a3    
    m = 2.0
    n = -2.0  
    a2m3 = np.zeros(ntheta)
    for i,t in enumerate(theta):
        d0, dneg1 = d_mn(m,n,0,t)
        pfunc = d0/(1j**(n-m)).real
        a2m3[i] = (gsf_coef22[0]-gsf_coef33[0])*pfunc 
        d1, d0 = d_mn(m,n,1,t,dm1=d0,dm2=dneg1)
        pfunc = d1/(1j**(n-m)).real
        a2m3[i] = a2m3[i] + (gsf_coef22[1]-gsf_coef33[1])*pfunc  

        dm2 = d0
        dm1 = d1
        for s in range(2,nMom):
            dfunc, dm2 = d_mn(m,n,s,t,dm1=dm1,dm2=dm2)
            pfunc = dfunc/(1j**(n-m)).real

            a2m3[i] = a2m3[i] + (gsf_coef22[s]-gsf_coef33[s])*pfunc  

            dm1 = dfunc

    # now combine a2 + a3 = a2p3 with a2 - a3 = a2m3
    p22 = 0.5*(a2p3 + a2m3)
    p33 = a2p3 - p22

    p22l = np.zeros(ntheta)
    p33l = np.zeros(ntheta)
    mu    = np.cos(np.radians(angle))
    for s in range(nMomL):
        P = legendre(s)
        p22l += leg_coef22[s]* P(mu)
        p33l += leg_coef33[s]* P(mu)

    plt.plot(angle,p22,label='GSF')
    plt.plot(angle,p22l,label='Legendre')
    plt.legend()
    plt.title('P22')
    plt.show()

    plt.plot(angle,p33,label='GSF')
    plt.plot(angle,p33l,label='Legendre')
    plt.legend()
    plt.title('P33')
    plt.show()    
