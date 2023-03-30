#!/usr/bin/env python3
"""
Checks the effects of temperature variance inside a fire.
"""

from .planck import *
from pylab import *


def do(sig=(50.,75.,100.)):

    T = linspace(400.,800.,400)

    L21 = B21(T)
    L31 = B31(T)

    nL21 = [nB21(T,s) for s in sig] 
    nL31 = [nB31(T,s) for s in sig]

    plot(L31,L21,'k',nL31[0],nL21[0],'b',nL31[1],nL21[1],'r',nL31[2],nL21[2],'g')
    xlabel('L$_{31}$')
    ylabel('L$_{21}$')
    title('Gaussian Convoluted Planck Functions: $\sigma$=(%3.0fK,%3.0fK,%3.0fK)'%sig)
    savefig('gaussian_fires_a.png')

    clf()
    plot(iB31(L31),iB21(L21),'k',
         iB31(nL31[0]),iB21(nL21[0]),'b',
         iB31(nL31[1]),iB21(nL21[1]),'r',
         iB31(nL31[2]),iB21(nL21[2]),'g')
    xlabel('T$_{31}$ (K)')
    ylabel('T$_{21}$ (K)')
    title('Gaussian Convoluted Planck Functions: $\sigma$=(%3.0fK,%3.0fK,%3.0fK)'%sig)
    savefig('gaussian_fires_b.png')
