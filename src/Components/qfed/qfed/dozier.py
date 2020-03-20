#!/usr/bin/env python

"""
Implements Dozier type algorihms for estimating fire size/temperature.

This software is hereby placed in the public domain.
Arlindo.daSilva@nasa.gov

"""

import sys

from mxd14  import *
from planck import *

from math               import pi
from pylab              import pcolor, plot, colorbar, axis, savefig, subplot, clf, \
                               title, xlabel, ylabel, prctile, median
from scipy.stats        import kde
from scipy.optimize     import fixed_point, brent, fmin
from matplotlib.mlab    import prctile, find
from numpy.core.numeric import isscalar

igbp_dir  = '/nobackup/Emissions/Vegetation/GL_IGBP_INPE'

P_SCALE = 1e4

class DOZIER(MxD14_L2):

    def classic_var(self,tau21=0.864,tau31=0.864,Verbose=False):
        """
        Implements the classic Dozier algorithm. On input, *tau_21* and
        *tau_31* are the atmospheric transmittances at 4 and 11 microns
        (MODIS channels 21 and 31). Variational version.
        """

        self.algo = 'dozier'

#       Compute radiances
#       -----------------
        L21 = B21(self.T21) 
        L31 = B31(self.T31) 
        E21 = B21(self.Tb21)
        E31 = B31(self.Tb31)
        N = L21.size

        if isscalar(tau21): tau21 = tau21 * ones(L21.shape)
        if isscalar(tau31): tau31 = tau31 * ones(L31.shape)
            
#       Use a variational approach - Needs vectorization
#       ------------------------------------------------
        sig21 = 1.
        sig31 = 1.
        x0 = [600.,P_SCALE*0.1] # [Tf,p]; p here is normalized; 10 ha/ 1 km2 = 0.1
        Tf = - ones(N)
        p = - ones(N)
        niter = 200
        for i in range(N):
            rvals = fmin(Jfunc2d, x0, ftol=0.001, maxiter=niter, disp=0, full_output=1, \
                     args=(L21[i],E21[i],tau21[i],sig21,L31[i],E31[i],tau31[i],sig31))
            x = rvals[0]
            iter = rvals[2]
            if iter < niter:
                Tf[i] = x[0]
                p[i] = 100. * x[1] / P_SCALE # units is %

#       Quality control
#       ---------------
        m = isnan(Tf) == False
        m = m & (Tf<1800.)
        m = m & (p>0) & (p<=100)

#       Add solution as attributes
#       --------------------------
        self.m = m
        self.Tf = Tf
        self.p = p

#       Replace fire size with median size for those fires that did not converge
#       ------------------------------------------------------------------------
        I = (m == False)
        self.p[I] = median(self.p[m])

        self.farea = (self.p/100.) * self.pixar    # km2
        self.hflux = 0.001 * self.pow / self.farea # kW/m2

#       Print out results
#       -----------------
        y = 100. * ( Tf[m].size ) / N + 0.05
        if Verbose:
            print_stats('__header__','Classic Dozier - Variational Results (Yield: %4.1f%%)'%y)
            print_stats('Tf (K)',Tf[m])
            print_stats('p (%)',p[m])
            print_stats('A (km2)',self.farea[m])
            print_stats('HF (kW/m2)',self.hflux[m])
            print_stats('__footer__')

#       Plot KDE
#       --------
###        if Verbose:
###            plot_dozier(Tf[m],p[m],L21[m],E21[m],tau21[m],L31[m],E31[m],tau31[m],\
###                    'Variational','var',pow=self.pow[m])

#........................................................................................

    classic = classic_var  # alias classic to clasic_var
    dozier = classic_var  # alias classic to clasic_var


    def classic_fp(self,tau21=0.864,tau31=0.864):
        """
        Implements the classic Dozier algorithm. On input, *tau_21* and
        *tau_31* are the atmospheric transmittances at 4 and 11 microns
        (MODIS channels 21 and 31). Fixed point version.
        """


#       Compute radiances
#       -----------------
        L21 = B21(self.T21) 
        L31 = B31(self.T31) 
        E21 = B21(self.Tb21)
        E31 = B31(self.Tb31)

        if isscalar(tau21): tau21 = tau21 * ones(L21.shape)
        if isscalar(tau31): tau31 = tau31 * ones(L31.shape)
            
#       The nonlinear equation to be solved is: 
#               B21(Tf) = a + b * B31(Tf)
#       or
#               T = iB21(a + b * B31(T))
#       ---------------------------------------
        r21 = (L21-E21)/(L31-E31)
        a21 = E21 - r21 * E31 
        b21 = r21 * tau31 / tau21

        if self.verb > 0:
            print_stats('__header__','Classic Fixed-point Dozier - Inputs')
            print_stats('DT21',self.T21-self.Tb21)
            print_stats('DT31',self.T31-self.Tb31)
            print_stats('__sep__')
            print_stats('b21',b21)
            print_stats('a21',a21)
            print_stats('__footer__')

#       Used fixed point algorithm to find solution
#       -------------------------------------------
        Tf = fixed_point(Tfunc21,self.T21,xtol=0.001,args=(a21,b21))
        p = 100. * (L21 - E21) / ( tau21 * B21(Tf) - E21 ) 

#       Quality control
#       ---------------
        m = isnan(Tf) == False
        m = m & (p>0)

#       Add solution as attributes
#       --------------------------
        self.m = m
        self.Tf = Tf
        self.p = p

#       Print out results
#       -----------------
        y = 100. * ( Tf[m].size ) / a21.size + 0.05
        print_stats('__header__','Classic Dozier - Fixed-point Results (yield: %4.1f%%)'%y)
        print_stats('Tf (K)',Tf[m])
        print_stats('p (%)',p[m])
        print_stats('__footer__')

#       Plot KDE
#       --------
        plot_dozier(Tf[m],p[m],L21[m],E21[m],tau21[m],L31[m],E31[m],tau31[m],\
                    'Fixed-point','fp',pow=self.pow[m])

 #........................................................................

    def bimodal_u(self,tau21=0.864,tau31=0.864,Verbose=False):
        """
        Implements the bi-modal Dozier algorithm. On input, *tau_21* and
        *tau_31* are the atmospheric transmittances at 4 and 11 microns
        (MODIS channels 21 and 31). 

        This is the unconstrained version, meaning that no additional
        MODIS channels are used. Instead, the most likely value of these
        fire properties are returned:

        a_F  ---  fire area (m2)
        h_F  ---  fire heat flux (kW/m2)
        r_F  ---  fraction of flamming energy
        frp_F  ---  flaming fre radiative power

        """

        self.algo = 'bimodal'

#       Setup Bayesian parameters
#       -------------------------
        self.bayes_setup(tau21=tau21,tau31=tau31)

#       Reserve space for output
#       ------------------------
        N = self.lon.size
        self.a_F = - 99.99 * ones(N)
        self.h_F = - 99.99 * ones(N)
        self.r_F = - 99.99 * ones(N)
        self.frp_F = - 99.99 * ones(N)
        self.m = zeros(N).astype('boolean')
                
        S21, S31, F21, F31 = (self.S21,self.S31,self.F21,self.F31)

        if Verbose:
            if N>100:
                Np = range(0,N,N/100)
                Np = range(N)
            elif N>10:
                Np = range(0,N,N/10)
            else:
                Np = range(N)
            print ""
            print "      Unconstrained Bimodal Dozier"
            print "      ----------------------------"
            print ""
            print "  %  |    Lon    Lat  b |    r_F     h_F"
            print "     |    deg    deg    |     %     kW/m2" 
            print "---- |  ------ ------ - | -------- --------"

#       Estimate parameters for each fire
#       ---------------------------------
        for n in range(N):

            L21,   L31   = (self.L21[n],   self.L31[n])
            E21,   E31   = (self.E21[n],   self.E31[n])
            tau21, tau31 = (self.tau21[n], self.tau31[n])
            pixar = self.pixar[n]
            pow = self.pow[n]

#           Estimate admissible solutions for (Ts,Tf) in range
#           --------------------------------------------------
            ps, pf, kappa = bayes_single(L21, E21, tau21, 
                                         L31, E31, tau31,
                                         S21, S31, F21, F31 )

#           Parameters in phase space
#           -------------------------
            r_F = pf * F21 / ( pf * F21 + ps * S21 ) # non-dimensional
            a_F = (pf/100.) * pixar                  # km2
            h_F = 0.001 * r_F * pow / a_F            # kW/m2
            pow_F = r_F * pow                        # MW
            
#           Kernel density estimates
#           ------------------------
            i = ((ps+pf)>=0)          # quality control
            if any(i):
                self.m[n]     = True
                self.a_F[n]   = mle_kde(a_F[i]) * 1e6 # m2
                self.r_F[n]   = mle_kde(r_F[i])       # %
                self.h_F[n]   = mle_kde(h_F[i])       # kW/m2
                self.pow_F[n] = mle_kde(frp_F[i])     # MW

            if Verbose:
                if n in Np:
                    ip = int(0.5+100.*n/N)
                    print "%3d%% | %7.2f %6.2f   | %8.2f %8.2f"%\
                         (ip,self.lon[n],self.lat[n], \
                          self.r_F[n],self.h_F[n])

 #........................................................................

    def bayes_setup(self,tau21=0.864,tau31=0.864, grid_type='T',
                    srange=[350.,650.],frange=[650.,1800.]):
        """
        Implements a "bayesian" version of the  Dozier algorithm. On input, 
          tau_21, tau_31  --  atmospheric transmittances at 4 and 11 microns
                              (MODIS channels 21 and 31).
          srange          --  Range of temperatures (K) for SMOLDERING fires
          frange          --  Range of temperatures (K) for FLAMING   fires
        This is a really simple minded, brute force algorithm for now.

        Approach: once all randiances are specified (including the
        fire ones), then we can solve for the smoldering/flaming area
        fractions ps/pf:
           
          DS21 * ps + DF21 * pf = L21 - E21
          DS31 * ps + DF31 * pf = L31 - E31

        and,

          DS21 = tau21 * B21(Ts) - E21
          DF21 = tau21 * B21(Tf) - E21
          etc.

        If *grid_type* is 'T' the 2x2 phase space will have increments constant
        in temperature, other it will be  constant in radiances.

        """


#       Compute radiances
#       -----------------
        self.L21 = B21(self.T21) 
        self.L31 = B31(self.T31) 
        self.E21 = B21(self.Tb21)
        self.E31 = B31(self.Tb31)

        self.tau21 = tau21
        self.tau31 = tau31
        if isscalar(tau21): self.tau21 = tau21 * ones(self.L21.shape)
        if isscalar(tau31): self.tau31 = tau31 * ones(self.L31.shape)

        if len(srange)==3: Ns = srange[2]
        else:              Ns = 250
        if len(frange)==3: Nf = frange[2]
        else:              Nf = 250

#       The grid is uniform Ts, Tf
#       --------------------------
        if grid_type == 'T':
            ds, df = ( (srange[1]-srange[0])/Ns, (frange[1]-frange[0])/Nf )
            self.Ts,  self.Tf  = mgrid[srange[0]:srange[1]:ds,frange[0]:frange[1]:df]
            self.S21, self.F21 = (B21(self.Ts),B21(self.Tf))

#       The grid is uniform in B21/F21
#       ------------------------------
        else:
            s21 = ( B21(srange[0]), B21(srange[1]) )
            f21 = ( B21(frange[0]), B21(frange[1]) )
            ds, df = ( (s21[1]-s21[0])/Ns, (f21[1]-f21[0])/Nf )
            self.S21, self.F21 = mgrid[s21[0]:s21[1]:ds,f21[0]:f21[1]:df]
            self.Ts,  self.Tf  = (iB21(self.S21), iB21(self.F21))

        self.S31, self.F31 = (B31(self.Ts),B31(self.Tf))

#       Evaludate radiances on the phase-space grid: smoldering vs. flaming temps
#       -------------------------------------------------------------------------
        if self.verb > 0:
            print_stats('__header__','Bayesian Dozier - Inputs')
            print_stats('L21',self.L21)
            print_stats('L31',self.L31)
            print_stats('__sep__')
            print_stats('E21',self.E21)
            print_stats('E31',self.E31)
            print_stats('__sep__')
            print_stats('Ts',iB21(self.S21))
            print_stats('Tf',iB21(self.F21))
            print_stats('__sep__')
            print_stats('S21',self.S21)
            print_stats('S31',self.S31)
            print_stats('__sep__')
            print_stats('F21/100',self.F21/100)
            print_stats('F31/100',self.F31/100)
            print_stats('__footer__')
            
        return

    def bayes_one(self,n):
        """
        Runs bayes_single() for a given observation with index "n"; you
        must call method bayesian() first. This is a convenience, short-hand
        method for development purposes only.
        """
        return bayes_single(self.L21[n],self.E21[n],self.tau21[n],
                            self.L31[n],self.E31[n],self.tau31[n],
                            self.S21,self.S31,self.F21,self.F31,Verb=True)
        

    def design(self,n,i,j,ka,kb=None,vmin=0.99,vmax=1.):
        """
        For an observation with index "n", and 2 additional walengths
        ka and kb (in microns), computes the likelihood function
        in "phase space" (meaning, the range of smoldering/flaming
        temperatures). The input (i,j) are used to create the synthetic
        observations corresponding to Ts[i,j] and Tf[i,j] at k1 and k2. 
        Requires setup by the bayes_setup() method.
        """
        
#       Find corresponding ps, pf in phase space (Ts,Tf)
#       ------------------------------------------------
        ps, pf, kappa = self.bayes_one(n)
        ps, pf = (ps/100.,pf/100.)
        Ts, Tf = (self.Ts, self.Tf)

#       Get Background radiances based on 4/11 micron Tb
#       ------------------------------------------------
        Tb_ = (self.Tb21[n]+self.Tb31[n])/2. # "b" for background here
        Ea_ = planck(-ka,Tb_)

#       Compute La in all of the phase space
#       ------------------------------------
        La = ps * planck(-ka,Ts) + pf * planck(-ka,Tf) + (1-ps-pf) * Ea_
 
#       The notional ground truth
#       -------------------------
        La_ = La[i,j]

#       Ok, now use Bayes theorem and compute the likelihood function
#       -------------------------------------------------------------
        siga = 0.1 * La_
        va = (La-La_)/siga

#       Do the same for second channel
#       ------------------------------
        if kb != None:
            Eb_ = planck(-kb,Tb_)
            Lb = ps * planck(-kb,Ts) + pf * planck(-kb,Tf) + (1-ps-pf) * Eb_ 
            Lb_ = Lb[i,j]
            sigb = siga
            vb = (Lb-Lb_)/sigb

#       Evaluate likelihood
#       -------------------
        if kb==None:
            P = exp( -va*va/2.0 ) # / (siga * sqrt(2*pi))
        else:
            P = exp( -(va*va + vb*vb)/2.0 ) # / ( siga * sigb * 2. * pi )

#       Plot it
#       -------
        clf()
        pcolor(Tf,Ts,P,vmin=0.99,vmax=1.)
        plot([self.Tf[i,j],],[self.Ts[i,j],],'wo')
        colorbar()
        xlabel('Flaming Temperature (K)')
        ylabel('Smoldering Temperature (K)')
        if kb==None:
            title('Likelihood given (3.959,11.03,%s) $\mu$m'%ka)
        else:
            title('Likelihood given (3.959,11.03,%s,%s) $\mu$m'%(ka,kb))

        return P
                 
#............................................................................

def bayes_single(L21,E21,tau21,L31,E31,tau31,S21,S31,F21,F31,Verb=False):
    """
    Given a *single* measurement of pixel radiances (L21,L31) background
    randiances (E21,E31) and atmospheric transmittances (tau21,Tau31)
    for MODIS channels 21 and 31 (4 and 11 microns) it evaluates the
    fractional smoldering/flaming areas corresponding to each
    smoldering/flamming radiance given by (S21,S31)/(F21,F31),
    respectifully. 

    Input
    -----
       Scalars: L21, E21, tau21, L31, E31, tau31
       Arrays:  S21, S31, F21, F31

    Output
    ------
       ps    --- fractional smoldering area
       pf    --- fractional flaming area
       kappa --- condition number

    Recall that it is not always possible to find a solution (ps,pf)
    for all possible input combinations. At those points where a
    physical solution was not possible, the fractional areas have been
    set to -1. In addition, the condition number of the 2x2 matrix
    used to compute (ps,pf) is also returned as a reliability 
    indicator of the *numerical* solution.

    """

#   Matrix elements
#   ---------------
    shape = S21.shape
    DS21 = tau21 * S21.ravel() - E21
    DF21 = tau21 * F21.ravel() - E21
    DS31 = tau31 * S31.ravel() - E31
    DF31 = tau31 * F31.ravel() - E31

#   RHS
#   ---
    DL21 = L21 - E21
    DL31 = L31 - E31

#   Determinant and condition number
#   --------------------------------
    det = DS21 * DF31 - DS31 * DF21
    chi = DS21 + DF31
    sqd = sqrt(chi*chi - 4 * det)
    kappa = abs((sqd + chi) / ( sqd - chi ))  # condition number
    kappa = where(kappa<1.,1./kappa,kappa)    # ensure kappa > 1 (just in case)

#   Solutions
#   ---------
    ps = (DF31 * DL21 - DF21 * DL31) / det
    pf = (DS21 * DL31 - DS31 * DL21) / det

#   Quality control
#   ---------------
    m = (isnan(ps) | isnan(pf) | isinf(ps) | isinf(pf) | (ps<0.) | (pf<0.) | (ps>1.) | (pf>1.) )
    n = (m==False)
    ps[m] = -1.
    pf[m] = -1.

    if Verb:
        y = 100. * ps[n].size / ps.size
        pr = 100 * pf / ( ps + pf )
        qf = pf * F21.ravel() 
        qs = ps * S21.ravel()
        qr = 100 * qf / (qf + qs)
        
        print_stats('__header__','Bayesian Dozier - Results (Yield: %4.1f%%)'%y)
        print_stats('ps',100*ps[n])
        print_stats('pf',100*pf[n])
        print_stats('pt',100*(pf[n]+ps[n]))
        print_stats('__sep__')
        print_stats('ps*S21',qs[n])
        print_stats('pf*F21',qf[n])
        print_stats('qs+qf',qs[n]+qf[n])
        print_stats('__sep__')
        print_stats('pr',pr[n])
        print_stats('qr',qr[n])
        print_stats('kappa',kappa[n])
        print_stats('__footer__')

#   Reshape
#   -------
    ps = 100 * reshape(ps,shape)
    pf = 100 * reshape(pf,shape)
    kappa = reshape(kappa,shape)

#   Return fractional area and condition number
#   -------------------------------------------
    return (ps,pf,kappa)

def plot_dozier(Tf,p,L21,E21,tau21,L31,E31,tau31,algo,prefix,pow=None):

        subplot(211)
        plot_kde(Tf,300.,1800.,500,'Fire Kinectic Temperature (K) - %s'%algo)
        subplot(212)
        plot_kde(p,0.,5.,500,'Fire Fractional Area (%)')
        savefig('%s.tf_kde.png'%prefix)
        subplot(111)

        p = p / 100 # units of fraction

        clf()
        L21_ = tau21*p*B21(Tf)+(1-p)*E21
        L31_ = tau31*p*B31(Tf)+(1-p)*E31
        plot(L21,L21_,'bo',L31,L31_,'ro')
        title('L21 (Blue) --- L31 (Red) --- %s'%algo)
        xlabel('Observed'), ylabel('Fitted')
        savefig('%s.fits.png'%prefix)

        if pow != None:
            clf()
            plot(pow,p*B21(Tf),'o')
            title('FRP vs p * B$_{21}$(T$_f$) - %s'%algo)
            xlabel('Fire Radiative Power')
            ylabel('p * B$_{21}$(T$_f$)')
            savefig('%s.pow_pB21.png'%prefix)

#........................................................................

def mle_kde(X,N=32):
    """
    Uses a Kernel Density Estimate to return to most likely value of X.
    """
    X = X.ravel()
    bins = linspace(X.min(),X.max(),N)
    kernel = kde.gaussian_kde(X)
    pdf = kernel(bins)
    j = pdf.argmax()
    return bins[j]

#........................................................................

def Jfunc(T,L21,E21,tau21,L31,E31,tau31):
    """Forces the 2 areas to be equal. Does not work too well (low yields)"""
    p21 = (L21 - E21) / ( tau21 * B21(T) - E21 ) 
    p31 = (L31 - E31) / ( tau31 * B31(T) - E31 ) 
    d = p21 - p31
    return d*d

def Jfunc2d(x,L21,E21,tau21,sig21,L31,E31,tau31,sig31):
    """A sum of 2 J_o kind of terms; could de-emphasize chanel 31"""
    Tf = x[0]
    #p = x[1] / 100000.
    p = x[1] / P_SCALE  # normalize so that control variables have same order of magnitude
    v21 = (L21 - ( p * tau21 * B21(Tf) + (1-p) * E21)) / sig21
    v31 = (L31 - ( p * tau31 * B31(Tf) + (1-p) * E31)) / sig31
    return v21 * v21 + v31 * v31

def Tfunc21(T,a,b):
    return iB21(a + b * B31(T))

def Tfunc31(T,a,b):
    return iB31(a + b * B21(T))

def fixed_point(func, x0, args=(), xtol=1e-4, maxiter=50):
    """Find the point where func(x) == x
    
    Given a function of one or more variables and a starting point, find a
    fixed-point of the function: i.e. where func(x)=x.

    Uses Steffensen's Method using Aitken's Del^2 convergence acceleration.
    See Burden, Faires, "Numerical Analysis", 5th edition, pg. 80

    This is a customized version of a SciPy function.

    """

    x0 = asarray(x0)
    p0 = x0                  
    for iter in range(maxiter):
        p1 = func(p0, *args)
        p2 = func(p1, *args)
        d = p2 - 2.0 * p1 + p0
        p = where(d == 0, p2, p0 - (p1 - p0)*(p1-p0) / d)
        relerr = where(p0 == 0, p, (p-p0)/p0)
        if all(relerr < xtol):
            return p
        p0 = p
        
#    print "Failed to converge after %d iterations, value is %s" % (maxiter,p)
    return p

#............................................................................

def utClassic():
    fires = DOZIER('data/182',Verb=1,qc_thresh=50.)
    fires.classic_fp()
    fires.classic_var()

def utBayesian():
    fires = DOZIER('182',Verb=1,qc_thresh=50.)
    fires.bayes_setup()
    return fires

def utBimodal():
    import VegType
    fires = DOZIER('182',Verb=1,qc_thresh=50.)
    fires.veg = VegType.getSimpleVeg(fires.lon,fires.lat,Path=igbp_dir)
    fires.bimodal_u(Verbose=True)
    return fires

def utDesign(tau21=0.864,tau31=0.864):
    f = DOZIER('data/182',Verb=1,qc_thresh=50.)
    f.bayes_setup(tau21=1.,tau31=1.);
    P = f.design(345,120,120,3.75,8.55)
    return (f,P)

def utAttach():
    f = DOZIER('data/182',Verb=1,qc_thresh=50.)
#    f.attach('http://thing4.gsfc.nasa.gov:9090/dods/GEOS-5/ARCTAS/0.5_deg/assim/tavg3d_dyn_v',Vars=('t','qv','o3','delp','ps'))
    f.attach('http://thing4.gsfc.nasa.gov:9090/dods/GEOS-5/ARCTAS/0.5_deg/assim/tavg3d_dyn_v',Vars=('ps'))


if __name__ == "__main__":
    fires = utBimodal()


