#!/usr/bin/env python
"""
 Kernel density estimate.

"""

import sys    
from pylab import *
from scipy import stats, mgrid, c_, reshape, random, rot90

from scipy.stats.stats import linregress

from numpy.lib.io import save, load

from grads import GrADS

#                     G L O B A L S
                                 
#   Hardwire these for now...
#   -------------------------

year = 2008
month = 06
yymm = year*100 + month

#   For titles
#   ----------
qname = {}
qname['modo'] = 'MODIS/TERRA Ocean'
qname['modl'] = 'MODIS/TERRA Land'
qname['mydo'] = 'MODIS/AQUA  Ocean'
qname['mydl'] = 'MODIS/AQUA  Land'
qname['deep'] = 'MODIS/AQUA  Deep-Blue'
qname['misr'] = 'MISR'
qname['anet'] = 'AERONET'
qname['parl'] = 'PARASOL Land'
qname['paro'] = 'PARASOL Ocean'
qname['omi']  = 'OMI'
    
#-----------------------------------------------------------------

def cat (X, Y):
    """
    Given 2 arrays of same shape, returns array of shape (2,N),
    where N = X.size = Y.size
    """
    xy = concatenate((X.ravel(),Y.ravel())) # shape is (N+N)
    return reshape(xy,(2,X.size))         # shape is (2,N)

def calculate(var,lev,units,bins=None):

#   Start grads and open relevant files
#   -----------------------------------
    ga = GrADS(Echo=False,Window=False)
    ga("sdfopen a0005.gritas_aod.omf.%s.nc"%yymm)
    ga("sdfopen a0005.gritas_aod.oma.%s.nc"%yymm)
    ga('set lev '+str(lev))

#   Read data, flatten it and remove undefs
#   ---------------------------------------
    print "Reading <%s> O-F..."%var
    omf = ga.expr("aod_"+var+'.1').data.ravel()
    print "Reading <%s> O-A..."%var
    oma = ga.expr("aod_"+var+'.2').data.ravel()
    i = find( (abs(omf)<1e3) & (abs(oma)<1e3) )
    omf = omf[i]
    amf = omf - oma[i]
    del oma, i, ga
    basename=var+'-'+str(lev)+'.npy'
    save('omf.'+basename,omf)
    save('amf.'+basename,amf)

#   Multivariate: omf in x, amf in y
#   --------------------------------
    x_values = omf
    y_values = amf
    if bins == None:
        bins = arange(-4., 4., 0.1 )

    x_bins = bins
    y_bins = bins

    N  = x_values.shape
    Nx = x_bins.size
    Ny = y_bins.size

    do_1d = True
    do_2d = True

    if N==0:
        print "Nothing to do for "+var
        return

#   1-D calculation
#   ---------------
    if do_1d:
        print "Starting 1D kernel density estimation with %d points..."%N
        figure()
        kernel = stats.kde.gaussian_kde(x_values)
        z = kernel(x_bins)
        plot(x_bins,z,label='O-F')
        kernel = stats.kde.gaussian_kde(y_values)
        z = kernel(y_bins)
        plot(y_bins,z,label='A-F')

        title('%dnm Log(eps+AOD) - %s [%d-%d]'%(lev,qname[var],year,month))
        xlabel(units)
        ylabel('PDF')
        legend()
        grid(b=True)
        savefig('pdf_1d-'+var+'-'+str(lev)+'.png')

#   2D calculation
#   --------------
    if do_2d:
        print "Starting the 2D kernel density estimation with %d points..."%N
        kernel = stats.kde.gaussian_kde(cat(x_values,y_values))

        print "Evaluating 2D kernel on grid with (Nx,Ny)=(%d,%d) ..."%(Nx,Ny)
        X, Y = meshgrid(x_bins,y_bins) # each has shape (Ny,Nx)
        Z = kernel(cat(X,Y))           # shape is (Ny*Nx)
        Z = reshape(Z,X.shape)
        
        print "X, Y, Z shapes: ", X.shape, Y.shape, Z.shape

        #   Save to file
        #   ------------
        fname = 'kde_2d-'+var+'-'+str(lev)+'.npy'
        print "Saving to file <"+fname+"> ..."
        save(fname,Z)

    return N
 
#---       
def plot_kde(var,lev,bins=None):

#   Load 2D KDE
#   -----------
    basename = 'kde_2d-'+var+'-'+str(lev)
    P = load(basename+'.npy')

    if bins == None:
        bins = arange(-2., 2., 0.05 )

    x = bins
    y = bins
   

    #   Plot results with 2 different colormaps
    #   ---------------------------------------
    print "Plotting..."

    fig = figure()
    ax = fig.add_axes([0.1,0.1,0.75,0.75])
    imshow(P, cmap=cm.gist_earth_r, origin='lower', 
           extent=(x[0],x[-1],y[0],y[-1]) )
    xlabel('O-F')
    ylabel('A-F')
###    title('%dnm Log(eps+AOD) - %s [%d-%d]'%(lev,qname[var],year,month))
    title('PDF - %dnm Log(eps+AOD) - %s'%(lev,qname[var]))
    grid()

#   Centroid
#   --------
    X, Y = meshgrid(x,y) # each has shape (Ny,Nx)
    yy = sum(P*Y,axis=0)/sum(P,axis=0)
###    plot(x,yy,'r',label='Centroid')

#   Regression line
#   ---------------
    basename2=var+'-'+str(lev)+'.npy'
    omf = load('omf.'+basename2)
    amf = load('amf.'+basename2)
    slope, intercept, r, prob, see = linregress(omf,amf)
    yy = slope * x + intercept
###    plot(x,yy,'k',label='Regression')
###    legend(loc='upper left')

#   Tighter colorbar
#   ----------------
    bbox = ax.get_position()
    l,b,w,h = bbox.bounds # mpl >= 0.98
    cax = axes([l+w+0.02, b, 0.04, h]) # setup colorbar axes.
    colorbar(cax=cax) # draw colorbar

    savefig(basename+'.png')
    
#................................................................

if __name__ == '__main__':
   
    bins = arange(-0.6, 0.6, 0.01 )
    var = 'modo'
    units = 'nm'
    lev = 470

    for var in ( 'modo', 'modl', 'mydo', 'mydl', 'misr', 'deep', 'parl', 'paro', 'omi' ):

        try:
            N = calculate(var,lev,units,bins=bins)
            if N>0:
                plot_kde(var,lev,bins=bins)
        except:
            print "Problems with "+var+", ignoring it"

    print "All done."
