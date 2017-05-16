"""
Explore Cholesky decomp of correlation functions.
"""
from mltune import *
from numpy import linspace, meshgrid

from matplotlib.pyplot import *

if __name__ == "__main__":

    # Create regular lat/lon grid and flatten it
    # ------------------------------------------
    x = linspace(-20,5,25)
    y = linspace(-20,-5,20)
    X, Y = meshgrid(x,y)
    lon = X.ravel()
    lat = Y.ravel()

    # Build Correlation matrix
    # ------------------------
    R = rdist(lon,lat)
    C = powerlaw(R,L=100E3)

    # Cholesky decomposition
    # ----------------------
    L = cholesky(C) # Lower triangular
    U = L.T         # Upper triangular

    # first row of you is the correlation function
    # --------------------------------------------    
    figure()
    plot(R[0]/1000,U[0]/U[0,0],'o',label='U First Row') 
    plot(R[0]/1000,C[0],'ro',label='Correlation Function')
    xlabel('R [km]')
    title('C = U.T * U')
    grid()
    legend()
    savefig('first_row.png')

    # Subsequent rows of U are shifted
    # --------------------------------
    figure()
    for i in range(99,500,100):
        plot(U[i,i:]/U[i].max(),label='Row %d'%i)
    title('U Rows')
    grid()
    legend()
    savefig('U_rows1.png')
    
    figure()
    for i in range(0,10):
        plot(U[i,i:i+20]/U[i].max(),label='Row %d'%i)
    title('U Rows')
    grid()
    legend()
    savefig('U_rows2.png')



    
