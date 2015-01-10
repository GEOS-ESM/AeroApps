#
# Howard U. Flux Tower at Beltsville
#

from numpy import *
from datetime import datetime, timedelta

from grads import GrADS

if __name__ == "__main__":

    # Read the file
    # -------------
    X = loadtxt('hu.flx.v1.2011-07.dat',skiprows=4)
    X = X.T

    I = (X==-9999.9)
    X[I] = NaN

    year, jday, hour,    u,    v,    w,    H,    LE, ustar = \
    X[0], X[1], X[2], X[5], X[7], X[8], X[9], X[10], X[11]

    time = []
    for y, j, h in zip(year,jday,hour):
        secs = (j * 24 + h) * 60 * 60
        date = datetime(int(y)-1,12,31) + timedelta(seconds=int(secs))
        time.append(date)


    H[abs(H)>1000.] = NaN
    LE[abs(LE)>1000.] = NaN

    H[H<-200] = NaN
    LE[LE<-200] = NaN

    # Interpolate GEOS-5 to Beltsville
    # --------------------------------
    #ga = GrADS(Window=False,Echo=False)
    #ga.open('/home/adasilva/opendap/daq/opendap/assim/tavg1_2d_flx_Nx')
    #lat = array([39.055302,])
    #lon = array([-76.878304,])



