"""
Test subpoint associated with L1 orbit using ISEE as reference.
"""

from pyobs import sgp4
from datetime import datetime, timedelta
from numpy import ones, array, sin, cos, pi, arcsin

def getCenterPoint(tyme,lon0=180.):
    """
    Return central (lon,lat) for ortographic projection associated with a
    L1 orbit. For simplicity a circular orbit around the sun is assumed. 
    """
    
    eps = 23.4 # obliquity in degrees, see http://en.wikipedia.org/wiki/Axial_tilt
     
    P_secs = 365.25 * 24 * 60 * 60 # about 1 year

    # Seconds of the day
    # ------------------
    sod = array([ t.hour*3600. + t.minute*60. + t.second for t in tyme])
    omega = 360. / ( 24. * 60. * 60 ) # earth rotation rate (sec-1)

    # Summer Solstice
    # ---------------
    ts = datetime(tyme[0].year,6,21) # good enough

    # Nomalized time by 1 year, having summer solstice as reference
    # -------------------------------------------------------------
    tn = array([(t-ts).total_seconds() for t in tyme])/P_secs # summer solsice is reference
    
    lon = lon0 - sod * omega # at 0Z sun is at 180E

    lat = (180./pi) * arcsin(sin(pi*eps/180.) * cos(2*pi*tn))
    

    return (lon,lat)
    
if __name__ == "__main__":

    # Create time
    # -----------
    t_beg = datetime(2005,1,1)
    t_end = datetime(2005,12,31)
    dt_secs = 30 * 60
    dt = timedelta(seconds=dt_secs)
    N = int((t_end-t_beg).total_seconds()/dt_secs) 
    tyme = array([t_beg + i * dt for i in range(N)])

    lon, lat = getCenterPoint(tyme)
    
