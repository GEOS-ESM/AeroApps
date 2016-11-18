"""
    Calculates the bidirecitonal surface reflectance given
    RTLS kernel weights
"""

from   numpy import pi, tan, cos, sin
from   numpy import arctan, sqrt, arccos
import numpy as np

def rtlsReflectance(Kiso,Kgeo,Kvol,sza,vza,saa,vaa):
    """ 
    Lucht et al. (2000) IEEE Transactions on Geoscience and Remote Sensing
    """


    h_b  = 2.0
    b_r  = 1.0

    raa  = vaa - saa  # raa = 0 is back scattering direction

    # degrees to radians
    raa = raa*pi/180.0
    sza = sza*pi/180.0
    vza = vza*pi/180.0

    Riso = 1.0

    # Calculate Rgeo
    szaP        = arctan(b_r*tan(sza))
    vzaP        = arctan(b_r*tan(vza))
    cos_thetaP  = cos(szaP)*cos(vzaP) + sin(szaP)*sin(vzaP)*cos(raa)
    dP          = sqrt(tan(vzaP)*tan(vzaP) + tan(szaP)*tan(szaP) - 2.0*tan(vzaP)*tan(szaP)*cos(raa))
    amfP        = (1./cos(vzaP)) + (1./cos(szaP))
    temp        = h_b*sqrt(dP*dP + (tan(szaP)*tan(vzaP)*sin(raa))*(tan(szaP)*tan(vzaP)*sin(raa)))/amfP
    cos_t       = np.array([min([1,t]) for t in temp])       
    t           = arccos(cos_t)
    sin_t       = sin(t)
    O           = (1.0-(1./pi)*(t-sin_t*cos_t))*amfP

    Rgeo        = 0.5*(1+cos_thetaP)*(1./cos(vzaP))*(1./cos(szaP)) - O

    # Calculate Rvol
    cos_theta   = cos(sza)*cos(vza) + sin(sza)*sin(vza)*cos(raa)
    theta       = arccos(cos_theta)

    Rvol        = (((0.5*pi - theta)*cos_theta + sin(theta))/(cos(sza) + cos(vza))) - 0.25*pi

    # Calculate bi-directional reflectance
    rho         = Kiso*Riso + Kgeo*Rgeo + Kvol*Rvol

    rho[Kiso < -900] = -999.

    return rho






