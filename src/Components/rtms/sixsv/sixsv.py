"""

Functions from Python interface to selected functions from "Second
Simulation of Satellite Signal in the Solar Spectrum (6S)" by
E. Vermote et al.

"""

from numpy import savez, zeros, ones, linspace, array

import sixsv_

PACE_channels = [  350, 360, 385, 412, 425, 443, 460, 475, 490, 510, 532, 555, 583,
                   617, 640, 655, 665, 678, 710, 748, 765, 820, 865, 1245, 1640, 2135,]

MODIS_channels = [470, 550, 660, 870, 1200, 1600, 2100]

def ocnAlbedo (speed, channels,
               del_azim=None, salinity=None, pigment=None):
    """
    Given,
    
       speed       ---   speed of wind (in m/s)
       wavelength  ---   wavelength of the computation (in nanometer)

    and optionally,
    
       del_azim    ---   azim. of sun - azim. of wind (in deg.), default: 0
       salinity    ---   salinity (in ppt), default = 35
       pigment     ---   pigment concentration (in mg.m-3), default = 0

    it returns,

      albedo      ---   the spherical albedo of the ocean 

    """

    N = len(speed)
    Nc = len(channels)
    if del_azim == None:
        del_azim = zeros(N)

    if salinity == None:
        salinity = 35 * ones(N)

    if pigment == None:
        pigment = zeros(N)

    
    albedo = zeros((Nc,N))
    for i in range(Nc):
        print(channels[i])
        albedo[i] = sixsv_.ocnalbedo(speed, channels[i], del_azim, salinity, pigment)
    
    return albedo.T # (nobs.nch)

#-------------------------------------------------------------------------------
if __name__ == "__main__":

    channels = array(MODIS_channels[:])

    speed = linspace(0.,50.,51)
    
    albedo = ocnAlbedo(speed,channels)
    albedo[0,:] = albedo[1,:]
    
    savez('cox-monk_lut.npz',channels=channels,speed=speed,albedo=albedo)
