#
# Playing around with ML tuning scripts.
#

from numpy import pi
import mltune as ml

if __name__ == "__main__":

    nymd_beg = 20080601
    nymd_end = 20080617
    fname = 'misr.200806.ods'
    fname = 'modo550.200806.ods'

    NYMD = list(range(nymd_beg,nymd_end+1))

#   Regional box
#   ------------
    NHMS = list(range(0,240000,30000))
    BBox = (-40.,6,-20.,26.)

#   Latitude band for a given orbit
#   -------------------------------
    NHMS = (120000,)
    BBox = (-180.,6,180.,26.)

    V = ml.getODSts(fname,NYMD,NHMS,N=500,BBox=BBox)

    sigO = 0.18
    sigF = 0.40
    L = 145.0 / 1000.
    omega = pi/4.
    alpha = ( omega, L )
    alpha = ( sigO, sigF, L )
    
    S = ml.cov_winplaw(alpha,V[0])
    J = ml.llfun(alpha,ml.cov_winplaw,V)

#    alpha_ = ml.estimate(alpha,V,ml.cov_winplaw,xtol=0.01)
    alpha_ = ml.estimate(alpha,V,ml.cov_winplaw)
