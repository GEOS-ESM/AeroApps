#
# Playing around with ML tuning scripts.
#

from numpy import pi, savez
import mltune as ml

def laod_estimate(expid,channel,nymd_beg,nymd_end,inst,BBOX=None,N=500,qch=None):

    yymm = nymd_beg / 100
    odsfn = '%s.%s-%d.obs.%d.ods'%(expid,inst,channel,yymm)
    if qch is None:
        npzfn = '%s.%s-%d.%d.npz'%(expid,inst,channel,yymm)
    else:
        npzfn = '%s-%d.%s-%d.%d.npz'%(expid,qch,inst,channel,yymm)

#   Date and time ranges
#   --------------------
    NYMD = range(nymd_beg,nymd_end+1)
    NHMS = range(0,240000,30000)

#   Regional box
#   ------------

#   Latitude bands
#   --------------
    if BBOX is None:
        BBOX = [ (-180.,-90., 180., -30.), 
                 (-180.,-30., 180., -15.), 
                 (-180.,-15., 180.,   0.),
                 (-180.,  0., 180.,  15.),
                 (-180., 15., 180.,  30.),
                 (-180., 30., 180.,  90.)  ]

#   Error estimates for each bounding box
#   -------------------------------------
    sigO = 0.13
    sigF = 0.40
    L = 145.0 / 1000.
    alpha = ( sigO, sigF, L )
    Alpha = ()
    sigAlpha = ()
    for BBox in BBOX:

        V = ml.getODSts(odsfn,NYMD,NHMS,N=N,BBox=BBox,qch=qch)

        alpha0 = alpha
        alpha, sigalpha = ml.estimate(alpha0,V,ml.cov_winplaw,xtol=0.01)
 
        Alpha = Alpha + (alpha,)
        sigAlpha = sigAlpha + (sigalpha,)

#   Save results
#   ------------
    savez(npzfn, expid=expid, channel=channel, inst=inst,
                 BBOX=BBOX, Alpha=Alpha, sigAlpha=sigAlpha )

#   Summary of results
#   ------------------
    k = 0
    alphan = ( 'sigO', 'sigF', '   L')
    print " "
    print "Input ODS File: ", odsfn
    for BBox in BBOX:
        alpha, sigalpha = (Alpha[k], sigAlpha[k])
        print "Latitude Band: %d %d"%(BBox[1],BBox[3])
        for j in range(alpha.size):
            print "    %s = %7.3f +/- %4.3f"%(alphan[j], alpha[j], sigalpha[j])
        k = k + 1

    return Alpha, sigAlpha

#........................................................................

if __name__ == "__main1__":

    expid = 'a0005'
    channel = 550
    nymd_beg = 20080601
    nymd_end = 20080630

#   AERONET
#   -------
    BBOX = [ (-140.,  15., -60., 60.),   # North America
             ( -10., -35.,  40., 65.),   # Europe
             (  40.,  20., 140., 65.) ]  # Asia

    BBOX = [ (-180.,  -90., 180., 90.), ] # North America


    for inst in ( 'misr', 'modo', 'mydo' ):
        for channel in ( 870, 660, 470 ):
            try:
                Alpha, sigAlpha = laod_estimate(expid,channel,nymd_beg,nymd_end,inst)
            except:
                print ">>>>>>> Something wrong with "+inst

#....................................................

if __name__ == "__main__":

    expid = 'a0005'
    channel = 550
    nymd_beg = 20080601
    nymd_end = 20080630

    for inst in ( 'mydo', ):
        for qch in ( 1, 3 ):
            try:
                Alpha, sigAlpha = laod_estimate(expid,channel,nymd_beg,
                                                nymd_end,inst,qch=qch)
            except:
                print ">>>>>>> Something wrong with "+inst
