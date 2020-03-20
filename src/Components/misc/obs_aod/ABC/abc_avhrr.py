"""

Playing with AVHRR.

"""

import os

from numpy import pi, sin, cos, log, zeros, ones, savez, arange, exp
from numpy import c_ as cat, random

from pyobs.sknet       import loadnet

from matplotlib.pyplot import title, savefig
from pyobs.avhrr       import AVHRR_L2B
from abc_modis         import NN, _plotKDE, aodFormat

MISSING = -1.0e20
d2r = pi / 180.

class ABC_AVHRR(NN,AVHRR_L2B):

    def __init__(self, Path='/nobackup/AVHRR/Level2/NPZ/2008/*.npz',
                 coxmunk_lut='/nobackup/NNR/Misc/avhrr.cox-munk_lut.npz',
                 N_bal=50*1000,Input=None,Target=None,Verb=False):

        self.verbose = Verb
        self.ident = 'avhrr'
        self.surface = 'ocean'
        
        AVHRR_L2B.__init__(self,Path,Verb=Verb)

        # Balance Observing system
        # ------------------------
        if N_bal is not None:
            I = self.balance(N_bal)
            self.reduce(I)

        # Define wind speed dependent ocean albedo
        # ----------------------------------------
        self.Wind = self.wind
        self.getCoxMunk(coxmunk_lut)

        # Air mass factor
        # ---------------
        self.amf = (1./cos(d2r*self.SolarZenith))+(1./cos(d2r*self.SensorZenith))  

        # Glint angle
        # -----------
        RelativeAzimuth = self.SensorAzimuth # = anchor_relative_azimuth
        cosGlintAngle = cos(self.SolarZenith*d2r) * cos(self.SensorZenith*d2r) + \
                        sin(self.SolarZenith*d2r) * sin(self.SensorZenith*d2r) * \
                        cos(RelativeAzimuth*d2r)

        # Angle transforms: for NN calculations we work with cosine of angles
        # -------------------------------------------------------------------
        self.SensorAzimuth = cos(self.SensorAzimuth*d2r)   
        self.SensorZenith  = cos(self.SensorZenith*d2r)
        self.SolarAzimuth  = cos(self.SolarAzimuth*d2r)
        self.SolarZenith   = cos(self.SolarZenith*d2r)
        self.GlintAngle    = cosGlintAngle

        # Sanity check
        # ------------
        self.iValid = (self.tau_630  > -0.01) &\
                      (self.tau_550  > -0.01) &\
                      (self.ref_630 >  0)    &\
                      (self.ref_860 >  0)

#                     (self.tau_860  > -0.01) &\

        # Log transforms
        # --------------
        self.ltau_550 = log(self.tau_550+0.01)
        self.ltau_630 = log(self.tau_630+0.01)
        self.ltau_860 = log(self.tau_860+0.01)
        self.lref_630 = log(self.ref_630)
        self.lref_860 = log(self.ref_860)
                      
        # Default input features
        # ----------------------
        if Input==None:
            self.Input = ( 'tau_630', 'tau_860', 'tpw')
        else:
            self.Input = Input

        # Default Target
        # --------------
        self.laod = True
        if Target==None:
            self.Target = ('tau_550',)
        else:
            self.Target = Target
         
#---

__Months__ = ['JAN','FEB','MAR','APR','MAY','JUN','JUL','AUG','SEP','OCT','NOV','DEC']

def dt2gat(t):
    """
    Convert datetime to grads time.
    """
    gat = "%d:%dZ%d%s%d"%(t.hour,t.minute,t.day,__Months__[t.month-1],t.year)
    return gat

#....................................................................................     

def _avhrrNNR(a):

    a.split()

    Input_Raw = [ 'albedo','tpw',
                  'SensorAzimuth',
                  'SensorZenith',
                  'SolarAzimuth',
                  'SolarZenith',
#                 'amf',
                  'lref_630',
                  'lref_860',
                  'fdu', 'fss', 'fcc', 'fsu',
                ]

    Input_Raw2 = [ 'albedo','tpw',
                  'SensorAzimuth',
                  'SensorZenith',
                  'SolarAzimuth',
                  'SolarZenith',
#                  'amf',
                  'ref_630',
                  'ref_860',
                  'fdu', 'fss', 'fcc', 'fsu',
                ]

    Input_Ret = [  'albedo', 'tpw',
#                  'SensorAzimuth',
#                  'SensorZenith',
#                  'SolarAzimuth',
#                  'SolarZenith',
                  'tau_630',
                  'tau_860',
                  'fdu', 'fss', 'fcc', 'fsu',
                ]

#    Experiments:
#    1  - retrieval correction
#    1b - reflectance based, NPZ files, balanced (N_bal=200K), speciation
#    1c - reflectance based, NPZ files, strongly typed
#         (a.fdu>0.77)|(a.fss>0.87)|(a.fcc>0.755)|(a.fsu>0.8)
#    1d - as 1b but no speciation
#
    ident = 'avhrr'
    expid = 'nnr_001e'
    Target = ['tau_550',]
    Input = Input_Raw
    a.laod = True

    # Training
    # --------
    nHidden = len(Input)
    topology = (len(Input), nHidden, 4, len(Target))
#   topology = (len(Input), nHidden, len(Target))
#    topology = (len(Input), 4, len(Target))
    biases = True
        
    print " "
    print "        AOD Neural Net Retrieval"
    print "        ------------------------"
    print " "
    print "  No. Valid Data:  ", len(a.lon[a.iValid]),  \
                                 int(100.*len(a.lon[a.iValid])/len(a.lon)),'%'
    print " No. Hidden Nodes: ", nHidden
    print "         Topology: ", topology
    print "   Input Features: ", Input[:]
    print "           Target: ", Target[:]
    print " "

    a.train(Input=Input,Target=Target,nHidden=nHidden,
            topology=topology,biases=biases,nproc=8)
    out, reg = a.test()

    a.savenet(expid+"."+ident+'_Tau.net')

    return (expid, ident, Target)

#---
def doPlots(a,expid,ident, Target, bins=None):

    if bins is not None:
        pass
    elif a.laod:
        bins = arange(-3, 0, 0.05 )
    else:
        bins = arange(0., 0.6, 0.01 )

    # Plot KDE of corrected AOD
    # -------------------------
    a.plotKDE(bins=bins,
              x_label='MERRAero',
              figfile=expid+"."+ident+"_kde-"+Target[0][1:]+"-corrected.png")

    # Plot KDE of uncorrected AOD
    # ---------------------------
    targets = a.getTargets(a.iValid).squeeze()
    if a.laod:
        formatter = aodFormat()
        original = a.ltau_630[a.iValid]
    else:
        formatter = None
        original = a.tau_630[a.iValid]

    _plotKDE(targets,original,
             x_label='MERRAero',
             y_label='Original AVHRR',
             x_bins=bins,
             formatter=formatter)
    if a.laod:
       title("Log("+Target[0][1:]+"+0.01)- "+ident)
    else:
       title(Target[0][1:]+"- "+ident)
    savefig(expid+"."+ident+"_kde-"+Target[0][1:]+'-orig.png')

    # Scatter diagram for testing
    # ---------------------------
    # a.plotScat(figfile=expid+"."+ident+"_scat-"+Target[0][1:]+'.png')

#....................................................................................     

if __name__ == "__main__":

    a = ABC_AVHRR(Path='/nobackup/AVHRR/Level2/NPZ/2008/*.npz',
                  N_bal=200*1000,Verb=True)

    bins = arange(-3, 0, 0.05 )
    #bins = arange(-3, -1.2, 0.02 )
    #bins = arange(0, 0.3, 0.3/128. )
    # a.iValid = a.iValid & ( (a.fdu>0.77)|(a.fss>0.87)|(a.fcc>0.755)|(a.fsu>0.8) )

    expid, ident, Target = _avhrrNNR(a)

    #doPlots(a,expid,ident, Target, bins=bins)

def hold():

    m, expid, ident, Target = _deepNNR()
    
    m = ABC_DEEP(verbose=True)
    a.getNNR()
    a.trainSVC()

    I = a.svcIndex
    J = I & (a.svcEval==1)&(a.albedo>0.15)&(a.ldTau550>-5)

    m, expid, ident, Target = _deepNNR()
    doPlots(m, expid, ident, Target)

    m = ABC_DEEP(verbose=True)
    a.getNNR()
    a.trainSVC()
    I = a.svcIndex
    J = I & (a.svcEval==1)&(a.albedo>0.15)&(a.ldTau550>-5)
