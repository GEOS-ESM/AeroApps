"""

Neural Net training for MISR.

"""

import os

from numpy import pi, cos, log, zeros, ones, savez, arange, exp
from numpy import c_ as cat
###from abc_modis         import NN, _plotKDE, aodFormat

#from ffnet        import loadnet
###from sknet        import loadnet

from pyobs             import NPZ
from pyobs.mapss       import MAPSS
from pyobs.mcd43gf     import MCD43GF
from pyobs.igbp        import IGBP

from matplotlib.pyplot import title, savefig

MISSING = -1.0e20
d2r = pi / 180.

#class ABC_MISR(NN):
class ABC_MISR(object):

    def __init__(self,
                 npzDir='/nobackup/MAPSS/Collocation/MISR',
                 tol=0.5,
                 Input=None,Target=None,
                 verbose=False):

        self.verbose = verbose
        self.ident = 'misr' 
        
        # Read in NPZ files written by collocation app
        # --------------------------------------------
        self.a = MAPSS(npzDir+'/mapss.anet_misr.??????.npz')
        self.m = MAPSS(npzDir+'/mapss.misr_maod.??????.npz')
        self.g = MAPSS(npzDir+'/mapss.misr_geom.??????.npz')
        self.r = MAPSS(npzDir+'/mapss.misr_mref.??????.npz')

        # Ancillary variables by prepAncillary below
        # ------------------------------------------
        self.x = MAPSS(npzDir+'/mapss.xtra_misr.??????.npz')

        # Inherit coordinates from AERONET file
        # -------------------------------------
        self.lon, self.lat, self.time, self.N = (self.a.lon, self.a.lat, self.a.time, self.a.nobs)

        # Air mass factor
        # ---------------
        self.amf1 = (1./cos(d2r*self.g.SolarZenith))+(1./cos(d2r*self.g.SensorZenith1))
        self.amf2 = (1./cos(d2r*self.g.SolarZenith))+(1./cos(d2r*self.g.SensorZenith2))
        self.amf3 = (1./cos(d2r*self.g.SolarZenith))+(1./cos(d2r*self.g.SensorZenith3))
        self.amf4 = (1./cos(d2r*self.g.SolarZenith))+(1./cos(d2r*self.g.SensorZenith4))
        self.amf5 = (1./cos(d2r*self.g.SolarZenith))+(1./cos(d2r*self.g.SensorZenith5))
        self.amf6 = (1./cos(d2r*self.g.SolarZenith))+(1./cos(d2r*self.g.SensorZenith6))
        self.amf7 = (1./cos(d2r*self.g.SolarZenith))+(1./cos(d2r*self.g.SensorZenith7))
        self.amf8 = (1./cos(d2r*self.g.SolarZenith))+(1./cos(d2r*self.g.SensorZenith8))
        self.amf9 = (1./cos(d2r*self.g.SolarZenith))+(1./cos(d2r*self.g.SensorZenith9))

        # Expose reflectances
        # -------------------
       # self.sRef412 = self.s.sRef412
       # self.sRef470 = self.s.sRef470
       # self.sRef660 = self.s.sRef660
       # self.dRef412 = self.r.dRef412
       # self.dRef470 = self.r.dRef470
       # self.dRef660 = self.r.dRef660
       # self.xRef412 = self.dRef412 - self.sRef412 
       # self.xRef470 = self.dRef470 - self.sRef470 
        #------------------------------------------------------------------
        "There were no analogous readings in the geom, mref or maod files."
        "So commented out Expose reflectances             ~Suniyya Waraich"
        #------------------------------------------------------------------

         
        # Expose AOD
        # ----------
        self.aTau440  = self.a.tau440
        self.aTau550  = self.a.tau550
        #self.dTau412  = self.d.tau550
        self.mTau446  = self.m.tau446
        #self.dTau470  = self.d.tau470
        #self.dTau550  = self.d.tau550
        self.mTau558  = self.m.tau558
        #self.dTau660  = self.d.tau660
        self.mTau672  = self.m.tau672
        #-------------------------------------------------------------------------
	"Changed AOD values. Replaced self.d tau values by available self.m "
	" tau values. tau550, tau470 and tau660 were unavailable. ~Suniyya Waraich "
	#-------------------------------------------------------------------------

        #Sanity Check
        #------------
	self.iValid = (self.a.tau550 >-0.01) &\
                      (self.albedo >0) &\
                      (self.m.tau558>-0.01)
                    # &\
                    # (self.d.qa_flag >0)

	
        #-------------------------------------------------------------------
        "Moved sanity check up here to prevent log errors. ~Suniyya Waraich"
        #-------------------------------------------------------------------        
	
	self.angstrom = -log(self.mTau672/self.mTau446)/log(672./446.)
        self.laTau550 = log(self.a.tau550+0.01)
        self.lmTau558 = log(self.m.tau558+0.01)
        #-------------------------------------------------------
        "Propagated changes made in Expose AOD"
        "by changing dTau470 to mTau446 and dTau660 to mTau 672"
        "and self.1mTau550 to self.1mTau558    ~Suniyya Waraich"
        #-------------------------------------------------------
        
        

        # Angle transforms: for NN calculations we work with cosine of angles
        # -------------------------------------------------------------------
        #print "self.g.ScatteringAngle1:",self.g.ScatteringAngle1
        self.ScatteringAngle1 = cos(self.g.ScatteringAngle1*pi/180.0)
        self.ScatteringAngle2 = cos(self.g.ScatteringAngle2*pi/180.0)
        self.ScatteringAngle3 = cos(self.g.ScatteringAngle3*pi/180.0)
        self.ScatteringAngle4 = cos(self.g.ScatteringAngle4*pi/180.0)
        self.ScatteringAngle5 = cos(self.g.ScatteringAngle5*pi/180.0)
        self.ScatteringAngle6 = cos(self.g.ScatteringAngle6*pi/180.0)
        self.ScatteringAngle7 = cos(self.g.ScatteringAngle7*pi/180.0)
        self.ScatteringAngle8 = cos(self.g.ScatteringAngle8*pi/180.0)
        self.ScatteringAngle9 = cos(self.g.ScatteringAngle9*pi/180.0)

        self.RelativeAzimuth1 = cos(self.g.RelativeAzimuth1*pi/180.0)
        self.RelativeAzimuth2 = cos(self.g.RelativeAzimuth2*pi/180.0)
        self.RelativeAzimuth3 = cos(self.g.RelativeAzimuth3*pi/180.0)
        self.RelativeAzimuth4 = cos(self.g.RelativeAzimuth4*pi/180.0)
        self.RelativeAzimuth5 = cos(self.g.RelativeAzimuth5*pi/180.0)
        self.RelativeAzimuth6 = cos(self.g.RelativeAzimuth6*pi/180.0)
        self.RelativeAzimuth7 = cos(self.g.RelativeAzimuth7*pi/180.0)
        self.RelativeAzimuth8 = cos(self.g.RelativeAzimuth8*pi/180.0)
        self.RelativeAzimuth9 = cos(self.g.RelativeAzimuth9*pi/180.0)

        self.SensorZenith1   = cos(self.g.SensorZenith1*pi/180.0)
        self.SensorZenith2   = cos(self.g.SensorZenith2*pi/180.0)
        self.SensorZenith3   = cos(self.g.SensorZenith3*pi/180.0)
        self.SensorZenith4   = cos(self.g.SensorZenith4*pi/180.0)
        self.SensorZenith5   = cos(self.g.SensorZenith5*pi/180.0)
        self.SensorZenith6   = cos(self.g.SensorZenith6*pi/180.0)
        self.SensorZenith7   = cos(self.g.SensorZenith7*pi/180.0)
        self.SensorZenith8   = cos(self.g.SensorZenith8*pi/180.0)
        self.SensorZenith9   = cos(self.g.SensorZenith9*pi/180.0)

        # Sanity check
        # ------------
       
        #----------------------------------------------------------------------------
        "Commented out last check and replaced d.tau550 by m.tau558 ~ Suniyya Waraich"
        #----------------------------------------------------------------------------
        
        # NNR AOD, on demand
        # ------------------
        self.lnTau550 = None  # See tranSVC()
        self.surface = None
        self.laod = True 

        # Default input features
        # ----------------------
        if Input==None:
            self.Input = ( 'albedo', 'amf', 'lnTau550', 'ldTau550' )
        else:
            self.Input = Input

        # Default Target
        # --------------
        if Target==None:
            self.Target = ('laTau550',)
        else:
            self.Target = Target
 
#-------------------------------------------------------------------
def prepAncillary(year, month, npzDir='/nobackup/MAPSS/Collocation/MISR',
                  aer_x = '/nobackup/MERRAero/inst2d_hwl_x.ddf',
                  slv_Nx = '/nobackup/MERRA/slv_Nx',
                  AlbedoGF_Root = '/nobackup/10/MODIS/005/Level3/Albedo/data',
                  igbp_dir='/nobackup/Emissions/Vegetation/GL_IGBP_INPE'):
    """
    Prepare ancillary data, saving it to NPZ file.
    """

    # Get coordinates
    # ---------------
    npzFile = 'mapss.xtra_misr.%d4%02d.npz'%(year,month)
    if os.path.exists(npzFile):
        print(">< Skipping ancillaries on ", year, month)
        return
    else:
        print("<> Creatig ancillaries on ", year, month)
    a = MAPSS(npzDir+'/mapss.anet_misr.%4d%02d.npz'%(year,month))

    a.sample = None

    # Vegetation type
    # ---------------
    print("   o Sampling vegetation type")
    veg = a.detailedVeg(Path=igbp_dir)

    # Speciate
    # --------
    print("   o Speciating aerosols...")
    a.speciate(aer_x)
    
    # Wind speed
    # ----------
    print("   o Sampling 10M wind...")
    a.sampleFile(slv_Nx,onlyVars=('U10M','V10M'))
    a.wind = a.sample.U10M**2 + a.sample.V10M**2 

    # Ocean Albedo (still neds to be water masked later)
    # --------------------------------------------------
    print("<> Doing Cox Munk...")
    a.getCoxMunk()

    # Land Albedo
    # -----------
    print("   o Sampling surface albedo...")
    a.AlbedoSample(Verbose = True,
                   root=AlbedoGF_Root)

    # Save to NPZ file
    # ----------------
    savez(npzFile,
          version=1,nobs=len(a.lon), lon=a.lon,lat=a.lat,time=a.time,
          u10m = a.sample.U10M, v10m=a.sample.V10M, ocnAlbedo=a.ocnAlbedo,
          fdu=a.fdu, fss=a.fss, fbc=a.fbc, foc=a.foc, fcc=a.fcc, fsu=a.fsu,
          lndAlbedo=a.lndAlbedo, veg = a.veg)

    return a
    
#-------------------------------------------------------------------

__Months__ = ['JAN','FEB','MAR','APR','MAY','JUN','JUL','AUG','SEP','OCT','NOV','DEC']

def dt2gat(t):
    """
    Convert datetime to grads time.
    """
    gat = "%d:%dZ%d%s%d"%(t.hour,t.minute,t.day,__Months__[t.month-1],t.year)
    return gat

#....................................................................................     

def _deepNNR():

    Input_Raw = [ 'albedo', 'angstrom',
             'ScatteringAngle',
             'SensorAzimuth',
             'SensorZenith',
             'SolarAzimuth',
             'SolarZenith',
             'sRef412',
             'sRef470',
#             'sRef660',
             'dRef412',
             'dRef470',
#             'dRef660',
            ]

    # Did not quite work
    Input_Del = [ 'albedo',
             'ScatteringAngle',
             'SensorAzimuth',
             'SensorZenith',
             'SolarAzimuth',
             'SolarZenith',
             'xRef412',
             'xRef470',
            ]

    Input_Ret = ['albedo','angstrom',
             'ScatteringAngle',
             'SensorAzimuth',
             'SensorZenith',
             'SolarAzimuth',
             'SolarZenith',
             'dTau412',
             'dTau470',
             'dTau660',
            ]

    Target = ['aTau550',]

    Input = Input_Raw

    # Read and split dataset in training/testing subsets
    # --------------------------------------------------
    m = AQC_DEEP(verbose=True)

    #m.laod = False
    
    # Q/C
    # ---
    m.iValid =    (m.sRef412>0) & (m.sRef412<0.5) & \
                  (m.sRef470>0) & (m.sRef470<0.5) & \
                  (m.dRef412>0)    & \
                  (m.dRef470>0)    & \
                  (m.aTau550>0)    & \
                  (m.albedo>0.15)

#                  (m.dTau550>0)    & \
#                  (m.lon>-20)&(m.lon<60)&(m.lat>15)&(m.lat<35) & \
#                  ((m.albedo>0.15)|(m.d.qa_flag>1))
#                  (m.d.qa_flag>0)  & \
#                  (m.albedo>0.15)
#                  (m.albedo>0)     & \
    
    m.split()
    
    ident = 'mydd'
    expid = 'nnr_002'
    doInflation = False # In production, inflate NNR using parameters defined below
    doAQC = False # In production, use NNR for Q/C only

    # Training
    # --------
    nHidden = 2*len(Input)
#    topology = (len(Input), nHidden, 2, len(Target))
    topology = (len(Input), nHidden, len(Target))
    biases = True
        
    print(" ")
    print("        AOD Neural Net Retrieval")
    print("        ------------------------")
    print(" ")
    print("  No. Valid Data:  ", m.aTau550[m.iValid].size,  \
                                 int(100.*m.aTau550[m.iValid].size/m.aTau550.size),'%')
    print(" No. Hidden Nodes: ", nHidden)
    print("         Topology: ", topology)
    print("   Input Features: ", Input[:])
    print("           Target: ", Target[:])
    print(" ")

    m.train(Input=Input,Target=Target,nHidden=nHidden,
            topology=topology,biases=biases)
    out, reg = m.test()

    # Record inflation
    # ----------------
    if doInflation:
        how_much,tau0,dtau = (0.4,0.35,0.2) # hard-wired
        m.net.inflation = (how_much,tau0,dtau)

    # Signal NNR is to be used for Q/C
    # --------------------------------
    if doAQC:
        m.net.doAQC = True

    m.savenet(expid+"."+ident+'_Tau.net')

    return (m, expid, ident, Target)

#---
def inflate(tau,inflation):
    how_much,tau0,dtau = inflation
    a = (tau-tau0)/dtau
    f = (1+how_much/(1.+exp(-a)))
    return tau*f

def doPlots(m,expid,ident, Target):

    # All bright surfaces
    # -------------------
    m.iValid =    (m.sRef412>0) & (m.sRef412<0.5) & \
                  (m.sRef470>0) & (m.sRef470<0.5) & \
                  (m.dRef412>0)    & \
                  (m.dRef470>0)    & \
                  (m.aTau550>0)    & \
                  (m.dTau550>0)    & \
                  (m.d.qa_flag>0)  & \
                  (m.albedo>0.15)

    # Plot KDE of corrected AOD
    # -------------------------
    m.plotKDE(figfile=expid+"."+ident+"_kde-"+Target[0][1:]+"-corrected.png")

    # Plot KDE of uncorrected AOD
    # ---------------------------
    targets = m.getTargets(m.iValid).squeeze()
    if m.laod:
        formatter = aodFormat()
        bins = None
        original = m.ldTau550[m.iValid]
    else:
        formatter = None
        bins = arange(0., 0.6, 0.01 )
        original = m.dTau550[m.iValid]
    _plotKDE(targets,original,y_label='Original Deep Blue',x_bins=bins,
             formatter=formatter)
    if m.laod:
       title("Log("+Target[0][1:]+"+0.01)- "+ident)
    else:
       title(Target[0][1:]+"- "+ident)
    savefig(expid+"."+ident+"_kde-"+Target[0][1:]+'-orig.png')

    # Plot KDE of NNR x DB AOD
    # ------------------------
    results = m.eval(m.iValid).squeeze()
    if m.laod:
        original = m.ldTau550[m.iValid]
    else:
        original = m.dTau550[m.iValid]
    _plotKDE(results,original,y_label='Original Deep Blue',
             x_label='NNR',x_bins=bins,formatter=formatter)
    if m.laod:
       title("Log("+Target[0][1:]+"+0.01)- "+ident)
    else:
       title(Target[0][1:]+"- "+ident)
    savefig(expid+"."+ident+"_kde-"+Target[0][1:]+'-orig-new.png')

    # Plot KDE of QC'd Deep Blue
    # --------------------------
    I = abs(results-original)<0.75
    _plotKDE(targets[I],original[I],y_label='Original Deep Blue with Q/C',
             x_label='AERONET',x_bins=bins,formatter=formatter)
    if m.laod:
       title("Log("+Target[0][1:]+"+0.01) - "+ident)
    else:
       title(Target[0][1:]+" - "+ident)
    savefig(expid+"."+ident+"_kde-"+Target[0][1:]+'-orig-qc.png')

    # Plot NNR Error vs NNR AOD
    # -------------------------
    if m.laod:
        nTau = exp(results) - 0.01
    else:
        nTau = results
    diff = nTau - m.aTau550[m.iValid]
    _plotKDE(nTau,diff,y_label='NNR Error',x_label='NNR',
             x_bins=arange(0., 0.30, 0.01 ),
             y_bins=arange(-0.15, 0.15, 0.01 ))
    title(Target[0][1:]+" - "+ident)
    savefig(expid+"."+ident+"_kde-"+Target[0][1:]+'-diff.png')

    # Plot Inflated Tau
    # -----------------
    try:
        lnTau = log(inflate(nTau,m.net.inflation)+0.01)
        _plotKDE(m.laTau550[m.iValid],lnTau,y_label='Inflated NNR',x_label='AERONET',
                 formatter=formatter)
        title(Target[0][1:]+" - "+ident)
        savefig(expid+"."+ident+"_kde-"+Target[0][1:]+'-infl.png')
    except:
        pass

    # North Africa
    # ------------
    m.iValid = m.iValid & \
               (m.lon>-20)&(m.lon<60)&(m.lat>15)&(m.lat<35)
    m.plotKDE(figfile=expid+"."+ident+"_kde-"+Target[0][1:]+"-nnr-sahara.png")

    # Inflated over North Africa
    # --------------------------
    try:
        results = m.eval(m.iValid).squeeze()
        if m.laod:
            nTau = exp(results) - 0.01
        else:
            nTau = results
        lnTau = log(inflate(nTau,m.net.inflation)+0.01)
        _plotKDE(m.laTau550[m.iValid],lnTau,y_label='Inflated NNR North Africa',
                 x_label='AERONET',
                 formatter=formatter)
        title(Target[0][1:]+" - "+ident)
        savefig(expid+"."+ident+"_kde-"+Target[0][1:]+'-infl-sahara.png')
    except:
        pass

    # Plot KDE of uncorrected AOD
    # ---------------------------
    targets = m.getTargets(m.iValid).squeeze()
    if m.laod:
        original = m.ldTau550[m.iValid]
    else:
        original = m.dTau550[m.iValid]
    _plotKDE(targets,original,y_label='Original Deep Blue',x_bins=bins,
             formatter=formatter)
    if m.laod:
        title("Log("+Target[0][1:]+"+0.01)- "+ident)
    else:
        title(Target[0][1:]+"- "+ident)
    savefig(expid+"."+ident+"_kde-"+Target[0][1:]+'orig-sahara.png')

    # Scatter diagram for testing
    # ---------------------------
    # m.plotScat(figfile=expid+"."+ident+"_scat-"+Target[0][1:]+'.png')

#....................................................................................     

if __name__ == "__main__":

    for year in range(2003,2015):
        for month in range(1,13):
            a = prepAncillary(year,month)

def hold():

    m, expid, ident, Target = _deepNNR()
    

    m = AQC_DEEP(verbose=True)
    m.getNNR()
    m.trainSVC()

    I = m.svcIndex
    J = I & (m.svcEval==1)&(m.albedo>0.15)&(m.ldTau550>-5)

    m, expid, ident, Target = _deepNNR()
    doPlots(m, expid, ident, Target)


    m = AQC_DEEP(verbose=True)
    m.getNNR()
    m.trainSVC()
    I = m.svcIndex
    J = I & (m.svcEval==1)&(m.albedo>0.15)&(m.ldTau550>-5) 
