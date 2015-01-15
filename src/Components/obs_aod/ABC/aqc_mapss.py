"""

Playing with Deep Blue.

"""

import os

from numpy import pi, cos, log, zeros, ones, savez, arange, exp
from numpy import c_ as cat

#from ffnet        import loadnet
from sknet        import loadnet
from pyobs        import NPZ

from matplotlib.pyplot import title, savefig
from pyobs.mapss       import MAPSS
from abc_modis         import NN, _plotKDE, aodFormat

MISSING = -1.0e20
d2r = pi / 180.

class AQC_DEEP(NN):

    def __init__(self,
                 npzDir='./NPZ',
                 albedoNPZ='mapss-deep_albedo.npz',
                 tol=0.5,
                 Input=None,Target=None,
                 verbose=False):

        self.verbose = verbose
        self.ident = 'mydd' # deeb blue land
        
        # Read in NPZ files written by collocation app
        # --------------------------------------------
        self.a = MAPSS(npzDir+'/mapss.anet.??????.npz')
        self.d = MAPSS(npzDir+'/mapss.dtau.??????.npz')
        self.r = MAPSS(npzDir+'/mapss.dref.??????.npz')
        self.s = MAPSS(npzDir+'/mapss.sref.??????.npz')

        # Inherit coordinates from AERONET file
        # -------------------------------------
        self.lon, self.lat, self.time, self.N = (self.a.lon, self.a.lat, self.a.time, self.a.nobs)

        # Albedo
        # ------
        if not os.path.exists(npzDir+'/'+albedoNPZ):
            self.getAlbedo(npzDir+'/'+albedoNPZ)
        else:
            self.albedo = NPZ(npzDir+'/'+albedoNPZ).var

        # Air mass factor
        # ---------------
        self.amf = (1./cos(d2r*self.r.SolarZenith))+(1./cos(d2r*self.r.SensorZenith))  

        # Expose reflectances
        # -------------------
        self.sRef412 = self.s.sRef412
        self.sRef470 = self.s.sRef470
        self.sRef660 = self.s.sRef660
        self.dRef412 = self.r.dRef412
        self.dRef470 = self.r.dRef470
        self.dRef660 = self.r.dRef660
        self.xRef412 = self.dRef412 - self.sRef412 
        self.xRef470 = self.dRef470 - self.sRef470 
 
        # Expose AOD
        # ----------
        self.aTau440  = self.a.tau440
        self.aTau550  = self.a.tau550
        self.dTau412  = self.d.tau550
        self.dTau470  = self.d.tau470
        self.dTau550  = self.d.tau550
        self.dTau660  = self.d.tau660
        self.angstrom = -log(self.dTau660/self.dTau470)/log(660./470.)
        self.laTau550 = log(self.a.tau550+0.01)
        self.ldTau550 = log(self.d.tau550+0.01)

        # Angle transforms: for NN calculations we work with cosine of angles
        # -------------------------------------------------------------------
        self.ScatteringAngle = cos(self.r.ScatteringAngle*pi/180.0) 
        self.SensorAzimuth   = cos(self.r.SensorAzimuth*pi/180.0)   
        self.SensorZenith    = cos(self.r.SensorZenith*pi/180.0)    
        self.SolarAzimuth    = cos(self.r.SolarAzimuth*pi/180.0)    
        self.SolarZenith     = cos(self.r.SolarZenith*pi/180.0)     

        # Sanity check
        # ------------
        self.iValid = (self.a.tau550>-0.01) &\
                      (self.albedo>0)        &\
                      (self.d.tau550>-0.01) &\
                      (self.d.qa_flag>0)
                     
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

        
    def addVar(self,ga,outfile,expr='ustar',vname=None, clmYear=None):
        """
        Given a grads object having the correct file as default,
        writes out a CSV file with the 3 variables. When *clmYear* is
        specified, the actual year in the time attribute is replaced
        with he climatological year *clmYear*.

        This algorithm uses a *nearest neighbor* interpolation.
        
        """

        N = self.N
        U = ones(N)
        U[:] = MISSING

        if vname == None:
            vname = expr

        for i in range(N):

            x = self.lon[i]
            y = self.lat[i]

            if clmYear == None:
                t = self.time[i]
            else:
                t = dt2gat(self.time[i])
                t = t[:-4] + str(clmYear) # replace year

            ga('set lon %f %f'%(x-1.,x+1.),Quiet=True)
            ga('set lat %f %f'%(y-1.,y+1.),Quiet=True)
            ga('set time %s'%t,Quiet=True)

            try:
                u, levs = ga.interp(expr, lons=(x,),lats=(y,))
                U[i] = u.data
            except:
                ga.flush()

            if U[i] >= 0.0:
                print self.time[i], "%8.3f %8.3f %6.3f  ...%8.3f%%"\
                      %(x,y,U[i],100.*i/float(N))

        self.__dict__[vname] = U

        version = 1
        meta = [ version, vname, expr ]
        savez(outfile,meta=meta,lon=self.lon,lat=self.lat,time=self.time,var=U)

    def getAlbedo(self,npzFile):
        from grads import GrADS
        ga = GrADS(Echo=False,Window=False)
        ga('open albedo_clim.ctl')
        self.addVar(ga,npzFile,expr='albedo',clmYear=2000)

    def trainNovel(self,nu=0.1,kernel='rbf',gamma=0.1,**kwopts):
        """
        Train Novelty SVM and apply it to "bad" data.                                           
        NoteL I did not work too well with refletances...
        """
        from sklearn import svm

        self.clf = svm.OneClassSVM(nu=nu, kernel=kernel, gamma=gamma,**kwopts)

        # Fit using good data
        # -------------------
        Xgood = self.getInputs(I=self.iGood)
        self.clf.fit(Xgood)
        self.NovelTrain = self.clf.predict(Xgood)
        self.NovelTrate = 100.*self.NovelTrain[self.NovelTrain==1].size/self.NovelTrain.size

        # Eval on bad data
        # ----------------
        Xbad = self.getInputs(I=self.iBad)
        self.NovelBad = self.clf.predict(Xbad)
        self.NovelBrate = 100.*self.NovelBad[self.NovelBad==-1].size/self.NovelBad.size

    def trainSVC(self,Verbose=True,tol=0.5,alb_min=0.15,doScores=False,**kwopts):
        """
        Train Novelty SVM and apply it to "bad" data.                                           
        """
        from sklearn import svm, cross_validation

        # Segregate retrievals into good/bad according to AERONET
        # -------------------------------------------------------
        diff = self.ldTau550-self.laTau550
        self.iGood = self.iValid&(abs(diff)<=tol)
        self.iBad  = self.iValid&(abs(diff)>tol)   # too large or too small
        self.iBad1 = self.iValid&(diff<-tol)       # too large
        self.iBad2 = self.iValid&(diff>+tol)       # too small

        # Instantiate classifier
        # ----------------------
        self.clf = svm.SVC(**kwopts)

        Input = self.Input

        # Features and class labels
        # -------------------------
        I = (self.iValid)&(self.d.qa_flag>0)&(self.lnTau550>-5)&(m.albedo>alb_min)
        X = self.getInputs(I=I,Input=Input)
        y = ones(self.N)  # Good = 1
        y[self.iBad] = -1 # Bad = -1
        y = y[I]

        # Train the classifier
        # --------------------
        if Verbose:
            print "- Training SVM classifier with %d samples"%y.size
        self.svc = svm.SVC(**kwopts)
        self.svc.fit(X,y)

        # Save for diagnostics
        # --------------------
        self.svcTarget   = MISSING * ones(self.N)
        self.svcEval     = MISSING * ones(self.N)

        self.svcIndex       = I       # Indices used for training
        self.svcTarget[I]   = y       # Correct outcome
        self.svcEval[I]     = self.svc.predict(X) # Trained outcome
        # self.svcFeature     = Input   # Inputs used

        # Calculate scores
        # ----------------
        if Verbose:
            print "- Calculating Cross Validated scores..."
        self.svcScores = cross_validation.cross_val_score(self.svc, X, y, n_jobs=1)
        if Verbose:
            print "- Cross Validated scores are ", self.svcScores

    def getNNR(self):
        """
        Apply pre-computed NNR to reflectances.
        """

        # Load network
        # ------------
        self.net = loadnet('nnr_001.mydd_Tau.net')

        Input = self.net.InputNames
        
        # Q/C
        # ---
        iValid =  (m.sRef412>0)    & \
                  (m.sRef470>0)    & \
                  (m.dRef412>0)    & \
                  (m.dRef470>0)    & \
                  (m.aTau550>0)    & \
                  (m.dTau550>0)    & \
                  (m.d.qa_flag>0)  & \
                  (m.albedo>0)

        # Evaluate
        # --------
        target = MISSING * ones(self.lon.shape)
        target[iValid] = self.net(self.getInputs(I=iValid,Input=Input))

        self.lnTau550 = target
 
#---

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
        
    print " "
    print "        AOD Neural Net Retrieval"
    print "        ------------------------"
    print " "
    print "  No. Valid Data:  ", m.aTau550[m.iValid].size,  \
                                 int(100.*m.aTau550[m.iValid].size/m.aTau550.size),'%'
    print " No. Hidden Nodes: ", nHidden
    print "         Topology: ", topology
    print "   Input Features: ", Input[:]
    print "           Target: ", Target[:]
    print " "

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

    m, expid, ident, Target = _deepNNR()
    
def hold():

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
