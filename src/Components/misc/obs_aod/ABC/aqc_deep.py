"""

Playing with Deep Blue.

     *** OLD CSV VERSION ****

"""


from numpy import pi, cos, log, zeros, ones
from numpy import c_ as cat

from ffnet import loadnet
from pyobs import NPZ

from anet        import *
from abc_modis   import NN, _plotKDE

d2r = pi / 180.

class AQC_DEEP(DEEP):

    def __init__(self,fname,csvVersion=1,tol=0.75,cloud_thresh=70,
                 albedoNPZ='NPZ/deep_albedo.npz',
                 Input=None,Target=None):
            
        DEEP.__init__(self,fname, csvVersion=csvVersion)


        # Albedo
        # ------
        if albedoNPZ==None:
            self.getAlbedo()
        else:
            self.albedo = NPZ(albedoNPZ).var
    
        self.laTau550 = log(self.aTau550+0.01)
        self.lmTau550 = log(self.mTau550+0.01)
        self.amf = (1./cos(d2r*self.SolarZenith))+(1./cos(d2r*self.SensorZenith))  

        # Sanity check
        # ------------
        self.iValid =(self.aTau550>0)&\
                     (self.mTau550>0)&\
                     (self.albedo>0)&\
                     (self.cloud<cloud_thresh)
                     
        # Segregate retrievals into good/bad according to AERONET
        # -------------------------------------------------------
        diff = self.lmTau550-self.laTau550
        self.iGood = self.iValid&(abs(diff)<=tol)
        self.iBad  = self.iValid&(abs(diff)<tol) # too large or too small
        self.iBad1 = self.iValid&(diff<-tol)     # too large
        self.iBad2 = self.iValid&(diff>+tol)     # too small

        # Input features
        # --------------
        if Input==None:
            self.Input = (   'albedo',
                             'amf',
                             'mean_mref0470_l',
                             'mean_mref0550_l',
                             'mean_mref0660_l',
                             'mean_mref0870_l',
                             'mean_mref1200_l',
                             'mean_mref1600_l',
                             'mean_mref2100_l',
                             'mean_surfre0470dpbl_l',
                             'mean_surfre0660dpbl_l',
                        )
        else:
            self.Input = Input

        # Target
        # ------
        if Target==None:
            self.Target = ('laTau550',)
        else:
            self.Target = Target

        # Angle transforms: for NN calculations we work with cosine of angles
        # -------------------------------------------------------------------
        self.ScatteringAngle = cos(self.ScatteringAngle*pi/180.0) 
        self.SensorAzimuth   = cos(self.SensorAzimuth*pi/180.0)   
        self.SensorZenith    = cos(self.SensorZenith*pi/180.0)    
        self.SolarAzimuth    = cos(self.SolarAzimuth*pi/180.0)    
        self.SolarZenith     = cos(self.SolarZenith*pi/180.0)     

        # NNR AOD, on demand
        # ------------------
        self.lnTau550 = None

    def getInputs(self,I,Input=None):
        """
        Given a set of indices *I*, returns the corresponding
        inputs for a neural net evaluation.
        Returns: inputs
        """
        if Input==None:
            Input = self.Input
        inputs = self.__dict__[Input[0]][I]
        for var in Input[1:]:
            inputs = cat[inputs,self.__dict__[var][I]]
        return inputs
    
    def getTargets(self,I):
        """
        Given a set of indices *I*, return the corresponding
        targets for a neural net evaluation:
        Returns: tagets
        """
        targets = self.__dict__[self.Target[0]][I]
        for var in self.Target[1:]:
            targets = cat[targets,self.__dict__[var][I]]
        return targets
 
    def getAlbedo(self):
        from grads import GrADS
        ga = GrADS(Echo=False,Window=False)
        ga('open albedo_clim.ctl')
        self.addVar(ga,'deep_albedo.npz',expr='albedo',clmYear=2000)

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

    def trainSVC(self,Verbose=True,**kwopts):
        """
        Train Novelty SVM and apply it to "bad" data.                                           
        """
        from sklearn import svm, cross_validation

        self.clf = svm.SVC(**kwopts)

        Input = (   'albedo',
                    'amf',
                    'lmTau550',
                    'lnTau550'
                )

        # Make sure NNR AOD is available
        # ------------------------------
        if self.lnTau550==None:
            if Verbose:
                print("- Evaluating NNR")
            self.getNNR()
        
        # Features and class labels
        # -------------------------
        I = (self.lnTau550>=log(0.01)) # Valid NNR retrievals
        X = self.getInputs(I=I,Input=Input)
        y = zeros(self.N)  # Good
        y[self.iBad1] = -1 # too small
        y[self.iBad2] = +1 # too large
        y = y[I]

        # Train the classifier
        # --------------------
        if Verbose:
            print("- Training SVM classifier with %d samples"%y.size)
        self.svc = svm.SVC(**kwopts)
        self.svc.fit(X,y)

        # Save for diagnostics
        # --------------------
        self.svcFeatures = Input
        self.svcTarget = y
        self.svcIndex = I
        self.svcEval = self.svc.predict(X)        

        # Calculate scores
        # ----------------
        if Verbose:
            print("- Calculating Cross Validated scores...")
        self.svc.Score = cross_validation.cross_val_score(self.svc, X, y, n_jobs=1)
        
    def getNNR(self):
        """
        Apply pre-computed NNR to reflectances.
        """

        Input_nnr = ['mRef550','mRef470','mRef660', 'mRef870',
                     'mRef1200','mRef1600','mRef2100',
                     'ScatteringAngle',
                     'SolarAzimuth', 'SolarZenith',
                     'SensorAzimuth','SensorZenith',
                     'cloud', 'albedo' ]

        # Load network
        # ------------
        self.net = loadnet('Net/nnr_001.mydl_Tau.net')

        # Evaluate
        # --------
        iGood =       self.iValid               & \
                      (self.mRef470 > 0.0)      & \
                      (self.mRef550 > 0.0)      & \
                      (self.mRef660 > 0.0)      & \
                      (self.mRef870 > 0.0)      & \
                      (self.mRef1200> 0.0)      & \
                      (self.mRef1600> 0.0)      & \
                      (self.mRef2100> 0.0)    

        target = MISSING * ones(self.lon.shape)
        target[iGood] = self.net(self.getInputs(I=iGood,Input=Input_nnr))

        self.lnTau550 = target
        
if __name__ == "__main__":

    d = AQC_DEEP('SUPER2_combo.Aqua.csv',csvVersion=2)
