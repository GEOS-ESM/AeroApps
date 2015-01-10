"""
This module implements the unbiasing of MODIS AOD retrievals based
on the (Neural Net) AOD Bias Correction model developed in module
abc_modis.

  This version works from MODIS MOD04/MYD04 Level 2 files.

"""

import warnings
from pyobs.mxd04 import MxD04_L2, MISSING, granules

from ffnet import loadnet
from numpy import  c_ as cat
from numpy import  copy, ones, sin, cos, exp, arccos, pi, any

# Translate Inputs between ANET and MODIS classes
# -----------------------------------------------
TranslateInput = dict ( mTau470 = ( 'aod', 470 ),
                        mTau550 = ( 'aod', 550 ),
                        mTau660 = ( 'aod', 660 ),
                        mTau870 = ( 'aod', 870 ),
                        mTAU470 = ( 'aod_coarse', 470),
                        mTAU550 = ( 'aod_coarse', 550),
                        mTAU660 = ( 'aod_coarse', 660),
                        mTAU870 = ( 'aod_coarse', 870),
                        mRef470  = ('reflectance',470),
                        mRef550  = ('reflectance',550),
                        mRef660  = ('reflectance',660),
                        mRef870  = ('reflectance',870),
                        mRef1200 = ('reflectance',1200),
                        mRef1600 = ('reflectance',1600),
                        mRef2100 = ('reflectance',2100),
                        )
for var in ( 'ScatteringAngle','GlintAngle',
             'SolarAzimuth', 'SolarZenith',
             'SensorAzimuth','SensorZenith',
             'cloud', 'wind', 'ustar', 'albedo' ):
    TranslateInput[var] = (var,)

# Translate Targets between ANET and MODIS classes
# ------------------------------------------------
TranslateTarget = dict ( aTau470 = ( 'aod_', 470 ),
                         aTau550 = ( 'aod_', 550 ),
                         aTau660 = ( 'aod_', 660 ),
                         aTau870 = ( 'aod_', 870 ),
                         )

class Unbias(MxD04_L2):
    """
    This class extends MODIS by adding methods for unbiasing
    the AOD retrievals based on the Neural Net model developed
    with class *abc_modis*.
    """

    def __init__(self,l2_path,prod,algo,syn_time,
                 ga,expr='mag(u10m,v10m)',vname='wind',
                 cloud_thresh=0.85,coll='051',verbose=0):
        """
        Contructs a MXD04 object from MODIS Aerosol Level 2
        granules. On input,

         l2_path --- top directory for the MODIS Level 2 files;
                      it must have subdirectories MOD04 and MYD04.
            prod --- either *MOD04* (Terra) or *MYD04* (Aqua)
            algo --- aerosol algorithm: LAND, OCEAN or DEEP (for
                      Deep Blue)
        syn_time --- synoptic time

              ga --- GrADS object for wind parameter
            expr --- *wind* variable expression
           vname --- name of wind variable
     cloud_tresh --- cloud fraction treshhold
              
        The following attributes are also defined:
           GlintAngle
           aod_coarse
           wind
           
        It also performs Q/C, setting attribute iGood. On,
        input, *cloud_thresh* is the cloud fraction limit.
        """

        self.verbose = verbose
        self.vname = vname
        
        # Initialize superclass
        # ---------------------
        Files = granules(l2_path,prod,syn_time,coll=coll)
        MxD04_L2.__init__(self,Files,algo,syn_time,only_good=True,Verb=verbose)

        if self.nobs < 1:
            return # no obs, nothing to do

        # Create attribute for holding the revised AOD
        # --------------------------------------------
        self.aod_ = MISSING * ones(self.aod.shape)
        
        # Compute coarse mode AOD
        # -----------------------
        if algo == "LAND":
            self.aod_fine = self.aod_fine[:,0:3] # Drop 2100 nm
        i = (self.aod<MISSING)&(self.aod_fine<MISSING)
        self.aod_coarse = MISSING * ones(self.aod.shape)
        self.aod_coarse[i] = self.aod[i] - self.aod_fine[i]

        # Add wind or albedo variable
        # ---------------------------
        if vname == 'albedo':
            self.addVar(ga,expr=expr,vname='albedo',clmYear=2000)
        else:
            self.addVar(ga,expr=expr,vname=vname)

        # Q/C
        # ---
        self.iGood = (self.qa_flag>0)          & \
                     (self.cloud<cloud_thresh) 

        if vname == 'albedo':
            print "ALBEDO: ", len(self.albedo), len(self.lat)
            print "ALBEDO: ", self.albedo.min(), self.albedo.max()
            self.iGood = self.iGood                   & \
                         (self.ScatteringAngle<170.)  & \
                         (self.albedo>0)              & \
                         (self.albedo<0.25)

        if vname == 'wind':
            self.iGood = self.iGood                   & \
                         (self.wind>0)             

        if any(self.iGood) == False:
            print "WARNING: Strange, no good obs left to work with; make sure wind speed/albedo is available."
            return

        # Make sure same good AOD is kept for gridding
        # --------------------------------------------
        self.aod[self.iGood==False,:] = MISSING

        # Glint Angle
        # -----------
        d2r = pi / 180.
        r2d = 180. / pi
        RelativeAzimuth = abs(self.SolarAzimuth - self.SensorAzimuth - 180.)
        cosGlintAngle = cos(self.SolarZenith*d2r) * cos(self.SensorZenith*d2r) + \
                        sin(self.SolarZenith*d2r) * sin(self.SensorZenith*d2r) * \
                        cos(RelativeAzimuth*d2r)
        
#        i = (abs(cosGlintAngle)<=1.0)
#        self.GlintAngle = MISSING * ones(self.SolarZenith.shape)
#        self.GlintAngle[i] = arccos(cosGlintAngle[i])*r2d

        # Angle transforms: for NN calculations we work with cosine of angles
        # -------------------------------------------------------------------
        self.ScatteringAngle = cos(self.ScatteringAngle*pi/180.0) 
        self.SensorAzimuth   = cos(self.SensorAzimuth*pi/180.0)   
        self.SensorZenith    = cos(self.SensorZenith*pi/180.0)    
        self.SolarAzimuth    = cos(self.SolarAzimuth*pi/180.0)    
        self.SolarZenith     = cos(self.SolarZenith*pi/180.0)     
        self.GlintAngle      = cosGlintAngle

    def _loadNet(self,nnFile):
        """
        Loads the Neural Net weights created with class ABC.
        """
        self.net = loadnet(nnFile)

    def _getInputs(self):
        """
        Get Inputs for Neural Net.
        """

        # Loop over inputs
        # ----------------
        first = True
        for inputName in self.net.InputNames:
            iName = TranslateInput[inputName]

            if self.verbose>0:
                print 'Getting NN input ',iName

            # Retrieve input
            # --------------
            if len(iName) == 2:
                name, ch = iName
                k = list(self.Channels).index(ch) # index of channel
                input = self.__dict__[name][:,k]
                self.iGood = self.iGood & (input!=MISSING) # Q/C
            elif len(iName) == 1:
                name = iName[0]
                input = self.__dict__[name][:]
                self.iGood = self.iGood & (input!=MISSING) # Q/C
            else:
                raise ValueError, "strange, len(iName)=%d"%len(iName)

            # Concatenate Inputs
            # ------------------
            if first:
                inputs = input
                first = False
            else:
                inputs = cat[inputs,input]

        # Keep only good observations
        # ---------------------------
        return inputs[self.iGood,:]

#--
    def apply(self,nnFile):
        """
        Apply bias correction to AOD.
        """

        # Stop here is no good obs available
        # ----------------------------------
        if self.nobs == 0:
            return # no data to work with
        if any(self.iGood) == False:
            return # no good data to work with

        # Load the Neural Net
        # -------------------
        self._loadNet(nnFile)

        # Evaluate NN on inputs
        # ---------------------
        targets = self.net(self._getInputs())

        # Replace targets with unbiased values
        # ------------------------------------
        i = 0
        self.channels_ = [] # channels being revised
        for targetName in self.net.TargetNames:
            name, ch = TranslateTarget[targetName]
            k = list(self.Channels).index(ch) # index of channel            
            self.channels_ = self.channels_ + [ch,]
            if self.verbose>0:
                print "Ubiasing ", name, ch, 'Log-AOD = ',self.net.laod 
            if self.net.laod:
                self.__dict__[name][self.iGood,k] = exp(targets[:,i]) - 0.01 # inverse
            else:
                self.__dict__[name][self.iGood,k] = targets[:,i]
            i += 1 

        return targets

    __call__= apply

#---

if __name__ == "__main__":

    pass
