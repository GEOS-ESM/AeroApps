"""
This module implements the unbiasing of MODIS AOD retrievals based
on the (Neural Net) AOD Bias Correction model developed in module
abc_modis.

  This version works from MODIS MOD04/MYD04 Level 2 files.

"""
import os, sys
sys.path.insert(0,'/home/pcastell/Enthought/Canopy_64bit/System/lib/python2.7/site-packages')
import warnings
from pyobs.mxd04 import MxD04_L2, MISSING, granules

from ffnet import loadnet
from numpy import  c_ as cat
from numpy import  copy, ones, sin, cos, exp, arccos, pi, any, log

# SDS to be read in
# ------------
SDS = dict( META =    ( "Scan_Start_Time",
                        "Latitude",
                        "Longitude",
                        "Solar_Zenith",
                        "Solar_Azimuth",
                        "Sensor_Zenith",
                        "Sensor_Azimuth",
                        "Scattering_Angle",
                        "Glint_Angle"),

            LAND =    ( 'Mean_Reflectance_Land',
                        'Surface_Reflectance_Land',
                        'Aerosol_Cloud_Fraction_Land',
                        'Quality_Assurance_Land'),

            OCEAN =   ( 'Mean_Reflectance_Ocean',
                        'Aerosol_Cloud_Fraction_Ocean',               
                        'Quality_Assurance_Ocean'),

            DEEP =    ( 'Deep_Blue_Spectral_TOA_Reflectance_Land',
                        'Deep_Blue_Spectral_Surface_Reflectance_Land',
                        'Deep_Blue_Cloud_Fraction_Land',
                        'Deep_Blue_Aerosol_Optical_Depth_550_Land_QA_Flag'))

# Translate Inputs between NNR and MODIS classes
# -----------------------------------------------
TranslateInput = dict ( OCEAN = dict( mRef470  = ('reflectance',470),
                                      mRef550  = ('reflectance',550),
                                      mRef660  = ('reflectance',660),
                                      mRef870  = ('reflectance',870),
                                      mRef1200 = ('reflectance',1200),
                                      mRef1600 = ('reflectance',1600),
                                      mRef2100 = ('reflectance',2100)),
                                      
                        LAND =  dict( mSre470  = ('sfc_reflectance',470),
                                      mSre660  = ('sfc_reflectance',660),
                                      mSre2100 = ('sfc_reflectance',2100),
                                      mRef412  = ('reflectance',412),
                                      mRef440  = ('reflectance',440),
                                      mRef470  = ('reflectance',470),
                                      mRef550  = ('reflectance',550),
                                      mRef660  = ('reflectance',660),
                                      mRef870  = ('reflectance',870),
                                      mRef1200 = ('reflectance',1200),
                                      mRef1600 = ('reflectance',1600),
                                      mRef2100 = ('reflectance',2100)),

                        DEEP =  dict( mSre412  = ('sfc_reflectance',412),
                                      mSre470  = ('sfc_reflectance',470),
                                      mSre660  = ('sfc_reflectance',660),
                                      mRef412  = ('reflectance',470),
                                      mRef470  = ('reflectance',550),
                                      mRef660  = ('reflectance',660)))

for var in ( 'ScatteringAngle','GlintAngle',
             'SolarAzimuth', 'SolarZenith',
             'SensorAzimuth','SensorZenith',
             'cloud','qa_flag'  ):
    TranslateInput[var] = (var,)

# Translate Targets between ANET and MODIS classes
# ------------------------------------------------
TranslateTarget = dict ( aTau470 = ( 'aod_', 470 ),
                         aTau550 = ( 'aod_', 550 ),
                         aTau660 = ( 'aod_', 660 ),
                         aTau870 = ( 'aod_', 870 ),
                         )

class MxD04_NNR(MxD04_L2):
    """
    This class extends MODIS by adding methods for producing
    NNR AOD retrievals based on the Neural Net model developed
    with class *abc_modis*.
    """

    def __init__(self,l2_path,prod,algo,syn_time,
                 ga,expr='mag(u10m,v10m)',vname='wind',
                 coxmunk_lut='coxmunk_lut.npz',
                 cloud_thresh=0.70,alb_min=0.25,coll='051',verbose=0):
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
         alb_min --- dark target will consider albedos < alb_min,
                     while deep blue will consider albedos > alb_bin
              
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
        if algo == "LAND" or algo == "DEEP":
            # self.aod_fine = self.aod_fine[:,0:3] # Drop 2100 nm
            self.aod_fine = MISSING * ones(self.aod.shape) # no longer provided in col. 6
        
        i = (self.aod<MISSING)&(self.aod_fine<MISSING)
        self.aod_coarse = MISSING * ones(self.aod.shape)
        self.aod_coarse[i] = self.aod[i] - self.aod_fine[i]

        # Add wind or albedo variable
        # ---------------------------
        if vname == 'albedo':
            self.addVar(ga,expr=expr,vname='albedo',clmYear=2000)
        else:
            self.addVar(ga,expr=expr,vname=vname)
            self.getCoxMunk(coxmunk_lut) # ocean albedo

        # Q/C
        # ---
        self.iGood = (self.qa_flag>0)             & \
                     (self.cloud<cloud_thresh) 

        if vname == 'albedo':
#            print "ALBEDO: ", len(self.albedo), len(self.lat)
#            print "ALBEDO: ", self.albedo.min(), self.albedo.max()
            if algo == 'LAND': # Dark Target (albedo<0.15)
                self.iGood = self.iGood           & \
                    (self.qa_flag==3)             & \
                    (self.ScatteringAngle<170.)   & \
                    (self.albedo>0)               & \
                    (self.albedo<=alb_min)
            else: # Deep Blue (Bright Target, albedo>0.15)
                sRef412 = self.Deep_Blue_Surface_Reflectance_Land[:,0]
                sRef470 = self.Deep_Blue_Surface_Reflectance_Land[:,1]
                dRef412 = self.Deep_Blue_Mean_Reflectance_Land[:,0]
                dRef470 = self.Deep_Blue_Mean_Reflectance_Land[:,1]
                self.orig = self.Deep_Blue_Aerosol_Optical_Depth_Land[:,2] # AOD 550
                self.iGood = self.iGood           & \
                    (self.ScatteringAngle<170.)   & \
                    (sRef412>0) & (sRef412<0.5) & \
                    (sRef470>0) & (sRef470<0.5) & \
                    (dRef412>0) & (dRef412<0.5) & \
                    (dRef470>0) & (dRef470<0.5) & \
                    (self.albedo>alb_min)
                    
        if vname == 'wind':
            self.iGood = self.iGood               & \
                         (self.wind>0)             

        if any(self.iGood) == False:
            print "WARNING: Strange, no good obs left to work with; make sure wind speed/albedo is available."
            return

        # Make sure same good AOD is kept for gridding
        # --------------------------------------------
        self.aod[self.iGood==False,:] = MISSING

        # Redefine Angstrom exponent for Deep Blue
        # ----------------------------------------
        if algo == "DEEP":
            I, J = (self.iGood==True, self.iGood==False)
            i470 = self.dChannels.index(470)
            i660 = self.dChannels.index(660)
            self.angstrom[I] = - log(self.aod[I,i660]/self.aod[I,i470])/log(660./470.)
            self.angstrom[J] = MISSING
            print "- Redefined angstrom exponent for 660/470: min=%4.2f, max=%4.2f"%\
                (self.angstrom[I].min(), self.angstrom[I].max())

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
                if 'mSre' in inputName: # LAND, surface reflectivity
                    k = list(self.sChannels).index(ch) # index of channel (all 7)
                elif 'mRef' in inputName: # LAND/OCEAN reflectances
                    k = list(self.Channels).index(ch) # index of channel (all 7)
                elif ('dRef' in inputName) or 'sRef' in inputName: # deep blue toa/sfc reflectance 
                    k = list(self.dChannels).index(ch) # index of channel (all 7)
                else:
                    k = list(self.channels).index(ch)  # algorithm specific channels
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
            if self.verbose>0:
                print "Ubiasing ", name, ch, 'Log-AOD = ',self.net.laod 
            k = list(self.channels).index(ch) # index of channel            
            self.channels_ = self.channels_ + [ch,]
            if self.net.laod:
                result = exp(targets[:,i]) - 0.01 # inverse
            else:
                result = targets[:,i]
            result = self.inflate(result) # inflate results, if desired
            self.__dict__[name][self.iGood,k] = result

            i += 1 

        return result

#---
    def inflate(self,result):
        """
        If desired, inflate the AOD estimate.
        """
        if 'inflation' in self.net.__dict__.keys():
            how_much, tau0, dtau = self.net.inflation
            result = inflate(result,how_much=how_much,tau0=tau0,dtau=dtau)
            print "- Inflating AOD by %4.2f with tau0=%4.2f and dtau=%4.2f"%(how_much,tau0,dtau)
        else:
            print "- No inflation being applied, skipping..."

        return result

#---
        
    __call__= apply

#---
def inflate(tau,how_much=0.5,tau0=0.35,dtau=0.1):
    a = (tau-tau0)/dtau
    f = (1+how_much/(1.+exp(-a)))
    return tau*f

#---

if __name__ == "__main__":

    from datetime import datetime
    from grads import GrADS

    l2_path = '/nobackup/MODIS/Level2/'
    algo = 'DEEP'
    prod = 'MYD04'
    coll = '051'
    syn_time = datetime(2008,6,30,12,0,0)
    albedo_file = '/nobackup/MODIS/Level3/ALBEDO/albedo_clim.ctl'

    ga = GrADS(Echo=False,Window=False)
    ga('open %s'%albedo_file)
    expr='albedo'
    vname = 'albedo'

    m = MxD04_NNR(l2_path,prod,algo.upper(),syn_time,
                  ga,expr=expr,vname=vname,coll=coll,
                  cloud_thresh=0.7,
                  verbose=True)

    laod = m.apply('/nobackup/NNR/Net/nnr_002.mydd_Tau.net')
    aod = exp(laod) - 0.01
