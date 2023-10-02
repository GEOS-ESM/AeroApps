"""
This module implements the MODIS NNR AOD retrievals.

This version works from MODIS MOD04/MYD04 Level 2 files.

Modified to work for VIIRS.  Feb 2023 P. Castellanos

"""
import os, sys
from   pyobs.vx04 import Vx04_L2, MISSING, granules, SDS 
from   ffnet       import loadnet
import numpy       as     np

ALIAS = dict( TOA_NDVI = 'ndvi',
              Total_Column_Ozone = 'colO3',
              Precipitable_Water = 'water',
              )

# Translate Inputs between NNR and MODIS classes
# -----------------------------------------------
def TranslateInput(key):
    if 'mRef' in key:
        prefix = 'reflectance'
        channel = int(key[4:])
        output = (prefix,channel)
    elif 'mSre' in key:
        prefix = 'sfc_reflectance'
        channel = int(key[4:])
        output = (prefix,channel)
    else:
        output = (key,)


    return output


# Translate Targets between ANET and MODIS classes
# ------------------------------------------------
def TranslateTarget(key):
    if 'aTau' in key:
        prefix = 'aod_'
        channel = int(key[4:])
        output = (prefix, channel)
    elif 'aAE' in key:
        prefix = 'aod_'
        channel = int(key[3:])
        output = (prefix,channel)

    return output


class Vx04_NNR(Vx04_L2):
    """
    This class extends VIIRS by adding methods for producing
    NNR AOD retrievals based on the Neural Net model developed
    with class *abc_viirs*.
    """

    def __init__(self,l2_path,sat,algo,syn_time,aer_x,
                 cloud_thresh=0.70,
                 glint_thresh=40.0,
                 scat_thresh=170.0,
                 cloudFree=None,
                 aodmax=1.0,
                 coll='002',
                 nsyn=8,
                 verbose=0):
        """
        Contructs a VX04 object from VIIRS Aerosol Level 2
        granules. On input,

        l2_path  --- top directory for the VIIRS Level 2 files;
                      it must have subdirectories AERDB and/or AERDT.
        sat      --- satellite, either SNPP or NOAAXX (e.g. NOAA20)
        algo     --- Algorithm: DT_LAND, DT_OCEAN, DB_LAND or DB_OCEAN
        syn_time --- synoptic time
        aer_x    --- GEOS control file of aerosol collection

        Optional parameters:
        glint_thresh --- glint angle threshhold
        scat_thresh  --- scattering angle thresshold
        cloud_tresh  --- cloud fraction threshhold
        cloudFree    --- cloud fraction threshhold for assuring no cloud contaminations when aod is > aodmax
                        if None, no cloud free check is made
        coll         --- VIIRS data collection
        nsyn         --- number of synoptic times
              
        The following attributes are also defined:
           fractions dust, sea salt, BC+OC, sulfate
           aod_coarse
           wind
           
        It also performs Q/C, setting attribute iGood. On,
        input, *cloud_thresh* is the cloud fraction limit.
        """

        self.verbose = verbose
        self.algo    = algo
        self.cloudFree = cloudFree
        self.aodmax = aodmax
        
        # Initialize superclass
        # set anet_wav to True so MODIS wavelengths align with AERONET
        # Needed for ODS files
        # -------------------------------------------------------------
        Files = granules(l2_path,algo,sat,syn_time,coll=coll,nsyn=nsyn)
        Vx04_L2.__init__(self,Files,algo,syn_time=syn_time,nsyn=nsyn,
                              only_good=True,
                              SDS=SDS,
                              alias=ALIAS,
                              Verb=verbose,
                              anet_wav=True)            

        if self.nobs < 1:
            return # no obs, nothing to do

        # Q/C
        # ---      
        self.iGood = self.cloud<cloud_thresh  

        for i,c in enumerate(self.rChannels):
            self.iGood = self.iGood & (self.reflectance[:,i]>0)

        if "LAND" in algo:
            self.iGood = self.iGood & (self.ScatteringAngle < scat_thresh)
            for i,c in enumerate(self.sChannels):
                self.iGood = self.iGood & (self.sfc_reflectance[:,i]>0)

        if "OCEAN" in algo:
            self.iGood = self.iGood & (self.GlintAngle > glint_thresh)

        if np.any(self.iGood) == False:
            print("WARNING: Strange, no good obs left to work with")
            return

        # Create attribute for holding NNR predicted AOD
        # ----------------------------------------------
        self.aod_ =  MISSING * np.ones((self.nobs,len(self.channels)))

        # Make sure same good AOD is kept for gridding
        # --------------------------------------------
        if len(self.aod.shape) == 1:
            self.aod.shape = self.aod.shape + (1,)
        self.aod[self.iGood==False,:] = MISSING


        # Angle transforms: for NN calculations we work with cosine of angles
        # -------------------------------------------------------------------
        self.ScatteringAngle = np.cos(self.ScatteringAngle*np.pi/180.0) 
        self.RelativeAzimuth   = np.cos(self.RelativeAzimuth*np.pi/180.0)   
        self.SensorZenith    = np.cos(self.SensorZenith*np.pi/180.0)    
        self.SolarZenith     = np.cos(self.SolarZenith*np.pi/180.0)     
        self.GlintAngle      = np.cos(self.GlintAngle*np.pi/180.0)

        # Get fractional composition
        # ------------------------------
        self.speciate(aer_x,Verbose=verbose)


    def speciate(self,aer_x,Verbose=False):
        """
        Use GAAS to derive fractional composition.
        """

        self.sampleFile(aer_x,onlyVars=('TOTEXTTAU',
                                        'DUEXTTAU',
                                        'SSEXTTAU',
                                        'BCEXTTAU',
                                        'OCEXTTAU',
                                        'SUEXTTAU',
                                        ),Verbose=Verbose)

        s = self.sample
        I = (s.TOTEXTTAU<=0)
        s.TOTEXTTAU[I] = 1.E30
        self.fdu  = s.DUEXTTAU / s.TOTEXTTAU
        self.fss  = s.SSEXTTAU / s.TOTEXTTAU
        self.fbc  = s.BCEXTTAU / s.TOTEXTTAU
        self.foc  = s.OCEXTTAU / s.TOTEXTTAU
        self.fcc  = self.fbc + self.foc
        self.fsu  = s.SUEXTTAU / s.TOTEXTTAU

        # Special handle nitrate (treat it as it were sulfate)
        # ----------------------------------------------------
        try:
            self.sampleFile(aer_x,onlyVars=('NIEXTTAU',),Verbose=Verbose)
            self.fsu += self.sample.NIEXTTAU / s.TOTEXTTAU
        except:
            pass   # ignore it for systems without nitrates

        # Handle brown carbon 
        # --------------------
        try:
            self.sampleFile(aer_x,onlyVars=('BRCEXTTAU',),Verbose=Verbose)
            self.fcc += self.sample.BRCEXTTAU / s.TOTEXTTAU
        except:
            pass   # ignore it for systems without brown carbon

        del self.sample

#---
    def sampleFile(self, inFile, npzFile=None, onlyVars=None, Verbose=False):
        """
        Interpolates all variables of inFile and optionally
        save them to file *npzFile*
        """
        from gfio import GFIO, GFIOctl, GFIOHandle

        # Instantiate grads and open file
        # -------------------------------
        name, ext = os.path.splitext(inFile)
        if ext in ( '.nc4', '.nc', '.hdf'):
          fh = GFIO(inFile)     # open single file
          if fh.lm == 1:
            timeInterp = False    # no time interpolation in this case
          else:
            raise ValueError("cannot handle files with more tha 1 time, use ctl instead")
        else:
          fh = GFIOctl(inFile)  # open timeseries
          timeInterp = True     # perform time interpolation
          tymes = np.array([self.syn_time]*self.nobs)

        self.sample = GFIOHandle(inFile)
        if onlyVars is None:
            onlyVars = fh.vname

        lons = self.lon
        lats = self.lat

        

        # Loop over variables on file
        # ---------------------------
        for v in onlyVars:
            if Verbose:
                print("<> Sampling ", v)
            if timeInterp:
              var = fh.sample(v,lons,lats,tymes,Verbose=Verbose)
            else:
              var = fh.interp(v,lons,lats)
            if (var.size == 1) & (len(var.shape) == 0):
                var.shape = (1,)  #protect against when only one value is returned and shape=()
            if len(var.shape) == 1:
                self.sample.__dict__[v] = var
            elif len(var.shape) == 2:
                var = var.T # shape should be (nobs,nz)
                self.sample.__dict__[v] = var
            else:
                raise IndexError('variable <%s> has rank = %d'%(v,len(var.shape)))

        if npzFile is not None:
            savez(npzFile,**self.sample.__dict__)            


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
            iName = TranslateInput(inputName)

            if self.verbose>0:
                print('Getting NN input ',iName)

            # Retrieve input
            # --------------
            if len(iName) == 2:
                name, ch = iName
                if 'mSre' in inputName: # LAND surface reflectivity
                    k = list(self.sChannels).index(ch) # index of channel 
                elif 'mRef' in inputName: # TOA reflectances
                    k = list(self.rChannels).index(ch) # index of channel 

                input = self.__dict__[name][:,k]
                
            elif len(iName) == 1:
                name = iName[0]
                input = self.__dict__[name][:]
                
            else:
                raise ValueError("strange, len(iName)=%d"%len(iName))

            # Concatenate Inputs
            # ------------------
            if first:
                inputs = input
                first = False
            else:
                inputs = np.c_[inputs,input]

        # Keep only good observations
        # ---------------------------
        return inputs[self.iGood,:]

#--
    def apply(self,nnFile):
        """
        Apply bias correction to AOD.
        """

        # Stop here if no good obs available
        # ----------------------------------
        if self.nobs == 0:
            return # no data to work with
        if np.any(self.iGood) == False:
            return # no good data to work with

        # Load the Neural Net
        # -------------------
        self._loadNet(nnFile)

        # Evaluate NN on inputs
        # ---------------------
        targets = self.net(self._getInputs())

        # If target is angstrom exponent
        # calculate AOD
        # ------------------------------
        doAE = False
        doAEfit = False
        for targetName in self.net.TargetNames:
            if 'AEfit' in targetName:
                doAEfit = True
            elif 'AE' in targetName:
                doAE = True

        if doAEfit:
            wavs = ['440','470','550','660','870']
            wav  = np.array(wavs).astype(float)
            nwav = len(wavs)
            for i,targetName in enumerate(self.net.TargetNames):
                    if 'AEfitm' in targetName:
                        AEfitm = targets[:,i]
                    if 'AEfitb' in targetName:
                        AEfitb = targets[:,i]
            nobs = targets.shape[0]
            targets_ = np.zeros([nobs,nwav])
            targetName = []
            for i in range(nwav):
                targets_[:,i] = -1.*(AEfitm*np.log(wav[i]) + AEfitb)
                targetName.append('aTau'+wavs[i])

            targets = targets_
            self.net.TargetNames = targetName

            # Save predicted angstrom exponent
            self.ae_ = MISSING*np.ones(self.nobs)
            self.ae_[self.iGood] = AEfitm

            # calculate standard retrieval AE
            # ------------------------------------
            I = np.array(self.channels) < 900 # only visible channels, this is relevant for ocean
            aechannels = np.array(self.channels)[I]
            aodT = self.aod[:,I].T
            fit = np.polyfit(np.log(aechannels),-1.*np.log(aodT[:,self.iGood]+0.01),1)
            self.ae = MISSING*np.ones(self.nobs)
            self.ae[self.iGood] = fit[0,:]
            bad = np.isnan(self.ae)
            self.ae[bad] = MISSING

        if doAE:
            for i,targetName in enumerate(self.net.TargetNames):
                if 'Tau' in targetName:
                    name, base_wav = TranslateTarget(targetName)
                    base_wav = np.float(base_wav)
                    base_tau = targets[:,i]
                    if self.net.laod:
                        base_tau = np.exp(base_tau) - 0.01 # inverse
            for i,targetName in enumerate(self.net.TargetNames):
                if 'AE' in targetName:
                    AE = targets[:,i]
                    name, wav = TranslateTarget(targetName)
                    wav = np.float(wav)
                    data = base_tau*np.exp(-1.*AE*np.log(wav/base_wav))
                    if self.net.laod:
                        targets[:,i] = np.log(data + 0.01)
                    else:
                        targets[:,i] = data

        # Targets do not have to be in VIIRS retrieval
        # ----------------------------------------------
        for i,targetName in enumerate(self.net.TargetNames):
            name, ch = TranslateTarget(targetName)
            try:
                k = list(self.channels).index(ch) # index of channel            
            except:
                # add new target channel to end
                self.channels = np.append(self.channels,ch)
                self.aod  = np.append(self.aod,MISSING*np.ones((self.nobs,1)),axis=1)
                self.aod_ = np.append(self.aod_,MISSING*np.ones((self.nobs,1)),axis=1)

        # Replace targets with unbiased values
        # ------------------------------------
        self.channels_ = [] # channels being revised
        for i,targetName in enumerate(self.net.TargetNames):
            name, ch = TranslateTarget(targetName)
            if self.verbose>0:
                if self.net.laod:
                    print("NN Retrieving log(AOD+0.01) at %dnm "%ch)
                else:
                    print("NN Retrieving AOD at %dnm "%ch)
            k = list(self.channels).index(ch) # index of channel            
            self.channels_ = self.channels_ + [ch,]
            if self.net.laod:
                result = np.exp(targets[:,i]) - 0.01 # inverse
            else:
                result = targets[:,i]

            self.__dict__[name][self.iGood,k] = result


        # Do extra cloud filtering if required
        if self.cloudFree is not None:                 
            cloudy = (self.cloud>=self.cloudFree)
    
            contaminated = np.zeros(np.sum(self.iGood)).astype(bool)
            for targetName in self.net.TargetNames:
                name, ch = TranslateTarget(targetName)
                k = list(self.channels).index(ch) # index of channel
                result = self.__dict__[name][self.iGood,k]
                contaminated = contaminated | ( (result > self.aodmax) & cloudy[self.iGood] )
            
            icontaminated = np.arange(self.nobs)[self.iGood][contaminated]
            for targetName in self.net.TargetNames:
                name, ch = TranslateTarget(targetName)
                k = list(self.channels).index(ch) # index of channel
                self.__dict__[name][icontaminated,k] = MISSING

            if doAEfit:
                self.ae[icontaminated] = MISSING
                self.ae_[icontaminated] = MISSING

            self.iGood[icontaminated] = False





#---
        
    __call__= apply


#---

if __name__ == "__main__":

    from datetime import datetime

    l2_path = '/nobackup/VIIRS/AERDB/SNPP'
    sat     = 'SNPP'
    algo    = 'DB_LAND'
    coll    = '002'
    aer_x   = '/nobackup/NNR/Misc/tavg1_2d_aer_Nx'

    syn_time = datetime(2012,0o3,0o1,00,0,0)

    if algo == 'DB_LAND':
        nn_file = '/nobackup/NNR/Net/VIIRS/nnr_001.vsnppdbl_Tau.net'

    m = Vx04_NNR(l2_path,sat,algo,syn_time,aer_x,
                  coll=coll,
                  cloud_thresh=0.7,
                  verbose=True)

    m.apply(nn_file)
    aod = m.aod_


