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
        channel = int(key[key.find('Ref')+3:])
        output = (prefix,channel)
    elif 'mSre' in key:
        prefix = 'sfc_reflectance'
        channel = int(key[key.find('Sre')+3:])
        output = (prefix,channel)
    elif 'mTau' in key:
        prefix = 'aod'
        channel = int(key[key.find('Tau')+3:])
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
                 aodmax=2.0,
                 aodSTD=3.0,
                 aodLength=0.5,
                 wavs=['440','470','550','660','870'],                 
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
        aodSTD      --- number of standard deviations for checking for outliers
        aodLength   --- length scale (degrees) to look for outliers                        
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
        self.aodSTD = aodSTD
        self.aodLength = aodLength
        if type(wavs) is str:
            self.wavs = wavs.split(',')
        else:
            self.wavs = wavs        
        
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
        if "pixel_elevation" in self.__dict__:
            self.pixel_elevation = self.pixel_elevation*1e-4

        if self.nobs < 1:
            return # no obs, nothing to do

        # Q/C
        # ---      
        self.iGood = self.cloud<cloud_thresh  

        for i,c in enumerate(self.rChannels):
            self.iGood = self.iGood & (self.reflectance[:,i]>0)

        if ("LAND" in algo) or ("DEEP" in algo):
            self.iGood = self.iGood & (self.ScatteringAngle < scat_thresh)
            if "LAND" in algo:
                # 412 surface reflectance not used for vegetated surfaces
                self.iGood = self.iGood & self.sfc_reflectance[:,0].mask & ~self.sfc_reflectance[:,1].mask & ~self.sfc_reflectance[:,2].mask

            elif "DEEP" in algo:
                for i,c in enumerate(self.sChannels):
                    self.iGood = self.iGood & ~self.sfc_reflectance[:,i].mask

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
        s.TOTEXTTAU[I] = 1.E-30
        self.fdu  = s.DUEXTTAU / s.TOTEXTTAU
        self.fss  = s.SSEXTTAU / s.TOTEXTTAU
        self.fbc  = s.BCEXTTAU / s.TOTEXTTAU
        self.foc  = s.OCEXTTAU / s.TOTEXTTAU
        self.fcc  = self.fbc + self.foc
        self.fsu  = s.SUEXTTAU / s.TOTEXTTAU

        for spc in ['fdu','fss','fbc','foc','fcc','fss']:
            i = np.isnan(self.__dict__[spc])
            self.__dict__[spc][i] = 0.0

            i = np.isinf(self.__dict__[spc])
            self.__dict__[spc][i] = 0.0

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
                elif 'mTau' in inputName: # Predicted AOD
                    k = list(self.channels).index(ch) 

                if inputName[0] == 'l':
                    feature = self.__dict__[name][:,k]
                    input = self.net.__dict__['scaler_'+inputName].transform(feature.reshape(-1,1)).squeeze()
                else:
                    input = self.__dict__[name][:,k]
                
            elif len(iName) == 1:
                name = iName[0]
                if inputName[0] == 'l':
                    feature = self.__dict__[name[1:]][:]
                    input = self.net.__dict__['scaler_'+inputName].transform(feature.reshape(-1,1)).squeeze()
                else:
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
        if hasattr(self.net,"scale"):
            if self.net.scale:
                if len(self.net.TargetNames) == 1:
                    targets = self.net.scaler.inverse_transform(targets.reshape(-1,1)).squeeze()
                else:
                    targets = self.net.scaler.inverse_transform(targets)        

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
            wav  = np.array(self.wavs).astype(float)
            nwav = len(self.wavs)
            AEfitb = None
            for i,targetName in enumerate(self.net.TargetNames):
                    if 'AEfitm' in targetName:
                        AEfitm = targets[:,i]
                    if 'AEfitb' in targetName:
                        AEfitb = targets[:,i]
                    if 'aTau550' in targetName:
                        tau550 = targets[:,i]

            if AEfitb is None:
                AEfitb = -1.*(tau550 + AEfitm*np.log(550.))
            nobs = targets.shape[0]
            targets_ = np.zeros([nobs,nwav])
            targetName = []
            for i in range(nwav):
                targets_[:,i] = -1.*(AEfitm*np.log(wav[i]) + AEfitb)
                targetName.append('aTau'+self.wavs[i])

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
            iIndex = np.arange(len(self.iGood))[self.iGood]
            aodT = aodT[:,iIndex] + 0.01
            mask = aodT.min(axis=0) > 0
            posIndex = iIndex[mask]            
            fit = np.polyfit(np.log(aechannels),-1.*np.log(aodT[:,mask]+0.01),1)
            self.ae = MISSING*np.ones(self.nobs)
            self.ae[posIndex] = fit[0,:]


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
            # start by checking the cloud masks
            cloudy = (self.cloud>=self.cloudFree)
    
            contaminated = np.zeros(np.sum(self.iGood)).astype(bool)
            for targetName in self.net.TargetNames:
                name, ch = TranslateTarget(targetName)
                k = list(self.channels).index(ch) # index of channel
                result = self.__dict__[name][self.iGood,k]
                contaminated = contaminated | ( (result > self.aodmax) & cloudy[self.iGood] )
            
            icontaminated = np.arange(self.nobs)[self.iGood][contaminated]

            if self.verbose:
                print('Filtering out ',np.sum(contaminated),' suspected cloud contaminated pixels')


            for targetName in self.net.TargetNames:
                name, ch = TranslateTarget(targetName)
                k = list(self.channels).index(ch) # index of channel
                self.__dict__[name][icontaminated,k] = MISSING

            if doAEfit:
                self.ae_[icontaminated] = MISSING

            self.iGood[icontaminated] = False

            # check for outliers
            # start with highest AOD550 value
            # find all the pixels within a 1 degree neighborhood
            # check if it is outside of mean + N*sigma of the other pixels
            # aodSTD parameter is equal to N
            # continue until no outliers are found
            find_outliers = True
            k = list(self.channels).index(550)
            aod550 = np.ma.array(self.aod_[self.iGood,k])
            aod550.mask = np.zeros(len(aod550)).astype(bool)
            Lon = self.Longitude[self.iGood]
            Lat = self.Latitude[self.iGood]
            gIndex = np.arange(self.nobs)[self.iGood]
            iOutliers = []
            count = 0
            while find_outliers & (count<len(aod550)):
                maxaod = aod550.max()
                imax   = np.argmax(aod550)
                aod550.mask[imax] = True
                lon = Lon[imax]
                lat = Lat[imax]

                # find the neighborhood of pixels
                iHood = (Lon<=lon+self.aodLength) & (Lon>=lon-self.aodLength) & (Lat<=lat+self.aodLength) & (Lat>=lat-self.aodLength)
                if (np.sum(iHood) <= 1) & (maxaod > self.aodmax):
                    #this pixel has no neighbors and is high. Filter it.
                    iOutliers.append(gIndex[imax])
                else:
                    aodHood = aod550[iHood]
                    if maxaod > (aodHood.mean() + self.aodSTD*aodHood.std()):
                        iOutliers.append(gIndex[imax])
                    else:
                        find_outliers = False  # done looking for outliers
                count +=1
            if self.verbose:
                print("Filtering out ",len(iOutliers)," outlier pixels")

            self.iOutliers = iOutliers

            if len(iOutliers) > 0:
                for targetName in self.net.TargetNames:
                    name, ch = TranslateTarget(targetName)
                    k = list(self.channels).index(ch) # index of channel
                    self.__dict__[name][iOutliers,k] = MISSING

                if doAEfit:
                    self.ae_[iOutliers] = MISSING

                self.iGood[iOutliers] = False



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


