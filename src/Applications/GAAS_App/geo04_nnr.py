"""
This module implements the Geostationary NNR AOD retrievals.

This version works from Dark Target Level 2 files.

"""
import os, sys
import warnings
from   pyobs.geo04 import GEO04_L2, MISSING, granules, BEST 
from   ffnet       import loadnet
from   numpy       import c_ as cat
from   numpy       import copy, ones, sin, cos, exp, arccos, pi, any, log
import numpy       as     np
from   pyobs.bits  import BITS

# SDS to be read in (usually less than in geo04)
# ----------------------------------------------
SDS = dict (
                     META = ('longitude',
                             'latitude',
                             'solar_zenith_angle',
                             'solar_azimuth_angle',
                             'sensor_zenith_angle',
                             'sensor_azimuth_angle',
                             'Scattering_Angle',
                             'Glint_Angle',
                             ),
                    LAND = ( 'Land_Ocean_Quality_Flag',
                             # u'Cloud_Pixel_Distance_Land_Ocean',
                             'Surface_Reflectance_Land',
                             'Corrected_Optical_Depth_Land',
                             'Mean_Reflectance_Land',
                             'Aerosol_Cloud_Fraction_Land',
                             ),
                    OCEAN = ('Land_Ocean_Quality_Flag',
                            #  u'Cloud_Pixel_Distance_Land_Ocean',
                             'Effective_Optical_Depth_Average_Ocean',
                            # u'Optical_Depth_Small_Average_Ocean',
                             'Aerosol_Cloud_Fraction_Ocean',
                            # u'Effective_Radius_Ocean',
                            # u'Asymmetry_Factor_Average_Ocean',
                            # u'Angstrom_Exponent_1_Ocean',
                            # u'Angstrom_Exponent_2_Ocean',
                             'Mean_Reflectance_Ocean',
                             ),
)



# Translate Inputs between NNR and MODIS classes
# -----------------------------------------------
TranslateInput = dict ( mRef470  = ('reflectance',470),
                        mRef640  = ('reflectance',640),
                        mRef860  = ('reflectance',860),
                        mRef1610 = ('reflectance',1610),
                        mRef2110 = ('reflectance',2110),
                        mSre412  = ('sfc_reflectance',412),
                        mSre470  = ('sfc_reflectance',470),
                        mSre660  = ('sfc_reflectance',660),
                        mSre2100 = ('sfc_reflectance',2100),                                    
                      )

for var in ( 'ScatteringAngle','GlintAngle',
             'SolarAzimuth', 'SolarZenith',
             'SensorAzimuth','SensorZenith',
             'cloud','qa_flag'  ):
    TranslateInput[var] = (var,)

# Translate Targets between ANET and MODIS classes
# ------------------------------------------------
TranslateTarget = dict ( aTau470 = ( 'aod_', 470 ),
                         aTau500 = ( 'aod_', 500 ),
                         aTau550 = ( 'aod_', 550 ),
                         aod_a   = ( 'aod_', 550 ),
                         aTau660 = ( 'aod_', 660 ),
                         aTau860 = ( 'aod_', 860 ),
                         )

class GEO04_NNR(GEO04_L2):
    """
    This class extends MODIS by adding methods for producing
    NNR AOD retrievals based on the Neural Net model developed
    with class *abc_c6*.
    """

    def __init__(self,l2_path,algo,syn_time,aer_x,
                 verbose=0):
        """
        Contructs a GEO04 object from GEO Aerosol Level 2
        granules. On input,

         l2_path --- top directory for the MODIS Level 2 files;
                      it must have subdirectories MOD04 and MYD04.
             algo --- aerosol algorithm: LAND, OCEAN or DEEP (for
                      Deep Blue)
        syn_time --- synoptic time

        The following attributes are also defined:
           fractions dust, sea salt, BC+OC, sulfate
           aod_coarse
           wind
           
        It also performs Q/C, setting attribute iGood.
        """

        self.verbose = verbose
        
        # Initialize superclass
        # ---------------------
        Files = granules(l2_path,syn_time)
        GEO04_L2.__init__(self,Files,algo,syn_time, SDS=SDS,                            
                               Verb=verbose)

        if self.nobs < 1:
            return # no obs, nothing to do



        # Filter observations (threshold hardwired for now)
        # -------------------------------------------------
        self.filter()

        if self.nobs < 1:
            return # no obs, nothing to do             

        # Create attribute for holding NNR predicted AOD
        # ----------------------------------------------
        self.aod_ = MISSING * ones(self.aod.shape)

        # Make sure same good AOD is kept for gridding
        # --------------------------------------------
        if len(self.aod.shape) == 1:
            self.aod.shape = self.aod.shape + (1,)
        self.aod[self.iGood==False,:] = MISSING

        # Angle transforms: for NN calculations we work with cosine of angles
        # -------------------------------------------------------------------
        self.ScatteringAngle = cos(self.ScatteringAngle*pi/180.0) 
        self.SensorAzimuth   = cos(self.SensorAzimuth*pi/180.0)   
        #self.SensorZenith    = cos(self.SensorZenith*pi/180.0)    
        #self.SolarAzimuth    = cos(self.SolarAzimuth*pi/180.0)    
        self.SolarZenith     = cos(self.SolarZenith*pi/180.0)     
        self.GlintAngle      = cos(self.GlintAngle*pi/180.0)

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
                                        'NIEXTTAU',
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
        self.fni  = s.NIEXTTAU / s.TOTEXTTAU

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
            try:
                iName = TranslateInput[inputName]
            except:
                iName = inputName

            if self.verbose>0:
                print('Getting NN input ',iName)

            # Retrieve input
            # --------------
            if type(iName) is str:
                input = self.__dict__[iName][:]

            elif len(iName) == 2:
                name, ch = iName
                if 'mSre' in inputName: # LAND or DEEP, surface reflectivity
                    k = list(self.sChannels).index(ch) # index of channel 
                elif 'mRef' in inputName: # MOD04 reflectances
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
                inputs = cat[inputs,input]

        # Keep only good observations
        # ---------------------------
        return inputs[self.iGood,:]

#--
    def apply(self,nnFile,aod_thresh=[2.,2.]):
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

        # Replace targets with NNR predicted values
        # -----------------------------------------
        self.channels_ = [] # channels being revised
        for i,targetName in enumerate(self.net.TargetNames):
            name, ch = TranslateTarget[targetName]
            if self.verbose>0:
                if self.net.laod:
                    print("NN Retrieving log(AOD+0.01) at %dnm "%ch)
                else:
                    print("NN Retrieving AOD at %dnm "%ch)
            k = list(self.aChannels).index(ch) # index of channel            
            self.channels_ += [ch,]
            if self.net.laod:
                result = exp(targets[:,i]) - 0.01 # inverse
            else:
                result = targets[:,i]

            self.__dict__[name][self.iGood,k] = result

        # Do extra cloud filtering if required
        # ------------------------------------
        if self.algo =="OCEAN":
               i = 0
        elif self.algo =="LAND":
               i = 1
        else:
               raise ValueError('Unknown algorithm <%s:'%self.algo) 
        if aod_thresh[i] is not None:
               aodBad = ((self.aod_[:,k]>aod_thresh[i])&(self.cloud>0.25))
               self.iGood = self.iGood & (aodBad==False)
               if self.verb>2:
                  print("-- Applying AOD-Cloud filter for NNR %s"%self.algo, aod_thresh[i], count_nonzero(self.iGood))
              

#---
        
    __call__= apply

#---

def _writeAll():
    """
    Write ODS files for a single, hardwired month. For testing.
    """

    l2_path = '/nobackup/1/GEO/ABI_Testing/Level2/'
    odsdir = '/nobackup/1/GEO/ABI_Testing/ODS/Y2018/M08'
    nc4dir = '/nobackup/1/GEO/ABI_Testing/Level3/Y2018/M08'
    for day in range(1,32):
        for hour in (0,3,6,9,12,15,18,21):

            from datetime import datetime
            syn_time = datetime(2018,8,day,hour,0)

            files = granules ( l2_path, syn_time, nsyn=8 )
            if len(files)==0: continue  # no granules, nothing to do.

            aer_x = '/nobackup/fp/das/Y2018/M08/f521_fp.inst1_2d_hwl_Nx.201808%02d_%02d00z.nc4'%(syn_time.day,syn_time.hour)

            print('-'*20)
            print('OCEAN', syn_time)
            nn_file = '/nobackup/NNR/Net/nnr_001.g16o_Tau.net'
            g = GEO04_NNR (l2_path,'OCEAN',syn_time,aer_x,verbose=True)
            g.apply(nn_file)
            g.writeODS(dir=odsdir,revised=True,channels=[550,],nsyn=8,Verb=1)
            g.writeg(dir=nc4dir,refine=16,channels=[550,],Verb=1)

            print('-'*20)
            print('LAND', syn_time)
            nn_file = '/nobackup/NNR/Net/nnr_001.g16l_Tau.net'
            g = GEO04_NNR (l2_path,'LAND',syn_time,aer_x,verbose=True)
            g.apply(nn_file)
            g.writeODS(dir=odsdir,revised=True,channels=[550,],nsyn=8,Verb=1)
            g.writeg(dir=nc4dir,refine=16,channels=[550,],Verb=1)

if __name__ == "__main__":

   _writeAll()


def hold():


    algo    = 'LAND'
    syn_time = datetime(2018,8,10,18,0,0)

    l2_path = '/nobackup/GEO/ABI_Testing/Level2/'
    aer_x = '/nobackup/fp/das/Y2018/M08/f521_fp.inst1_2d_hwl_Nx.201808%02d_%02d00z.nc4'%\
             (syn_time.day,syn_time.hour)

    if algo == 'OCEAN':
        nn_file = '/nobackup/NNR/Net/nnr_001.g16o_Tau.net'
    elif algo == 'LAND':
        nn_file = '/nobackup/NNR/Net/nnr_001.g16l_Tau.net'

    g = GEO04_NNR(l2_path,algo.upper(),syn_time,aer_x,
                  verbose=True)


    g.apply(nn_file)
    aod = g.aod_

    g.writeg(filename=None,dir='.',expid=None,refine=16,res=None,
           channels=[550,],Verb=1)
