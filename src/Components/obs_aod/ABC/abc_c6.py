"""
   This module implements a Neural Net based MODIS Collection 6 Neural Net Retrieval.

   Arlindo da Silva, June 2015.

"""

import os

from matplotlib.pyplot import  cm, imshow, plot, figure
from matplotlib.pyplot import  xlabel, ylabel, title, grid, savefig, legend
from matplotlib        import  ticker
from numpy             import  c_ as cat
from numpy             import  random, sort, pi, load, cos, log, std, exp
from numpy             import  reshape, arange, ones, zeros, interp
from numpy             import  meshgrid, concatenate
from scipy             import  stats, mgrid
from giant             import  MISSING, LAND, OCEAN
from nn                import NN

#..............................................................
class aodFormat(ticker.Formatter):
    def __call__(self,x,pos=None):
        y = exp(x)-0.01
        return '%4.2f'%y

#----------------------------------------------------------------------------    

class ABC_Ocean (OCEAN,NN):

    def __init__ (self,fname, Wind=None,
                  coxmunk_lut='/nobackup/NNR/Misc/coxmunk_lut.npz',
                  outliers=3., laod=True, verbose=0,
                  cloud_thresh=0.70, csvVersion = 1,
                  Input = ['mTAU550','mTAU470','mTAU660','mTAU870',
                           'ScatteringAngle','GlintAngle',
                           'SolarAzimuth','SolarZenith',
                           'SensorAzimuth','SensorZenith',
                           'cloud', 'wind'],
                  Target = [ 'aTau470', 'aTau550', 'aTau660', 'aTau870' ]):
        """
        Initializes the AOD Bias Correction (ABC) for the MODIS Ocean algorithm.

        On Input,

        fname   ---  file name for the CSV file with the co-located MODIS/AERONET
                     data (see class OCEAN)
        Wind    ---  Type of wind related parameter to be read from a NPZ file. Typical
                     values are:
                                  merra_ustar
                                  merra_wind
                     Requires a NPZ file with the data.
        outliers --  number of standard deviations for outlinear removal.
        laod    ---  if True, targets are log-transformed AOD, log(Tau+0.01)
        """

        self.Input = Input
        self.Target = Target
        self.verbose = verbose
        self.laod = laod
        self.Wind = Wind

        OCEAN.__init__(self,fname,csvVersion=csvVersion) # initialize superclass

        # Read in wind if desired
        # ------------------------
        if Wind != None:
            if csvVersion==1:
                self.wind = load(self.ident + "_" + Wind + ".npz")['var']
            else:
                self.wind = load(self.ident + "2_" + Wind + ".npz")['var']
        else:
            self.wind = zeros(self.N)
            print "WARNING: No longer using *ncep_windspd* because of too many undefs" 
            print "WARNING: wind attribute being set to zero"

        # Define wind speed dependent ocean albedo
        # ----------------------------------------
        self.getCoxMunk(coxmunk_lut)

        # Q/C
        # ---
        self.iValid =  (self.qa>0) & \
                      (self.aTau470 > 0.0) & (self.aTau550 > 0.0) &  \
                      (self.aTau660 > 0.0) & (self.aTau870 > 0.0) &  \
                      (self.mTau470 > 0.0) & (self.mTau550 > 0.0) &  \
                      (self.mTau660 > 0.0) & (self.mTau870 > 0.0) &  \
                      (self.mtau470 > 0.0) & (self.mtau550 > 0.0) &  \
                      (self.mtau660 > 0.0) & (self.mtau870 > 0.0) &  \
                      (self.cloud <cloud_thresh) & (self.wind>=0.0) &  \
                      (self.GlintAngle != MISSING )

        # Outlier removal based on log-transformed AOD
        # --------------------------------------------
        if outliers > 0.:
            d = log(self.mTau550+0.01) - log(self.aTau550+0.01)
            dg = d[self.iValid]
            if self.verbose>0:
                print "Outlier removal: %d   sig_d = %f  nGood=%d "%(-1,std(dg),dg.size)
            for iter in range(3):
                self.iValid = self.iValid & (abs(d)<outliers*std(d[self.iValid]))
                if self.verbose>0:
                    dg = d[self.iValid]
                    print "Outlier removal: %d   sig_d = %f  nGood=%d "%(iter,std(dg),dg.size)
            
        # Angle transforms: for NN work we work with cosine of angles
        # -----------------------------------------------------------
        self.ScatteringAngle = cos(self.ScatteringAngle*pi/180.0) 
        self.SensorAzimuth   = cos(self.SensorAzimuth*pi/180.0)   
        self.SensorZenith    = cos(self.SensorZenith*pi/180.0)    
        self.SolarAzimuth    = cos(self.SolarAzimuth*pi/180.0)    
        self.SolarZenith     = cos(self.SolarZenith*pi/180.0)     
        self.GlintAngle      = cos(self.GlintAngle*pi/180.0)      

#----------------------------------------------------------------------------    

class ABC_Land (LAND,NN):

    def __init__ (self, fname,
                  Albedo='albedo',
                  alb_min = 0.25,
                  outliers=3.,
                  laod=True,
                  verbose=0,
                  cloud_thresh=0.70,
                  csvVersion = 1,
                  Input = ['mTau550','mTau470','mTau660',
                           'mTAU550','mTAU470','mTAU660',
                           'ScatteringAngle',
                           'SolarAzimuth','SolarZenith',
                           'SensorAzimuth','SensorZenith',
                           'cloud' ],
                  Target = [ 'aTau470', 'aTau550', 'aTau660' ]):
        """
        Initializes the AOD Bias Correction (ABC) for the MODIS Land algorithm.

        On Input,

        fname   ---  file name for the CSV file with the co-located MODIS/AERONET
                     data (see class OCEAN)

        Albedo  ---  albedo file name identifier; albedo file will be created
                     from this identifier (See below).
        outliers --  number of standard deviations for outlinear removal.
        laod    ---  if True, targets are log-transformed AOD, log(Tau+0.01)
        """

        if csvVersion != 1:
            raise ValueError, 'must use CVS Version 1 for land because of QA flags'

        self.Input = Input
        self.Target = Target
        self.verbose = verbose
        self.laod = laod

        LAND.__init__(self,fname,csvVersion=csvVersion) # initialize superclass

        # Read in wind if desired
        # ------------------------
        self.albedo = load(self.ident + "_" + Albedo + ".npz")['var']

        # Q/C: enforce QA=3 and albedo in (0,0.25), scattering angle<170
        # --------------------------------------------------------------
        self.iValid =  (self.qa==3)                & \
                      (self.aTau470 > -0.01)      & \
                      (self.aTau550 > -0.01)      & \
                      (self.aTau660 > -0.01)      & \
                      (self.mTau470 > -0.01)      & \
                      (self.mTau550 > -0.01)      & \
                      (self.mTau660 > -0.01)      & \
                      (self.mTau2100> -0.01)      & \
                      (self.mSre470 >  0.0)       & \
                      (self.mSre660 >  0.0)       & \
                      (self.mSre2100>  0.0)       & \
                      (self.cloud<cloud_thresh)   & \
                      (self.ScatteringAngle<170.) & \
                      (self.albedo>0)             & \
                      (self.albedo<alb_min)

        print self.qa[self.iValid].shape
        
        # Outlier removal based on log-transformed AOD
        # --------------------------------------------
        if outliers > 0.:
            d = log(self.mTau550+0.01) - log(self.aTau550+0.01)
            dg = d[self.iValid]
            if self.verbose>0:
                print "Outlier removal: %d   sig_d = %f  nGood=%d "%\
                      (-1,std(dg),dg.size)
            for iter in range(3):
                self.iValid = self.iValid & (abs(d)<outliers*std(d[self.iValid]))
                if self.verbose>0:
                    dg = d[self.iValid]
                    print "Outlier removal: %d   sig_d = %f  nGood=%d "\
                          %(iter,std(dg),dg.size)
            
        # Angle transforms: for NN work we work with cosine of angles
        # -----------------------------------------------------------
        self.ScatteringAngle = cos(self.ScatteringAngle*pi/180.0) 
        self.SensorAzimuth   = cos(self.SensorAzimuth*pi/180.0)   
        self.SensorZenith    = cos(self.SensorZenith*pi/180.0)    
        self.SolarAzimuth    = cos(self.SolarAzimuth*pi/180.0)    
        self.SolarZenith     = cos(self.SolarZenith*pi/180.0)     
        self.GlintAngle      = cos(self.GlintAngle*pi/180.0)      

    def getAlbedo(self,npz_file):
        from grads import GrADS
        ga = GrADS(Echo=False,Window=False)
        ga('open albedo_clim.ctl')
        self.addVar(ga,npz_file,expr='albedo',clmYear=2000)

#---------------------------------------------------------------------------------
def _cat2 (X, Y):
    """
    Given 2 arrays of same shape, returns array of shape (2,N),
    where N = X.size = Y.size
    """
    xy = concatenate((X.ravel(),Y.ravel())) # shape is (N+N)
    return reshape(xy,(2,X.size))         # shape is (2,N)

def _plotKDE(x_values,y_values,x_bins=None,y_bins=None,
             x_label='AERONET', y_label='MODIS',formatter=None):
        """
        Plot Target vs Model using a 2D Kernel Density Estimate.
        """

        if x_bins == None: x_bins = arange(-5., 1., 0.1 )
        if y_bins == None: y_bins = x_bins

        Nx = len(x_bins)
        Ny = len(y_bins)

        print "Evaluating 2D kernel on grid with (Nx,Ny)=(%d,%d) ..."%(Nx,Ny)
        kernel = stats.kde.gaussian_kde(_cat2(x_values,y_values))
        X, Y = meshgrid(x_bins,y_bins)   # each has shape (Ny,Nx)
        Z = kernel(_cat2(X,Y))           # shape is (Ny*Nx)
        Z = reshape(Z,X.shape)

        fig = figure()
        # ax = fig.add_axes([0.1,0.1,0.75,0.75])
        ax = fig.add_axes([0.1,0.1,0.75,0.75])
        if formatter != None:
            ax.xaxis.set_major_formatter(formatter)
            ax.yaxis.set_major_formatter(formatter)
        imshow(Z, cmap=cm.gist_earth_r, origin='lower', 
               extent=(x_bins[0],x_bins[-1],y_bins[0],y_bins[-1]) )
        plot([x_bins[0],x_bins[-1]],[y_bins[0],y_bins[-1]],'k')
        xlabel(x_label)
        ylabel(y_label)
        grid()

#--------------------------------------------------------------------------------------

def _remove1():

    Input_all = ['mTau550','mTau470','mTau660','mTau870',
                 'mTAU550','mTAU470','mTAU660','mTAU870',
                 'ScatteringAngle','GlintAngle',
                 'SolarAzimuth', 'SolarZenith',
                 'SensorAzimuth','SensorZenith',
                 'cloud', 'ustar']

    Target=['aTau550',]

    mydo = ABC_Ocean('SUPER_ocean.Aqua.csv',ustar=False,verbose=1)
    mydo.split()

    for i in [-1,]+range(len(Input_all)):

        print "------------------------------------------------------------------------------"
        Input = Input_all[:] # make a copy of it
        if i<0:
            print "--> Excluding: (nothing)"
        else:
            print "--> Excluding: ", Input_all[i]
            del Input[i]         # delete ith item

        nHidden = len(Input)

        print "--> nHidden = ", nHidden
        print "-->  Inputs = ", Input
        
        mydo.train(Input=Input,Target=Target,nHidden=nHidden)
        out, reg = mydo.test()

#--------------------------------------------------------------------------------------

def _testOcean(filename):

    Input_nnc = ['mTau550','mTau470','mTau660','mTau870',
                 'mTAU550','mTAU470','mTAU660','mTAU870',
                 'ScatteringAngle','GlintAngle',
                 'SolarAzimuth', 'SolarZenith',
                 'SensorAzimuth','SensorZenith',
                 'cloud', 'wind' ]

    Input_nnr1 = ['mRef470','mRef550','mRef660', 'mRef870',
                 'mRef1200','mRef1600','mRef2100',
                 'ScatteringAngle', 'GlintAngle',
                 'SolarAzimuth', 'SolarZenith',
                 'SensorAzimuth','SensorZenith',
                 'cloud', 'wind' ]

    Input_nnr2 = ['mRef470','mRef550','mRef660', 'mRef870',
                 'mRef1200','mRef1600','mRef2100',
                 'ScatteringAngle', 'GlintAngle',
                 'SolarAzimuth', 'SolarZenith',
                 'SensorAzimuth','SensorZenith',
                 'cloud', 'albedo' ]

    Input_min = ['mTau550','mTAU550',
                 'GlintAngle', 'cloud', 'wind' ]

#    Target=['aTau_c',]
    Target=['aTau550',]
#    Target = [ 'aTau550','aTau470','aTau660', 'aTau870' ]

    # Read and split dataset in training/testing subsets
    # --------------------------------------------------
    mxdo = ABC_Ocean(filename,Wind='merra_wind',
                     verbose=1,csvVersion=2)
    mxdo.split()
    ident = mxdo.ident

    expid = 'nnr_002'
    for Input in (Input_nnr2,):

        nHidden = len(Input)
        topology = (len(Input), nHidden, len(Target))
        
        print "-"*80
        print "--> nHidden = ", nHidden
        print "-->  Inputs = ", Input
        
        mxdo.train(Input=Input,Target=Target,nHidden=nHidden,topology=topology)
        out, reg = mxdo.test()

        mxdo.savenet(expid+"."+ident+'_Tau.net')

        # Plot KDE of corrected AOD
        # -------------------------
        mxdo.plotKDE(figfile=expid+"."+ident+"_kde-"+Target[0][1:]+"-corrected.png")

        # Plot KDE of uncorrected AOD
        # ---------------------------
        targets = mxdo.getTargets(mxdo.iValid).squeeze()
        original = log(mxdo.mTau550[mxdo.iValid]+0.01)
        _plotKDE(targets,original,y_label='Original MODIS')
        title("Log("+Target[0][1:]+"+0.01)- "+ident)
        savefig(expid+"."+ident+"_kde-"+Target[0][1:]+'.png')

        # Scatter diagram for testing
        # ---------------------------
        mxdo.plotScat(figfile=expid+"."+ident+"_scat-"+Target[0][1:]+'.png')

    return mxdo

#---------------------------------------------------------------------
def _testLand(filename):

    Input_nnc = ['mTau550','mTau470','mTau660', 'mTau2100',
                 'mSre470','mSre660', 'mSre2100',
                 'ScatteringAngle',
                 'SolarAzimuth', 'SolarZenith',
                 'SensorAzimuth','SensorZenith',
                 'cloud', 'albedo' ]

    Input_nnr1 = ['mRef550','mRef470','mRef660', 'mRef870',
                 'mRef1200','mRef1600','mRef2100',
#                 'mSre470','mSre660', 'mSre2100',
                 'ScatteringAngle',
                 'SolarAzimuth', 'SolarZenith',
                 'SensorAzimuth','SensorZenith',
                 'cloud', 'albedo' ]

    Input_nnr2 = ['mRef550','mRef470','mRef660', 'mRef870',
                 'mRef1200','mRef1600','mRef2100',
                 'mSre470','mSre660', 'mSre2100',
                 'ScatteringAngle',
                 'SolarAzimuth', 'SolarZenith',
                 'SensorAzimuth','SensorZenith',
                 'cloud', 'albedo' ]

    Input_min = ['mTau550',
                 'ScatteringAngle', 'cloud', 'albedo' ]

#    Target=['aTau_c',]
    Target=['aTau550',]
#    Target = [ 'aTau550', 'aTau470', 'aTau660', ]


    # Read and split dataset in training/testing subsets
    # --------------------------------------------------
    mxdl = ABC_Land(filename,alb_min=0.25,verbose=1,csvVersion=1)
    mxdl.split()

    ident = mxdl.ident
    expid = 'nnr_002'
    for Input in (Input_nnr2,):

        nHidden = len(Input)
        
        print "-"*80
        print "--> nHidden = ", nHidden
        print "-->  Inputs = ", Input
        
        mxdl.train(Input=Input,Target=Target,nHidden=nHidden)
        out, reg = mxdl.test()

        mxdl.savenet(expid+"."+ident+'_Tau.net')

        # Plot KDE of corrected AOD
        # -------------------------
        mxdl.plotKDE(figfile=expid+"."+ident+"_kde-"+Target[0][1:]+"-corrected.png")

        # Plot KDE of uncorrected AOD
        # ---------------------------
        targets = mxdl.getTargets(mxdl.iValid).squeeze()
        original = log(mxdl.mTau550[mxdl.iValid]+0.01)
        _plotKDE(targets,original,y_label='Original MODIS')
        title("Log("+Target[0][1:]+"+0.01)- "+ident)
        savefig(expid+"."+ident+"_kde-"+Target[0][1:]+'.png')

        # Scatter diagram for testing
        # ---------------------------
        mxdl.plotScat(figfile=expid+"."+ident+"_scat-"+Target[0][1:]+'.png')

    return mxdl

#---------------------------------------------------------------------
def _svrLand(filename):
    from sklearn.svm import SVR
    
    Input_nnc = ['mTau550','mTau470','mTau660', 'mTau2100',
                 'mSre470','mSre660', 'mSre2100',
                 'ScatteringAngle',
                 'SolarAzimuth', 'SolarZenith',
                 'SensorAzimuth','SensorZenith',
                 'cloud', 'albedo' ]

    Input_nnr1 = ['mRef550','mRef470','mRef660', 'mRef870',
                 'mRef1200','mRef1600','mRef2100',
#                 'mSre470','mSre660', 'mSre2100',
                 'ScatteringAngle',
                 'SolarAzimuth', 'SolarZenith',
                 'SensorAzimuth','SensorZenith',
                 'cloud', 'albedo' ]

    Input_nnr2 = ['mRef550','mRef470','mRef660', 'mRef870',
                 'mRef1200','mRef1600','mRef2100',
                 'mSre470','mSre660', 'mSre2100',
                 'ScatteringAngle',
                 'SolarAzimuth', 'SolarZenith',
                 'SensorAzimuth','SensorZenith',
                 'cloud', 'albedo' ]

    Input_min = ['mTau550',
                 'ScatteringAngle', 'cloud', 'albedo' ]

#    Target=['aTau_c',]
    Target=['aTau550',]
#    Target = [ 'aTau550', 'aTau470', 'aTau660', ]


    # Read and split dataset in training/testing subsets
    # --------------------------------------------------
    mxdl = ABC_Land(filename,alb_min=0.25,verbose=1,csvVersion=1)
    mxdl.split()

    ident = mxdl.ident
    expid = 'svr_002'
    for Input in (Input_nnr2,):

        X = mxdl.getInputs(I=mxdl.iValid,Input=Input)
        y = log(mxdl.aTau550[mxdl.iValid]+0.01)

        print X.shape, y.shape

        # SVR fitting
        # -----------
        svr = SVR(kernel='rbf', C=1e3, gamma=0.1,verbose=True)
        svr = svr.fit(X, y)

        # Plot KDE of corrected AOD
        # -------------------------
        targets = y
        results = svr.predict(X)
        _plotKDE(targets,results,y_label='SVR Fit')
        title("SVR Log("+Target[0][1:]+"+0.01)- "+ident)
        savefig(expid+"."+ident+"_kde-"+Target[0][1:]+'.png')

        # Scatter diagram for testing
        # ---------------------------
        #mxdl.plotScat(figfile=expid+"."+ident+"_scat-"+Target[0][1:]+'.png')

    return mxdl

def doAlbedo():
    from anet import LAND
    from grads import GrADS
    ga = GrADS(Echo=False,Window=False)
    ga('open albedo_clim.ctl')

    modl = LAND('SUPER_land.Terra.csv')
    modl.addVar(ga,'modl_albedo.npz',expr='albedo',clmYear=2000)
        
    mydl = LAND('SUPER_land.Aqua.csv')
    mydl.addVar(ga,'mydl_albedo.npz',expr='albedo',clmYear=2000)
        
def doWind():
    from anet import OCEAN
    from grads import GrADS
    ga = GrADS(Echo=False,Window=False)

    ga('xdfopen merra_slv-hourly.ddf')

    modo = OCEAN('SUPER2_combo.Terra.csv',)
    modo.addVar(ga,'modo2_merra_wind.npz',expr='mag(u10m,v10m)',vname='wind')

#    mydo = OCEAN('SUPER2_combo.Aqua.csv')
#    mydo.addVar(ga,'mydo2_merra_wind.npz',expr='mag(u10m,v10m)',vname='wind')

        
#------------------------------------------------------------------
  
if __name__ == "__main__":

    modl = _svrLand('SUPER_land.Terra.csv')

    # modo = _testOcean('SUPER2_combo.Terra.csv')
    # mydo = _testOcean('SUPER2_combo.Aqua.csv')
    # mydl = _testLand('SUPER_land.Aqua.csv')
    # modl = _testLand('SUPER_land.Terra.csv')

def hold():

    doWind()
    doAlbedo()

    mxdx = _testOcean('SUPER2_combo.Terra.csv')

    mydo = _testOcean('SUPER2_combo.Aqua.csv')
        
    mxdl = _testLand('SUPER_land.Aqua.csv')

#    mxdx = _testOcean('SUPER2_combo.Aqua.csv')
    mxdx = _testOcean('SUPER2_combo.Terra.csv')

    mxdx = _testLand('SUPER_land.Aqua.csv')

    
