"""
   Classes for reading Rob Levy's CSV files and create classes for
   each sat/algorithm type.
"""

MISSING = -99.99

from numpy import loadtxt, ones, savez, pi, cos, sin, arccos, zeros, interp

from pyobs.npz import NPZ

Months = ('jan','feb','mar','apr','may','jun','jul','aug','sep','oct','nov','dec')

def date2nymd(s):
    S = s.split('-')
    yy = int(S[0]) 
    mm = int(S[1])
    dd = int(S[2])
    return (yy,mm,dd)

def time2nhms(s):
    S = s.split(':')
    h = int(S[0])
    m = int(S[1])
    return (h,m)

def gatime(date,time):
    yy,mm,dd = date2nymd(date)
    return "%sZ%d%s%d"%(time,dd,Months[mm-1],yy)

def _getCols1():
    """
    Robs Version 1 CSV files
    """
    cols = 'Date,DOY,Time,Location,Satellite,Collection,Longitude,Latitude,SolarZenith,SolarAzimuth,SensorZenith,SensorAzimuth,ScatteringAngle,nval_AOT_1020_l20,mean_AOT_1020_l20,mean_AOT_870_l20,mean_AOT_675_l20,sdev_AOT_675_l20,mean_AOT_500_l20,mean_AOT_440_l20,mean_AOT_380_l20,mean_AOT_340_l20,mean_Water_cm_l20,nval_AOT_1020_l15,mean_AOT_1020_l15,mean_AOT_870_l15,mean_AOT_675_l15,sdev_AOT_675_l15,mean_AOT_500_l15,mean_AOT_440_l15,mean_AOT_380_l15,mean_AOT_340_l15,mean_Water_cm_l15,npix_AOT0550,mean_AOT0550,sdev_AOT0550,mean_rAOTse0550,sdev_rAOTse0550,mean_AOT0470corr_l,npix_AOT0550corr_l,pval_AOT0550corr_l,mean_AOT0550corr_l,sdev_AOT0550corr_l,mean_AOT0660corr_l,mean_AOT2100corr_l,mean_rAOTse0550_l,pval_rAOTse0550_l,mean_AOT0550sm_l,pval_AOT0550sm_l,mean_Aexp0470_0670_l,mean_surfre0470_l,mean_surfre0660_l,mean_surfre2100_l,mean_fiterr_l,mean_atype_l,mean_cfrac_l,mean_mconc_l,QA0470_l,mean_mref0470_l,mean_mref0550_l,mean_mref0660_l,mean_mref0870_l,mean_mref1200_l,mean_mref1600_l,mean_mref2100_l,pval_mref0470_l,pval_mref0550_l,pval_mref0660_l,pval_mref0870_l,pval_mref1200_l,pval_mref1600_l,pval_mref2100_l,mean_AOT0470ea_o,npix_AOT0550ea_o,pval_AOT0550ea_o,mean_AOT0550ea_o,sdev_AOT0550ea_o,mean_AOT0660ea_o,mean_AOT0870ea_o,mean_AOT1200ea_o,mean_AOT1600ea_o,mean_AOT2100ea_o,mean_AOT0470sa_o,npix_AOT0550sa_o,pval_AOT0550sa_o,mean_AOT0550sa_o,sdev_AOT0550sa_o,mean_AOT0660sa_o,mean_AOT0870sa_o,mean_AOT1200sa_o,mean_AOT1600sa_o,mean_AOT2100sa_o,mean_rAOTse0550a_o,mean_effr0550a_o,sdev_effr0550a_o,mean_solindx_sa_o,mean_solindx_la_o,mean_lsqerr_a_o,mean_cfrac_o,sdev_cfrac_o,QAavg_o,mean_mref0470_o,mean_mref0550_o,mean_mref0660_o,mean_mref0870_o,mean_mref1200_o,mean_mref1600_o,mean_mref2100_o,sdev_mref0470_o,sdev_mref0550_o,sdev_mref0660_o,sdev_mref0870_o,sdev_mref1200_o,sdev_mref1600_o,sdev_mref2100_o,mean_wni,mean_wir,pval_wni,pval_wir,mean_pathrad0470_l,mean_pathrad0660_l,mean_critref0470_l,mean_critref0660_l,mean_errprad0470_l,mean_errprad0660_l,mean_errcref0470_l,mean_errcref0660_l,mean_qwtprad0470_l,mean_qwtprad0660_l,mean_qwtcref0470_l,mean_qwtcref0660_l,npix_AOT0550dpbl_l,pval_AOT0550dpbl_l,mean_AOT0550dpbl_l,sdev_AOT0550dpbl_l,mean_AOT0412dpbl_l,mean_AOT0470dpbl_l,mean_AOT0660dpbl_l,mean_Aext0412_0470dpbl_l,mean_SSA0412dpbl_l,mean_SSA0470dpbl_l,mean_SSA0660dpbl_l,mean_surfre0412dpbl_l,mean_surfre0470dpbl_l,mean_surfre0660dpbl_l,tau_550_norm,eta_norm,tau_f,tau_c,alpha_norm,alpha_f,Deta,tau_466,tau_553,tau_644,tau_866,tau_2119,Angs_466_644,exp_errorO_pct,exp_errorL_pct,ncep_pwat,ncep_O3,ncep_pres,ncep_windspd,ncep_winddir'
    return cols

def _getCols2():
    """
    Robs Version 2 CSV files
    """
    cols='Date,DOY,Time,Location,Satellite,Collection,Longitude,Latitude,SolarZenith,SolarAzimuth,SensorZenith,SensorAzimuth,ScatteringAngle,nval_AOT_1020_l20,mean_AOT_1020_l20,mean_AOT_870_l20,mean_AOT_675_l20,sdev_AOT_675_l20,mean_AOT_500_l20,mean_AOT_440_l20,mean_AOT_380_l20,mean_AOT_340_l20,mean_Water_cm_l20,nval_AOT_1020_l15,mean_AOT_1020_l15,mean_AOT_870_l15,mean_AOT_675_l15,sdev_AOT_675_l15,mean_AOT_500_l15,mean_AOT_440_l15,mean_AOT_380_l15,mean_AOT_340_l15,mean_Water_cm_l15,npix_AOT0550,mean_AOT0550,sdev_AOT0550,mean_rAOTse0550,sdev_rAOTse0550,mean_AOT0470corr_l,npix_AOT0550corr_l,pval_AOT0550corr_l,mean_AOT0550corr_l,sdev_AOT0550corr_l,mean_AOT0660corr_l,mean_AOT2100corr_l,mean_rAOTse0550_l,pval_rAOTse0550_l,mean_AOT0550sm_l,pval_AOT0550sm_l,mean_Aexp0470_0670_l,mean_surfre0470_l,mean_surfre0660_l,mean_surfre2100_l,mean_fiterr_l,mean_atype_l,mean_cfrac_l,mean_mconc_l,QAdark_l,mean_mref0470_l,mean_mref0550_l,mean_mref0660_l,mean_mref0870_l,mean_mref1200_l,mean_mref1600_l,mean_mref2100_l,pval_mref0470_l,pval_mref0550_l,pval_mref0660_l,pval_mref0870_l,pval_mref1200_l,pval_mref1600_l,pval_mref2100_l,mean_AOT0470ea_o,npix_AOT0550ea_o,pval_AOT0550ea_o,mean_AOT0550ea_o,sdev_AOT0550ea_o,mean_AOT0660ea_o,mean_AOT0870ea_o,mean_AOT1200ea_o,mean_AOT1600ea_o,mean_AOT2100ea_o,mean_AOT0470sa_o,npix_AOT0550sa_o,pval_AOT0550sa_o,mean_AOT0550sa_o,sdev_AOT0550sa_o,mean_AOT0660sa_o,mean_AOT0870sa_o,mean_AOT1200sa_o,mean_AOT1600sa_o,mean_AOT2100sa_o,mean_rAOTse0550a_o,mean_effr0550a_o,sdev_effr0550a_o,mean_solindx_sa_o,mean_solindx_la_o,mean_lsqerr_a_o,mean_cfrac_o,sdev_cfrac_o,QAavg_o,mean_mref0470_o,mean_mref0550_o,mean_mref0660_o,mean_mref0870_o,mean_mref1200_o,mean_mref1600_o,mean_mref2100_o,sdev_mref0470_o,sdev_mref0550_o,sdev_mref0660_o,sdev_mref0870_o,sdev_mref1200_o,sdev_mref1600_o,sdev_mref2100_o,mean_wni,mean_wir,pval_wni,pval_wir,mean_pathrad0470_l,mean_pathrad0660_l,mean_critref0470_l,mean_critref0660_l,mean_errprad0470_l,mean_errprad0660_l,mean_errcref0470_l,mean_errcref0660_l,mean_qwtprad0470_l,mean_qwtprad0660_l,mean_qwtcref0470_l,mean_qwtcref0660_l,QAdpbl_l,npix_AOT0550dpbl_l,pval_AOT0550dpbl_l,mean_AOT0550dpbl_l,sdev_AOT0550dpbl_l,mean_AOT0412dpbl_l,mean_AOT0470dpbl_l,mean_AOT0660dpbl_l,mean_Aext0412_0470dpbl_l,mean_SSA0412dpbl_l,mean_SSA0470dpbl_l,mean_SSA0660dpbl_l,mean_surfre0412dpbl_l,mean_surfre0470dpbl_l,mean_surfre0660dpbl_l,tau_550_norm,eta_norm,tau_f,tau_c,alpha_norm,alpha_f,Deta,tau_466,tau_553,tau_644,tau_866,tau_2119,Angs_466_644,exp_errorO_pct,exp_errorL_pct,ncep_pwat,ncep_O3,ncep_pres,ncep_windspd,ncep_winddir,eta_MODIS, tau_f_MODIS, tau_c_MODIS, alpha_MODIS'
    return cols

def _getNames(cols):
    
    Names = cols.split(',')
    NewName = {}
    for name in Names:
        NewName[name] = name # by default same name as Rob's

    NewName['Longitude'] = 'lon'
    NewName['Latitude'] = 'lat'
        
    NewName['tau_f'] = 'aTau_f'
    NewName['tau_c'] = 'aTau_c'

    NewName['tau_466'] = 'aTau470'
    NewName['tau_553'] = 'aTau550'
    NewName['tau_644'] = 'aTau660'
    NewName['tau_866'] = 'aTau870'

    NewName['mean_AOT0470ea_o'] = 'mTau470'
    NewName['mean_AOT0550ea_o'] = 'mTau550'
    NewName['pval_AOT0550ea_o'] = 'pTau550'
    NewName['mean_AOT0660ea_o'] = 'mTau660'
    NewName['mean_AOT0870ea_o'] = 'mTau870'

    NewName['mean_AOT0470sa_o'] = 'mtau470' # lower case "tau" means "small"
    NewName['mean_AOT0550sa_o'] = 'mtau550'
    NewName['pval_AOT0550sa_o'] = 'ptau550'
    NewName['mean_AOT0660sa_o'] = 'mtau660'
    NewName['mean_AOT0870sa_o'] = 'mtau870'
    NewName['QAavg_o'] = 'qa'

    NewName['mean_mref0470_o']  = 'mRef470'
    NewName['mean_mref0550_o']  = 'mRef550'
    NewName['mean_mref0660_o']  = 'mRef660'
    NewName['mean_mref0870_o']  = 'mRef870'
    NewName['mean_mref1200_o']  = 'mRef1200'
    NewName['mean_mref1600_o']  = 'mRef1600'
    NewName['mean_mref2100_o']  = 'mRef2100'

    NewName['mean_AOT0470corr_l'] = 'mTau470'
    NewName['mean_AOT0550corr_l'] = 'mTau550'
    NewName['pval_AOT0550corr_l'] = 'pTau550'
    NewName['mean_AOT0660corr_l'] = 'mTau660'
    NewName['mean_AOT2100corr_l'] = 'mTau2100'
    NewName['mean_surfre0470_l']  = 'mSre470'
    NewName['mean_surfre0660_l']  = 'mSre660'
    NewName['mean_surfre2100_l']  = 'mSre2100'

    NewName['mean_mref0470_l']  = 'mRef470'
    NewName['mean_mref0550_l']  = 'mRef550'
    NewName['mean_mref0660_l']  = 'mRef660'
    NewName['mean_mref0870_l']  = 'mRef870'
    NewName['mean_mref1200_l']  = 'mRef1200'
    NewName['mean_mref1600_l']  = 'mRef1600'
    NewName['mean_mref2100_l']  = 'mRef2100'

    NewName['mean_AOT0412dpbl_l'] = 'mTau412'
    NewName['mean_AOT0470dpbl_l'] = 'mTau470'
    NewName['mean_AOT0550dpbl_l'] = 'mTau550'
    NewName['mean_AOT0560dpbl_l'] = 'mTau660'
    NewName['mean_SSA0412dpbl_l']  = 'mSsa412'
    NewName['mean_SSA0470dpbl_l']  = 'mSsa470'
    NewName['mean_SSA0660dpbl_l']  = 'mSsa660'
    NewName['mean_surfre0412dpbl_l']  = 'mSre412'
    NewName['mean_surfre0470dpbl_l']  = 'mSre470'
    NewName['mean_surfre0660dpbl_l']  = 'mSre660'

    NewName['QA0470_l'] = 'qa'
    NewName['QAdark_l'] = 'qa'
    NewName['QAdpbl_l'] = 'qa'

    NewName['mean_atype_l'] = 'mAtype'
    NewName['mean_atype_o'] = 'mAtype'

    NewName['mean_cfrac_l'] = 'cloud'
    NewName['mean_cfrac_o'] = 'cloud'
    NewName['ncep_windspd'] = 'speed'

    return (Names, NewName)

#----
class ANET(object):
    """Base class for colocation objects."""

    def __init__ (self,fname,xVars=(),csvVersion=1):

        if csvVersion == 1:
            cols = _getCols1()
        elif csvVersion == 2:
            cols = _getCols2()
        else:
            raise ValueError, 'Invalid csvVersion %d'%csvVersion

        Names, NewName = _getNames(cols)

        Meta = ( 'Longitude',    'Latitude', 'Date', 'Time',
                 'SolarZenith',   'SolarAzimuth', 'SensorZenith', 
                 'SensorAzimuth', 'ScatteringAngle' )

        self.filename = fname
        if 'Aqua' in fname:     self.sat = 'Aqua'
        elif 'Terra' in fname:  self.sat = 'Terra'
        else:                   self.sat = 'Unknown'
            
        Vars = Meta + xVars

        iVars = ()
        formats = ()
        converters = {}
        for name in Vars:
            try:
                i = Names.index(name)
            except:
                raise ValueError, "cannot find <%s> in file"%name
            iVars = iVars + (i,)
            if name=='Date':
                formats = formats + ('S10',)
            elif name=='Time':
                formats = formats + ('S5',)
#            elif name=='QAdark_l':
#                formats = formats + ('S5',)
#            elif name=='QAavg_o':
#                formats = formats + ('S5',)
            else:
                converters[i] = lambda s: float(s or MISSING)
                formats = formats + ('f4',)

#       Read the data
#       -------------
        data = loadtxt(fname, delimiter=',',
                       dtype={'names':Vars,'formats':formats},
                       converters = converters,
                       skiprows=1, usecols=iVars)
        N = len(data)

        self.N = N

#       Save data columns as attributes, with nicer names
#       -------------------------------------------------
        for i in range(len(Vars)):
            if formats[i]=='f4':
                v = ones(N)
                for j in range(N):
                    v[j] = data[j][i]
            else:
                v = []
                for j in range(N):
                    v.append(data[j][i])

            self.__dict__[Vars[i]] = v

            if NewName[Vars[i]] != Vars[i]:
                self.__dict__[NewName[Vars[i]]] = v # alias

        # Add glint angle (TODO: Check whether SensorAzimuth is relative
        # sensor Azimulth)
        # --------------------------------------------------------------
        d2r = pi / 180.
        r2d = 180. / pi
        RelativeAzimuth = abs(self.SolarAzimuth - self.SensorAzimuth - 180.)
        cosGlintAngle = cos(self.SolarZenith*d2r) * cos(self.SensorZenith*d2r) + \
                        sin(self.SolarZenith*d2r) * sin(self.SensorZenith*d2r) * \
                        cos(RelativeAzimuth*d2r)
        
        i = (abs(cosGlintAngle)<=1.0)
        self.GlintAngle = MISSING * ones(self.SolarZenith.size)
        self.GlintAngle[i] = arccos(cosGlintAngle[i])*r2d

#        The code below was used to check the glint angle calculation;
#          using this formula I was able to reproduce the scattering
#          angle on file.
#        cosScatAngle  = -cos(self.SolarZenith*d2r) * cos(self.SensorZenith*d2r) + \
#                        sin(self.SolarZenith*d2r) * sin(self.SensorZenith*d2r) * \
#                        cos(RelativeAzimuth*d2r)
#        i = (abs(cosScatAngle)<=1.0)
#        self.ScatAngle = MISSING * ones(self.SolarZenith.size)
#        self.ScatAngle[i] = arccos(cosScatAngle[i])*r2d

#       Create grads time
#       -----------------
        self.time = []
        for i in range(N):
               self.time.append(gatime(self.Date[i],self.Time[i]))
   
#---
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
                t = self.time[i][:-4] + str(clmYear) # replace year

            ga('set lon %f %f'%(x-1.,x+1.),Quiet=True)
            ga('set lat %f %f'%(y-1.,y+1.),Quiet=True)
            ga('set time %s'%t,Quiet=True)

            try:
                u, levs = ga.interp(expr, lons=(x,),lats=(y,))
                U[i] = u.data
            except:
                ga.flush()

            if U[i] >= 0.0:
                print "%s %s %8.3f %8.3f %6.3f  ...%8.3f%%"\
                      %(self.Date[i],self.Time[i],x,y,U[i],100.*i/float(N))

        self.__dict__[vname] = U

        version = 1
        meta = [ version, self.filename, self.sat, self.surface, vname, expr ]
        savez(outfile,meta=meta,lon=self.lon,lat=self.lat,
              date=self.Date,time=self.Time,var=U)

#---
    def addVarNear(self,ga,outfile,expr='tau',vname=None, clmYear=None):
        """
        Given a grads object having the correct file as default,
        writes out a CSV file with the specified expression/variable. When *clmYear* is
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
                t = self.time[i][:-4] + str(clmYear) # replace year

            ga('set lon '+str(x))
            ga('set lat '+str(y))
            ga('set time %s'%t,Quiet=True)

            U[i] = ga.expr(expr).data # nearest neighbor

            if U[i]>=-1:
                print "%s %s %8.3f %8.3f %6.3f  ...%8.3f%%"\
                      %(self.Date[i],self.Time[i],x,y,U[i],100.*i/float(N))
            else:
                U[i] = MISSING

        self.__dict__[vname] = U

        version = 1
        meta = [ version, self.filename, vname, expr ]
        savez(outfile,meta=meta,lon=self.lon,lat=self.lat,
              date=self.Date,time=self.Time,var=U)

#---
class OCEAN(ANET):
    """Loads a MODIS over Ocean."""

    def __init__(self,fname,csvVersion=1):
        ANET.__init__(self,fname, csvVersion=csvVersion,
                      xVars=('tau_466', 'tau_553', 'tau_644', 'tau_866',
                             'tau_f', 'tau_c',
                             'mean_AOT0470ea_o',
                             'mean_AOT0550ea_o',
                             'pval_AOT0550ea_o', 
                             'mean_AOT0660ea_o',
                             'mean_AOT0870ea_o',
                             'mean_AOT0470sa_o',
                             'mean_AOT0550sa_o',
                             'pval_AOT0550sa_o', 
                             'mean_AOT0660sa_o',
                             'mean_AOT0870sa_o',
                             'mean_mref0470_o',
                             'mean_mref0550_o',
                             'mean_mref0660_o',
                             'mean_mref0870_o',
                             'mean_mref1200_o',
                             'mean_mref1600_o',
                             'mean_mref2100_o',
                             'QAavg_o', 'mean_cfrac_o', 'ncep_windspd'))

        self.surface = 'ocean'
        if self.sat == 'Aqua':
            self.ident = 'mydo'
        elif self.sat == 'Terra':
            self.ident = 'modo'
        else:
            self.ident = 'xxxo'

#       Derive large particle AOD
#       -------------------------
        self.mTAU470 = MISSING * ones(self.mTau550.size)
        self.mTAU550 = MISSING * ones(self.mTau550.size)
        self.mTAU660 = MISSING * ones(self.mTau550.size)
        self.mTAU870 = MISSING * ones(self.mTau550.size)

        self.mTAU470[self.mtau470>=0] = self.mTau470 - self.mtau470 # subtract small from total
        self.mTAU550[self.mtau550>=0] = self.mTau550 - self.mtau550 # subtract small from total
        self.mTAU660[self.mtau660>=0] = self.mTau660 - self.mtau660 # subtract small from total
        self.mTAU870[self.mtau870>=0] = self.mTau870 - self.mtau870 # subtract small from total

#---
    def getCoxMunk(self,filename='/nobackup/NNR/Misc/coxmunk_lut.npz',channel=550):
        """
        Returns ocean albedo.
        """
        
        # Get precomputed albedo LUT
        # --------------------------
        lut = NPZ(filename)
        
        # Trimmed wind speed
        # ------------------
        w10m = self.wind.copy()
        w10m[w10m<0] = 0
        w10m[w10m>50.] = 50.

        j = list(lut.channels).index(channel)

        # Interpolate albedo
        # ------------------
        albedo = zeros(len(w10m))
        albedo[:] = interp(w10m,lut.speed,lut.albedo[:,j])

        self.albedo = albedo

#---

class LAND(ANET):
    """Loads a MODIS over Land."""

    def __init__(self,fname,csvVersion=1):
        if csvVersion==1:
            qa_name = 'QA0470_l'
        else:
            qa_name = 'QAdark_l'
#            qa_name = 'QAdpbl_l'
            
        ANET.__init__(self,fname, csvVersion=csvVersion,
                      xVars=('tau_466', 'tau_553', 'tau_644', 
                             'tau_f', 'tau_c',
                             'mean_AOT0470corr_l',
                             'mean_AOT0550corr_l', 
                             'pval_AOT0550corr_l',
                             'mean_AOT0660corr_l',
                             'mean_AOT2100corr_l',
                             'mean_surfre0470_l',
                             'mean_surfre0660_l',
                             'mean_surfre2100_l',
                             'mean_mref0470_l',
                             'mean_mref0550_l',
                             'mean_mref0660_l',
                             'mean_mref0870_l',
                             'mean_mref1200_l',
                             'mean_mref1600_l',
                             'mean_mref2100_l',
                             'mean_atype_l',
                             'mean_cfrac_l',qa_name))
        self.surface = 'land'
        if self.sat == 'Aqua':
            self.ident = 'mydl'
        elif self.sat == 'Terra':
            self.ident = 'modl'
        else:
            self.ident = 'xxxl'
        
#---

class DEEP(ANET):
    """Loads MODIS Deep Blue."""

    def __init__(self,fname,csvVersion=1):
        if csvVersion==1:
            qa_name = 'QA0470_l'
        else:
            qa_name = 'QAdark_l'
#            qa_name = 'QAdpbl_l'
            
        ANET.__init__(self,fname, csvVersion=csvVersion,
                      xVars=('tau_466', 'tau_553', 'tau_644', 
                             'mean_AOT0412dpbl_l',
                             'mean_AOT0470dpbl_l',
                             'mean_AOT0550dpbl_l', 
                             'mean_AOT0660dpbl_l', 
                             'mean_mref0470_l',
                             'mean_mref0550_l',
                             'mean_mref0660_l',
                             'mean_mref0870_l',
                             'mean_mref1200_l',
                             'mean_mref1600_l',
                             'mean_mref2100_l',
                             'mean_surfre0412dpbl_l',
                             'mean_surfre0470dpbl_l',
                             'mean_surfre0660dpbl_l',
                             'mean_SSA0412dpbl_l',
                             'mean_SSA0470dpbl_l',
                             'mean_SSA0660dpbl_l',
                             'mean_cfrac_l',qa_name))
        self.surface = 'land'
        if self.sat == 'Aqua':
            self.ident = 'mydl'
        elif self.sat == 'Terra':
            self.ident = 'modl'
        else:
            self.ident = 'xxxl'
        
#....................................................................

if __name__ == "__main__":

    from grads import GrADS
    ga = GrADS(Echo=False,Window=False)

    ga('xdfopen merra_slv-hourly.ddf')

    mxdo = OCEAN('SUPER2_combo.Aqua.csv')
    mxdo.addVar(ga,'mydo2_merra_wind.npz',expr='mag(u10m,v10m)',vname='wind')

def _savethis():

    ga('open albedo_clim.ctl')
    mxdl = LAND('SUPER_land.Aqua.csv',csvVersion=1)
    mxdl.addVar(ga,mxdl.ident+'_albedo.npz',expr='albedo',clmYear=2000)

    misr = ANET('SUPER2_combo.Terra.csv',csvVersion=2,
                xVars=('tau_466', 'tau_553', 'tau_644', 'tau_866'))

    ga('xdfopen misr.ddf')
    ga('set lev 550')
    misr.addVarNear(ga,'misr.npz',expr='tauext',vname='misrTau550')

    ga('open albedo_clim.ctl')
    mxdl = LAND('SUPER_land.Terra.csv',csvVersion=1)
    mxdl.addVar(ga,mxdl.ident+'_albedo.npz',expr='albedo',clmYear=2000)

    ga('xdfopen merra_slv-hourly.ddf')

    mxdo = OCEAN('SUPER2_ocean.Aqua.csv')
    mxdo.addVar(ga,'mydo2_merra_wind.npz',expr='mag(u10m,v10m)',vname='wind')
    mxdo = OCEAN('SUPER2_ocean.Terra.csv')
    mxdo.addVar(ga,'modo2_merra_wind.npz',expr='mag(u10m,v10m)',vname='wind')

    modo = OCEAN('SUPER_ocean.Terra.csv')
    modo.addVar(ga,'modo_merra_ustar.npz',expr='ustar')
#    modo.addVar(ga,'modo_merra_speed.npz',expr='speed')
    del modo
    
    mydo = OCEAN('SUPER_ocean.Aqua.csv')
    mydo.addVar(ga,'mydo_merra_ustar.npz',expr='ustar')

#    mydo.addVar(ga,'mydo_merra_speed.npz',expr='speed')
        
#    mydo = OCEAN('SUPER_ocean.Aqua.csv')
#    modl = LAND('SUPER_land.Terra.csv')
#    mydl = LAND('SUPER_land.Aqua.csv')



