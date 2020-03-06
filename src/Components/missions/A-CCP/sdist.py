#!/usr/bin/env python
"""
class to do size distribution calculations
"""
import os
from netCDF4 import Dataset
import numpy as np
from   MAPL  import config
from scipy import interpolate
from   py_leo_vlidort.copyvar  import _copyVar

class SDIST(object):
    #---
    def sizeDistribution(self):
        """
        Get aerosol size distribution from GEOS aerosol mixing ratios
        Specifically written for optics tables generated for A-CCP simulations
        made by Oksu.  Will not work with Pete's old tables
        """

        # Create master bins for all levels
        # ---------------------------------
        rmin      = 0.001e-6   #meters 
        rmax      = self.getRMAX()
        nbins     = 100
        RE        = np.linspace(rmin,1.1*rmax,101)
        DR        = RE[1:] - RE[0:-1]
        R         = RE[0:-1] + 0.5*DR
        RLOW      = RE[0:-1]
        RUP       = RE[1:]
        
        self.R = R
        self.DR = DR
        self.RLOW = RLOW
        self.RUP  = RUP
        self.RMIN = rmin
        self.RMAX = rmax
        self.nbins = nbins

        # Lognormal Species: BC, OC, and SU
        # ---------------------------------
        logspc = ['BCPHOBIC','BCPHILIC','OCPHOBIC','OCPHILIC','SU']
        for spc in logspc:
            if spc == 'SU':
                spc_ = 'SO4'
            else:
                spc_ = spc
            if spc_ in self.AERNAMES:
                self.logNormalDistribution(spc)

        # Dust
        # ----------
        if 'DU001' in self.AERNAMES:
           self.dustDistribution()

        # Sea Salt
        # ------------
        if 'SS001' in self.AERNAMES:
           self.seasaltDistribution()

        # Mixture Distribution
        # -------------
        #self.TOTdist = self.BCPHILICdist + self.BCPHOBICdist + self.OCPHILICdist + self.OCPHOBICdist + self.SUdist + self.DUdist + self.SSdist
        for spc in self.AERdistNAMES:
            spcdist = spc+'dist'     
            
            if hasattr(self,'TOTdist'):
                self.TOTdist += self.__dict__[spcdist]
            else:
                self.TOTdist = self.__dict__[spcdist]

        # Column size distribution
        nobs,nlev,nR = self.TOTdist.shape
        for spc in self.AERdistNAMES+['TOT']:
            spcdist = spc + 'dist'
            self.__dict__['col'+spcdist] = np.zeros([nobs,nR])
            for t in range(nobs):
                dz = self.ze[:-1,t] - self.ze[1:,t]
                for r in range(nR):
                    self.__dict__['col'+spcdist][t,r] = np.sum(self.__dict__[spcdist][t,:,r]*dz)

#        # convert to dV/dlnR [microns^3/microns^2] 
#        spc = 'colTOTdist'
#        temp = np.zeros(self.__dict__[spc].shape)  
#        for t in range(nobs):
#            temp[t,:] = self.__dict__[spc][t,:]*self.R*1e6

#        self.__dict__[spc] = temp


        # effective radius for total only
        self.TOTreff = np.zeros([nobs,nlev])
        self.colTOTreff = np.zeros([nobs])
        R = self.R
        DR = self.DR
        R3 = R**3
        R2 = R**2
        for t in range(nobs):
            dndr = self.colTOTdist[t,:]/R3
            self.colTOTreff[t] = np.sum(R3*dndr*DR)/np.sum(R2*dndr*DR)
            for k in range(nlev):
                dndr = self.TOTdist[t,k,:]/R3
                self.TOTreff[t,k] = np.sum(R3*dndr*DR)/np.sum(R2*dndr*DR)


    def logNormalDistribution(self,spc):
        # master bins
        R    = self.R
        DR   = self.DR
        RLOW = self.RLOW
        RUP  = self.RUP

        # Read optics table
        cf = config.Config(self.rcFile)
        optable = cf('filename_optical_properties_{}'.format(spc[0:2]))

        if spc == 'SU':
            irad = 0
        elif 'PHOBIC' in spc:
            irad = 0
        elif 'PHILIC' in spc:
            irad = 1

        # Read in size distribution information
        nc = Dataset(optable)
        rhTable    = nc.variables['rh'][:]
        rminTable  = nc.variables['rMin'][irad,:]
        rmaxTable  = nc.variables['rMax'][irad,:]
        rmodeTable = nc.variables['rMode'][irad,:]
        sigma      = nc.variables['sigma'][irad]
        rhopTable  = nc.variables['rhop'][irad,:]
        gfTable    = nc.variables['growth_factor'][irad,:]
        
        nc.close()

        rminTable = interpolate.interp1d(rhTable, rminTable)
        rmaxTable = interpolate.interp1d(rhTable, rmaxTable)
        rmodeTable = interpolate.interp1d(rhTable, rmodeTable)
        rhopTable  = interpolate.interp1d(rhTable, rhopTable)
        gfTable    = interpolate.interp1d(rhTable, gfTable)
      
        # rh dims are [ntyme,nlev]
        # mr (mixing ratio) dims are [ntyme,nlev]
        # air (air density) dims are [ntyme,nlev]
        if spc == 'SU':
            spc_ = 'SO4'
        else:
            spc_ = spc
        rh  = self.RH.copy()
        mr  = self.__dict__[spc_]
        air = self.AIRDENS

        nobs,nlev = air.shape

        # set up empty column size distribution array
        SPCdist = np.empty([nobs,nlev,len(R)])

        # 0 <= rh <= 0.99
        # this is what is done in Chem_MieMod.F90
        rh[rh < 0] = 0
        rh[rh > 0.99] = 0.99   

        # loop through orbit
        for t in range(nobs):
            rmode = rmodeTable(rh[t,:])
            rmin  = rminTable(rh[t,:])
            rmax  = rmaxTable(rh[t,:])
            rhop0 = rhopTable(0)
            gf    = gfTable(rh[t,:])

            # loop through layers
            # get size distribution for each layer
            for k in range(nlev):
                # get the aerosol number distribution
                rNum = rmode[k]
                lsigma = np.log(sigma)
                C      = np.sqrt(2.*np.pi)
                dndr = (1./(R*lsigma*C))*np.exp(-(np.log(R/rNum)**2.)/(2.*lsigma**2.)) 

                # Truncate distribution according to rlow and rup
                ii = RUP <= rmin[k]
                dndr[ii] = 0

                ii = RLOW >= rmax[k]
                dndr[ii] = 0

                # deal with lowest bin
                # number concentration is scaled to the
                # fraction of the bin that is covered by
                # the aerosol distribution
                if any(RLOW < rmin[k]):
                    bini = np.arange(len(R))
                    i = bini[RLOW < rmin[k]][-1]
                    drtilda = RUP[i] - rmin[k]
                    dndr[i] = dndr[i]*drtilda/DR[i]

                #deal with the highest bin
                # number concentration is scaled to the
                # fraction of the bin that is covered by
                # the aerosol distribution
                if any(RLOW < rmax[k]):
                    bini = np.arange(len(R))
                    i = bini[RLOW < rmax[k]][-1]
                    drtilda = rmax[k] - RLOW[i]
                    dndr[i] = dndr[i]*drtilda/DR[i]

                # Now get the volume distribution
                # dvdr
                dvdr = 4./3.*np.pi*R**3.*dndr
                
                # Get aerosol DRY! volume concentration
                M0 = mr[t,k]*air[t,k]
                V0 = M0/rhop0

                # Get the Wet volume concentration
                Vwet = V0*gf[k]**3

                # normalize dvdr so the integral is equal to the 
                # wet volume concentration    
                dvdr = dvdr*Vwet/np.sum(dvdr*DR)

                # save this layer aerosol distribution to the column
                SPCdist[t,k,:] = dvdr

        self.__dict__[spc+'dist'] = SPCdist

    def seasaltDistribution(self):

        #constants for adjusting size bins for RH
        c1 = 0.7674
        c2 = 3.079
        c3 = 2.573e-11
        c4 = -1.424        

        # master bins
        R    = self.R
        DR   = self.DR
        RLOW = self.RLOW
        RUP  = self.RUP

        # Read optics table
        cf = config.Config(self.rcFile)
        optable = cf('filename_optical_properties_SS')

        # Read in growth factors and rh tables
        nc = Dataset(optable)
        rhTable = nc.variables['rh'][:]
        gfTable = nc.variables['growth_factor'][:]

        # Density of dry particles [kg m-3]
        # dimensions [radisu,rh]
        rhop = nc.variables['rhop'][:]
        rhop0 = rhop[:,0]

        # Major bins
        rMaxMaj = nc.variables['rUp'][:]
        rMinMaj = nc.variables['rLow'][:]

        nc.close()

        # number of major bins
        nbinMaj = len(rMaxMaj)

        # rh (relative humidity) dims are [ntyme,nlev]
        # mr (mixing ratio) dims are [ntyme,nlev]
        # air (air density) dims are [ntyme,nlev]
        rh = self.RH.copy()
        ss001,ss002,ss003,ss004,ss005 = self.SS001,self.SS002,self.SS003,self.SS004,self.SS005
        air = self.AIRDENS

        # set up empty column size distribution array
        nobs,nlev = air.shape
        SPCdist = np.empty([nobs,nlev,len(R)])

        # put all sea-salt mixing ratios in one array for convenience
        SS = np.zeros([nobs,nlev,5])
        SS[:,:,0] = ss001
        SS[:,:,1] = ss002
        SS[:,:,2] = ss003
        SS[:,:,3] = ss004
        SS[:,:,4] = ss005


        # 0 <= rh <= 0.95
        # this is what is done in Chem_MieMod.F90
        rh[rh < 0] = 0
        rh[rh > 0.95] = 0.95   

        # loop throuhg time steps
        for t in range(nobs):    
            # loop through major bins
            for iBin in range(nbinMaj):
                # get the radii of the bin
                rmin = rMinMaj[iBin]
                rmax = rMaxMaj[iBin]
                rMinCM = rmin*100.
                rMaxCM = rmax*100.                    

                # get growth factors table for this bin
                fTable = interpolate.interp1d(rhTable, gfTable[iBin,:])
                gf = fTable(rh[t,:])  

                # loop through layers 
                # get size distribution for each layer
                for k in range(nlev):

                    #adjust bin edges for humidified particles
                    rhUse = rh[t,k]
                    rMinUse = (c1*rMinCM**c2 /(c3*rMinCM**c4 - np.log10(rhUse))+rMinCM**3.)**(1./3.)/100.                                      
                    rMaxUse = (c1*rMaxCM**c2 /(c3*rMaxCM**c4 - np.log10(rhUse))+rMaxCM**3.)**(1./3.)/100.  

                    # Determine the dNdr of the particle size distribution using the
                    # Gong 2003 particle sub-bin distribution
                    rrat = rmin/rMinUse
                    r80Rat = 1.65*rrat      # ratio of the r80 radius to the wet radius
                    r80  = R*r80Rat * 1.e6  # radius in r80 space in um
                    dr80 = DR*r80Rat * 1.e6

                    aFac = 4.7*(1.+30.*r80)**(-0.017*r80**(-1.44))
                    bFac = (0.433-np.log10(r80))/0.433
                    dndr80 = 1.373*r80**(-aFac)*(1.+0.057*r80**3.45)*10.**(1.607*np.exp(-bFac**2.))
                    dndr = dndr80 * r80Rat

                    # Truncate distribution according to rlow and rup
                    ii = RUP <= rMinUse
                    dndr[ii] = 0

                    ii = RLOW >= rMaxUse
                    dndr[ii] = 0

                    # deal with lowest bin
                    # number concentration is scaled to the
                    # fraction of the bin that is covered by
                    # the aerosol distribution
                    if any(RLOW < rMinUse):
                        bini = np.arange(len(R))
                        i = bini[RLOW < rMinUse][-1]
                        drtilda = RUP[i] - rMinUse
                        dndr[i] = dndr[i]*drtilda/DR[i]

                    #deal with the highest bin
                    # number concentration is scaled to the
                    # fraction of the bin that is covered by
                    # the aerosol distribution
                    if any(RLOW < rMaxUse):
                        bini = np.arange(len(R))
                        i = bini[RLOW < rMaxUse][-1]
                        drtilda = rMaxUse - RLOW[i]
                        dndr[i] = dndr[i]*drtilda/DR[i]

                    # Now get the volume distribution
                    # dvdr
                    dvdr = 4./3.*np.pi*R**3.*dndr

                    # Get aerosol DRY! volume concentration
                    mr = SS[t,k,iBin]
                    M0 = mr*air[t,k]
                    V0 = M0/rhop0[iBin]

                    # Get the Wet volume
                    Vwet = V0*gf[k]**3

                    # normalize dvdr so the integral is equal to the 
                    # wet volume
                    dvdr = dvdr*Vwet/np.sum(dvdr*DR)

                    # add this aerosol distribution to the master
                    SPCdist[t,k,:] = SPCdist[t,k,:] + dvdr

        self.__dict__['SSdist'] = SPCdist

    def dustDistribution(self):
        # master bins convert to microns
        R    = self.R*1e6
        DR   = self.DR*1e6
        RLOW = self.RLOW*1e6
        RUP  = self.RUP*1e6

        # Read optics table
        cf = config.Config(self.rcFile)
        optable = cf('filename_optical_properties_DU')

        # Open Table
        nc = Dataset(optable)

        # lognormal volume median radius in microns
        rv = nc.variables['rv'][:]

        # lognormal sigma
        sigma = nc.variables['sigma'][:]


        # Density of dry particles [kg m-3]
        # hard coded for dust tables developed from GRASP
        rhop0 = np.array([ 2500.,  2650.,  2650.,  2650.,  2650.])

        # Major bins
        # hard coded for dust tables developed from GRASP
        # microns
        RMAX = 20
        RMIN = 0.08
        rMaxMaj = np.array([RMAX,RMAX,RMAX,RMAX,RMAX])
        rMinMaj = np.array([RMIN,RMIN,RMIN,RMIN,RMIN])

        nc.close()

        # number of major bins
        nbinMaj = len(rMaxMaj)

        # mr (mixing ratio) dims are [ntyme,nlev]
        # air (air density) dims are [ntyme,nlev]
        du001,du002,du003,du004,du005 = self.DU001,self.DU002,self.DU003,self.DU004,self.DU005
        air = self.AIRDENS

        # set up empty column size distribution array
        nobs,nlev = air.shape
        SPCdist = np.empty([nobs,nlev,len(R)])

        # put all dust mixing ratios in one array for convenience
        DU = np.zeros([nobs,nlev,5])
        DU[:,:,0] = du001
        DU[:,:,1] = du002
        DU[:,:,2] = du003
        DU[:,:,3] = du004
        DU[:,:,4] = du005

        # loop throuhg time steps
        for t in range(nobs):    
            # loop through major bins
            for iBin in range(nbinMaj):
                # get the radii of the bin
                rmin = rMinMaj[iBin]
                rmax = rMaxMaj[iBin]
                rvUse = rv[iBin]
                sigmaUse = sigma[iBin]                

                # loop through layers 
                # get size distribution for each layer
                for k in range(nlev):

                    # get the aerosol volume distribution
                    C      = np.sqrt(2.*np.pi)
                    dvdr = (1./(R*sigmaUse*C))*np.exp(-((np.log(R)-rvUse)**2.)/(2.*sigmaUse**2.)) 

                    # Truncate distribution according to rmin and rmax
                    ii = RUP <= rmin
                    dvdr[ii] = 0

                    ii = RLOW >= rmax
                    dvdr[ii] = 0

                    # deal with lowest bin
                    # number concentration is scaled to the
                    # fraction of the bin that is covered by
                    # the aerosol distribution
                    if any(RLOW < rmin):
                        bini = np.arange(len(R))
                        i = bini[RLOW < rmin][-1]
                        drtilda = RUP[i] - rmin
                        dvdr[i] = dvdr[i]*drtilda/DR[i]

                    #deal with the highest bin
                    # number concentration is scaled to the
                    # fraction of the bin that is covered by
                    # the aerosol distribution
                    if any(RLOW < rmax):
                        bini = np.arange(len(R))
                        i = bini[RLOW < rmax][-1]
                        drtilda = rmax - RLOW[i]
                        dvdr[i] = dvdr[i]*drtilda/DR[i]

                    # Get aerosol DRY! volume concentration
                    mr = DU[t,k,iBin]
                    M0 = mr*air[t,k]
                    V0 = M0/rhop0[iBin]

                    # normalize dvdr so the integral is equal to the 
                    # volume concentration
                    dvdr = dvdr*V0/np.sum(dvdr*DR)

                    # convert dr from microns to meters
                    dvdr = dvdr*1e6

                    # add this aerosol distribution to the master
                    SPCdist[t,k,:] = SPCdist[t,k,:] + dvdr

        self.__dict__['DUdist'] = SPCdist


    def getRMAX(self):
        cf = config.Config(self.rcFile)

        #lognormals
        spclist = 'BC','OC','SU'
        RMAX    = 0
        for spc in spclist:
            # Read optics table            
            optable = cf('filename_optical_properties_{}'.format(spc))
            nc      = Dataset(optable)

            if spc == 'SU':
                rmax    = nc.variables['rMax'][0,:]
            else:
                rmax    = nc.variables['rMax'][1,:]

            nc.close()
            rmax  = rmax.max()
            if rmax > RMAX: RMAX = rmax

        # Sea-salt
        c1 = 0.7674
        c2 = 3.079
        c3 = 2.573e-11
        c4 = -1.424        
        rhUse = 0.95
        # read optics table
        optable = cf('filename_optical_properties_SS')
        nc      = Dataset(optable)
        rUp     = nc.variables['rUp'][:]
        nc.close()
        rMaxUse = rUp[-1]
        rMaxCM  = rMaxUse*100.
        rmax    = (c1*rMaxCM**c2 /(c3*rMaxCM**c4 - np.log10(rhUse))+rMaxCM**3.)**(1./3.)/100.
        if rmax > RMAX: RMAX = rmax

        # Dust
        # hard coded 20 microns for tables generated from GRASP
        rmax    = 20
        rmax    = rmax*1e-6   #meters
        if rmax > RMAX: RMAX = rmax

        return RMAX

    def writenc(self):
        """
        write a netcdf File of size distributions
        """
        if not os.path.exists(os.path.dirname(self.outFile)):
            os.makedirs(os.path.dirname(self.outFile))

        # Open NC file
        # ------------
        nc = Dataset(self.outFile,'w')	

        # Set global attributes
        # ---------------------
        nc.title = 'GEOS Size Distribution'
        nc.institution = 'NASA/Goddard Space Flight Center'
        nc.source = 'Global Model and Assimilation Office'
        nc.history = ''
        nc.references = 'n/a'
        nc.contact = 'Patricia Castellanos <patricia.castellanos@nasa.gov>'
        nc.Conventions = 'CF'


        # Open extFile for reading
        nctrj = Dataset(self.inFile.replace('%col','aer_Nv'))

        # Create dimensions
        # -----------------
        ntime,nlev = self.AIRDENS.shape
        nt = nc.createDimension('time',ntime)
        ls = nc.createDimension('ls',19)
        x  = nc.createDimension('x',1)
        y  = nc.createDimension('y',1)
        r  = nc.createDimension('r',len(self.R))
        l  = nc.createDimension('lev',nlev)

        _copyVar(nctrj,nc,u'trjLon',dtype='f4',zlib=False,verbose=self.verbose)
        _copyVar(nctrj,nc,u'trjLat',dtype='f4',zlib=False,verbose=self.verbose)
        _copyVar(nctrj,nc,u'time', dtype='i4',zlib=False,verbose=self.verbose)
        _copyVar(nctrj,nc,u'isotime', dtype='S1',zlib=False,verbose=self.verbose)
        _copyVar(nctrj,nc,u'x',dtype='f4',zlib=False,verbose=self.verbose)
        _copyVar(nctrj,nc,u'y',dtype='f4',zlib=False,verbose=self.verbose)

        nctrj.close()

        # Create Variables
        # ------------------
        dim = ('r',)
        this = nc.createVariable('radius','f4',dim,zlib=True)
        this.units = 'm'
        this[:] = self.R

        this = nc.createVariable('dradius','f4',dim,zlib=True)
        this.long_name = 'radius bin width'
        this.units = 'm'
        this[:] = self.DR

        this = nc.createVariable('rlow','f4',dim,zlib=True)
        this.long_name = 'radius bin lower edge'
        this.units = 'm'
        this[:] = self.RLOW

        this = nc.createVariable('rup','f4',dim,zlib=True)
        this.long_name = 'radius bin upper edge'
        this.units = 'm'
        this[:] = self.RUP

        dim = ('time','lev','r',)
        for spc in self.AERdistNAMES+['TOT']:
            spcdist = spc+'dist'
            print spcdist
            this = nc.createVariable(spcdist,'f4',dim,zlib=True)
            this.long_name = spc + ' size distribution (dV/dr)'
            this.units = 'm^3/m'
            this[:] = self.__dict__[spcdist]

        dim = ('time','r',)
        for spc in self.AERdistNAMES+['TOT']:
            spcdist = 'col'+spc+'dist'
            this = nc.createVariable(spcdist,'f4',dim,zlib=True)
            this.long_name = spc + 'column size distribution (dV/dr)'
            this.units = 'm^3/m'
            this[:] = self.__dict__[spcdist]


        dim = ('time','lev',)
        this = nc.createVariable('reff_profile','f4',dim,zlib=True)
        this.long_name = 'profile effective radius'
        this.units = 'm'
        this[:] = self.TOTreff

        dim = ('time',)
        this = nc.createVariable('reff_column','f4',dim,zlib=True)
        this.long_name = 'column effective radius'
        this.units = 'm'
        this[:] = self.colTOTreff

        nc.close()

