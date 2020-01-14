#!/usr/bin/env python

"""
    Gets geometry for the polarimeter swath

    Patricia Castellanos, Dec, 2019

"""

import os
import argparse
from   datetime        import datetime, timedelta
from   dateutil.parser import parse         as isoparser
from   MAPL            import Config
import numpy  as np
from netCDF4 import Dataset
from  py_leo_vlidort import LidarAngles_ 
from mpl_toolkits.basemap import Basemap
from   py_leo_vlidort.copyvar  import _copyVar
from multiprocessing import Pool

# Generic Lists of Varnames and Units
SDS_AER    = ['LONGITUDE','LATITUDE','isotime']


ncALIAS = {'LONGITUDE': 'trjLon',
           'LATITUDE': 'trjLat'}

#radius of the earth km
Re = 6378.0
hgtssRef = 500.


# --
def distance(lon0,lat0,lon,lat,Re):
    # get index of closest lon and lat to lon0, lat0

    # use haversine formula for distance
    phi0, phi, lam0, lam = map(np.radians,[lat0,lat,lon0,lon])
    dphi = phi - phi0
    dlam = lam - lam0

    a = np.sin(dphi*0.5)**2 + np.cos(phi0)*np.cos(phi)*np.sin(dlam*0.5)**2
    if np.any(a > 1.0):
        a[a> 1.0] = 1.0
    c = 2.0*np.arcsin(np.sqrt(a))
    d = Re*c

    return d

def get_imin(args):
    lon,lat,lon_granule,lat_granule,Re,Istyme = args
    #lon, lat have dims ntyme
    # lon_granule have dims ntymetotal, ncross 
    ntyme = len(lon)
    ntymeTotal, ncross = lon_granule.shape

    IMIN = np.zeros([ntyme,2]).astype(int)
    for ityme in range(ntyme):
        lon0 = lon[ityme]
        lat0 = lat[ityme]
        # pick a window to search around to speed things up
        styme = np.max([ityme+Istyme-200,0])
        etyme = np.min([ityme+Istyme+200,ntymeTotal])
        
        d = distance(lon0,lat0,lon_granule[styme:etyme,:],lat_granule[styme:etyme,:],Re)
        imin = np.argmin(d)

        itymemin, icrossmin = np.unravel_index(imin,d.shape)
        itymemin = itymemin + styme

        IMIN[ityme,0] = itymemin
        IMIN[ityme,1] = icrossmin

    return IMIN

def SolarAngles(args):
    Tyme,lat_granule,lon_granule = args
    ntyme,ncross = lon_granule.shape
    angles = np.zeros([ntyme,ncross,2])
    
    for ityme,tyme in enumerate(Tyme):
        for icross in range(ncross):
            CLAT = lat_granule[ityme,icross]
            CLON = lon_granule[ityme,icross]

            solar_angles = LidarAngles_.solarangles(tyme.year,tyme.month,tyme.day,
                                                tyme.hour,tyme.minute,tyme.second,
                                                CLAT,CLON,
                                                0.0)
            angles[ityme,icross,0] = solar_angles[1][0]  #sza
            angles[ityme,icross,1] = solar_angles[0][0]  #saa
    return angles

class SWATH(object):
    """
    get polarimeter swath
    """
    def __init__(self,starttyme,inFilem1,inFile,inFilep1,outFile,hgtss,along_track_deg,
                 verbose=True,
                 cross_track_km=None,cross_track_dkm=None):
        self.starttyme = starttyme
        self.inFilem1 = inFilem1
        self.inFilep1 = inFilep1
        self.inFile   = inFile
        self.outFile  = outFile
        self.hgtss    = hgtss
        self.SDS_AER  = SDS_AER
        self.verbose  = verbose
        self.cross_track_km = cross_track_km
        self.cross_track_dkm = cross_track_dkm
        self.Re       = Re
        self.hgtssRef = hgtssRef

        # Parse along_track_deg variable
        if ',' in along_track_deg:
            # list
            self.vza_along = np.array(along_track_deg.split(',')).astype(float)
            self.nalong    = len(self.vza_along)
        elif ':' in along_track_deg:
            vzamin,vzamax,nalong = np.array(along_track_deg.split(':')).astype(float)

            self.vza_along = np.linspace(vzamin,vzamax,nalong)
            self.nalong    = int(nalong)

        # initialize empty lists
        for sds in self.SDS_AER:
            self.__dict__[sds] = []

        # Read in model data
        self.readSampledGEOS()
        self.ntyme = len(self.LONGITUDE[0])
        if self.inFilem1 is None:
            self.Istyme = 0
        else:
            self.Istyme = self.ntyme

        # Make lists into arrays
        for sds in self.SDS_AER:
            self.__dict__[sds] = np.concatenate(self.__dict__[sds])

        # convert isotime to datetime
        self.tyme = []
        for isotime in self.isotime:
            self.tyme.append(isoparser(''.join(isotime)))

        self.tyme = np.array(self.tyme)
        self.ntymeTotal = len(self.tyme)

        # Calculate Nadir Scene Geometry
        self.hgtss = hgtss
        self.nadirAngles()

        # Get cross track swath deg full angle
        if self.cross_track_km is not None:
            self.get_vza_cross()

        # Calculate Calculate Full Pixel Locations and view angles
        self.granulePixels()
        self.granuleSolarAngles()
        istart = self.Istyme 
        iend   = self.Istyme + self.ntyme
        self.scatAngle_granule = self.getScatAngle(self.sza_granule,self.saa_granule,self.vza_granule[istart:iend,:,:],self.vaa_granule[istart:iend,:,:])


        # figure out when satellite pixels overlap with each satellite sub-point
        self.subpointView()

        # Put together arrays of angles for when satellite points at the satellite sub-point
        # and write to file
        self.subpointViewAngles()
        self.writenc()

    def writenc(self):
        """
        write a netcdf File of nadir view angles
        """
        if not os.path.exists(os.path.dirname(self.outFile)):
            os.makedirs(os.path.dirname(self.outFile))

        # Open NC file
        # ------------            
        nc = Dataset(self.outFile,'w')

        # Set global attributes
        # ---------------------
        nc.title = 'Multiangle viewing geometry for satellite nadir pixel'
        nc.institution = 'NASA/Goddard Space Flight Center'
        nc.source = 'Global Model and Assimilation Office'
        nc.history = ''
        nc.references = 'n/a'
        nc.contact = 'Patricia Castellanos <patricia.castellanos@nasa.gov>'
        nc.Conventions = 'CF'
        nc.satellite_alt = self.hgtss
        nc.cross_track_full_deg = self.cross_track_deg

        # Open extFile for reading
        nctrj = Dataset(inFile.replace('%col','aer_Nv'))
        ntime = len(nctrj.dimensions['time'])
     
        # Create dimensions
        # -----------------
        nt = nc.createDimension('time',ntime)
        ls = nc.createDimension('ls',19)
        x  = nc.createDimension('x',1)
        y  = nc.createDimension('y',1)
        nalong = nc.createDimension('nalong',self.nalong)
        ncross = nc.createDimension('ncross',self.ncross)

        _copyVar(nctrj,nc,u'trjLon',dtype='f4',zlib=False,verbose=self.verbose)
        _copyVar(nctrj,nc,u'trjLat',dtype='f4',zlib=False,verbose=self.verbose)
        _copyVar(nctrj,nc,u'time', dtype='i4',zlib=False,verbose=self.verbose)
        _copyVar(nctrj,nc,u'isotime', dtype='S1',zlib=False,verbose=self.verbose)
        _copyVar(nctrj,nc,u'x',dtype='f4',zlib=False,verbose=self.verbose)
        _copyVar(nctrj,nc,u'y',dtype='f4',zlib=False,verbose=self.verbose)   

        nctrj.close()

        # Create Variables
        # ------------------
        dim = ('time','nalong',)
        this = nc.createVariable('time_ss','i4',dim)
        this[:] = self.time
        this.long_name = 'time when polarimeter views satellite subpoint'
        this.units = 'seconds since ' + self.starttyme.isoformat()

        this = nc.createVariable('sza_ss','f4',dim)
        this[:] = self.sza
        this.long_name = 'SZA when polarimeter views satellite subpoint'
        this.units = 'degrees, >90 is below horizon'

        this = nc.createVariable('saa_ss','f4',dim)
        this[:] = self.saa
        this.long_name = 'SAA when polarimeter views satellite subpoint'
        this.units = 'degrees, clockwise from north 0-360'

        this = nc.createVariable('vza_ss','f4',dim)
        this[:] = self.vza
        this.long_name = 'VZA when polarimeter views satellite subpoint'
        this.units = 'degrees'

        this = nc.createVariable('vaa_ss','f4',dim)
        this[:] = self.vaa
        this.long_name = 'VAA when polarimeter views satellite subpoint'
        this.units = 'degrees, clockwise from north 0-360'

        this = nc.createVariable('scatAngle_ss','f4',dim)
        this[:] = self.scatAngle
        this.long_name = 'scattering angle when polarimeter view satellite subpoint'
        this.units = 'degrees, 0 = forward scattering'

        dim = ('time','nalong','ncross',)
        istart = self.Istyme 
        iend   = self.Istyme + self.ntyme

        this = nc.createVariable('longitude','f4',dim)
        this[:] = self.lon_granule[istart:iend,:,:]
        this.long_name = 'longitude of granule pixel centers'
        this.units  = 'degrees_east'

        this = nc.createVariable('latitude','f4',dim)
        this[:] = self.lat_granule[istart:iend,:,:]
        this.long_name = 'latitude of granule pixel centers'
        this.units = 'degrees_north'

        this = nc.createVariable('vza','f4',dim)
        this[:] = self.vza_granule[istart:iend,:,:]
        this.long_name = 'view zenith angle'
        this.units = 'degrees'

        this = nc.createVariable('vaa','f4',dim)
        this[:] = self.vaa_granule[istart:iend,:,:]
        this.long_name = 'view azimuth angle'
        this.units = 'degrees, clockwise from north 0-360'

        this = nc.createVariable('sza','f4',dim)
        this[:] = self.sza_granule
        this.long_name = 'solar zenith angle'
        this.units = 'degrees'

        this = nc.createVariable('saa','f4',dim)
        this[:] = self.saa_granule
        this.long_name = 'solar azimuth angle'
        this.units = 'degrees, clockwise from north 0-360'

        this = nc.createVariable('scatAngle','f4',dim)
        this[:] = self.scatAngle_granule
        this.long_name = 'scattering angle'
        this.units = 'degrees, 0 = forward scattering'        

        nc.close()



    def getScatAngle(self,sza,SAA,vza,VAA):
        """
        The angle between the sun, the pixel and sensor.  A value of
        180.0 would imply pure backscatter.  A value of zero would
        imply pure forward scatter (not seen in satellite remote sensing).
        Values less than 90 imply forward scattering.
        Range is 0 to 180.0 
        """       
        # do this so you don't change original values
        saa = SAA.copy()
        vaa = VAA.copy()
        if np.any(vaa > 180):
            I = vaa > 180
            vaa[I] = vaa[I] - 360.0
        if np.any(saa > 180):
            I = saa > 180
            saa[I] = saa[I] - 360.0

        raa = np.abs(saa-vaa)
        if np.any(raa > 180.0):
            I = raa > 180.0
            raa[I] = 360.0 - raa[I]
        raa = 180.0 - raa

        szar = np.radians(sza)
        vzar = np.radians(vza)
        raar = np.radians(raa)
        scattering_angle = -1.0*np.cos(szar)*np.cos(vzar) - np.sin(szar)*np.sin(vzar)*np.cos(raar)
        if np.any(scattering_angle) > 1.0:
            I = scattering_angle > 1.0
            scattering_angle[I] = 1.0

        if np.any(scattering_angle) < -1.0:
            I = scattering_angle < -1.0
            scattering_angle[I] = -1.0

        scattering_angle = np.arccos(scattering_angle)
        scattering_angle = np.degrees(scattering_angle)

        return scattering_angle

    # --
    def subpointViewAngles(self):
        sza = np.zeros([self.ntyme,self.nalong])
        saa = np.zeros([self.ntyme,self.nalong])
        vza = np.zeros([self.ntyme,self.nalong])
        vaa = np.zeros([self.ntyme,self.nalong])
        scatAngle = np.zeros([self.ntyme,self.nalong])
        time = np.zeros([self.ntyme,self.nalong]).astype(int)

        for ialong in range(self.nalong):
            sza[:,ialong] = self.SZA[self.Itymeview[:,ialong]]
            saa[:,ialong] = self.SAA[self.Itymeview[:,ialong]]
            dtyme = self.tyme[self.Itymeview[:,ialong]] - self.starttyme
            for ityme in range(self.ntyme):
                vza[ityme,ialong] = self.vza_granule[:,ialong,:][self.Itymeview[ityme,ialong],self.Icrossview[ityme,ialong]]
                vaa[ityme,ialong] = self.vaa_granule[:,ialong,:][self.Itymeview[ityme,ialong],self.Icrossview[ityme,ialong]]
                time[ityme,ialong] = dtyme[ityme].total_seconds()

        self.sza = np.ma.array(sza)
        self.saa = np.ma.array(saa)
        self.vza = np.ma.array(vza)
        self.vaa = np.ma.array(vaa)
        self.time = np.ma.array(time)

        scatAngle = self.getScatAngle(sza,saa,vza,vaa)
        self.scatAngle = np.ma.array(scatAngle)
        
        mask = np.zeros(self.sza.shape).astype(bool)
        if self.inFilem1 is None:
            mask[0:100,:] = True
        if self.inFilep1 is None:
            mask[-100:,:] = True

        
        self.sza.mask = mask
        self.saa.mask = mask
        self.vza.mask = mask
        self.vaa.mask = mask
        self.scatAngle.mask = mask
        self.time.mask = mask
    # --
    def min_distance(self,lon0,lat0,lon,lat):
        # get index of closest lon and lat to lon0, lat0

        # use haversine formula for distance
        phi0, phi, lam0, lam = map(np.radians,[lat0,lat,lon0,lon])
        dphi = phi - phi0
        dlam = lam - lam0

        a = np.sin(dphi*0.5)**2 + np.cos(phi0)*np.cos(phi)*np.sin(dlam*0.5)**2
        if np.any(a > 1.0):
            a[a> 1.0] = 1.0
        c = 2.0*np.arcsin(np.sqrt(a))
        distance = self.Re*c

        imin = np.argmin(distance)

        return imin

    # --
    def subpointView(self):
        # figure out which satellite pixel gets closest to each satellite sub-point
        # where and when is the satellite when it points to the satellite sub-point
        # with each along track viewing angle
        self.Itymeview    = np.zeros([self.ntyme,self.nalong]).astype(int)
        self.Icrossview    = np.zeros([self.ntyme,self.nalong]).astype(int)
        p = Pool(self.nalong)
        
        lon = self.LONGITUDE[self.Istyme:self.Istyme+self.ntyme]
        lat = self.LATITUDE[self.Istyme:self.Istyme+self.ntyme]
        args = [(lon,lat,self.lon_granule[:,ialong,:],self.lat_granule[:,ialong,:],self.Re,self.Istyme) for ialong in range(self.nalong)]
        result = p.map(get_imin,args)
        for ialong,r in enumerate(result):
            self.Itymeview[:,ialong] = r[:,0]
            self.Icrossview[:,ialong] = r[:,1]
        
        p.close()
    # --
    def pixel_loc(self,lon0,lat0,vza,vaa):
        # givent satellite sub-point and viewing direction, get pixel lat/lon
        lon0, lat0 = np.radians(lon0),np.radians(lat0)
        slat0 = np.sin(lat0)
        clat0 = np.cos(lat0)

        # Get angular distance
        eta  = np.radians(vza)
        seta = np.sin(eta)
        srho = self.Re/(self.Re+self.hgtss)
        cepsilon = seta/srho
        epsilon  = np.arccos(cepsilon)
        # eta + epsilon + lambr =  90 degress
        lambr = np.pi*0.5 - eta - epsilon
        clamb = np.cos(lambr)
        slamb = np.sin(lambr)
        
        az = np.radians(vaa)
        caz = np.cos(az)
        saz = np.sin(az)
        
        slat = slat0*clamb + clat0*slamb*caz
        lat = np.degrees(np.arcsin(slat))

        lon = lon0 + np.arctan2(saz*slamb*clat0,clamb-slat0*slat)
        lon = np.degrees(lon)

        return lon,lat


    # --
    def granulePixels(self):
        # Calculate full Granule Pixel Locations and View Angles
        self.lon_granule = np.zeros([self.ntymeTotal,self.nalong,self.ncross])
        self.lat_granule = np.zeros([self.ntymeTotal,self.nalong,self.ncross])
        self.vza_granule = np.zeros([self.ntymeTotal,self.nalong,self.ncross])
        self.vaa_granule = np.zeros([self.ntymeTotal,self.nalong,self.ncross])

        srho = self.Re/(self.Re+self.hgtss)
        for ialong,along_deg in enumerate(self.vza_along):
            eta_along = np.radians(abs(along_deg))
            seta_along = np.sin(eta_along)
            cepsilon_along   = seta_along/srho
            epsilon_along    = np.arccos(cepsilon_along)
            #subtended angle in along track direction
            lambr_along = np.pi*0.5 - eta_along - epsilon_along
            for icross,cross_deg in enumerate(self.vza_cross):
                eta_cross = np.radians(abs(cross_deg))
                seta_cross = np.sin(eta_cross)
                cepsilon_cross = seta_cross/srho
                epsilon_cross = np.arccos(cepsilon_cross)
                # subtended angle in cross track direction
                lambr_cross = np.pi*0.5 - eta_cross - epsilon_cross

                # use spherical pythagorean theorem to get
                # subtended angle along view
                clambr_view = np.cos(lambr_cross)*np.cos(lambr_along)
                lambr_view = np.arccos(clambr_view)
                slambr_view = np.sin(lambr_view)

                # view zenith angle
                teta = srho*slambr_view/(1-srho*clambr_view)
                eta  = np.arctan(teta)
                self.vza_granule[:,ialong,icross] = np.degrees(eta)

                # get view azimuth angle using spherical law of cosines
                clambr_cross = np.cos(lambr_cross)
                clambr_along = np.cos(lambr_along)
                slambr_along = np.sin(lambr_along)
                caz = (clambr_cross - clambr_view*clambr_along)/(slambr_view*slambr_along)
                azr  = np.arccos(caz)
                az   = np.degrees(azr)
                if along_deg >= 0:
                    # looking forwards
                    vaa = self.VAA
                    if cross_deg >= 0:
                        vaa = vaa + az
                    else:
                        vaa = vaa - az                        
                else:
                    # looking backwards
                    vaa = self.VAA + 180.0           
                    if any(vaa > 360.0):
                        I = vaa > 360.0
                        vaa[I] = 360.0 - vaa[I]
                    if cross_deg >= 0:
                        vaa = vaa - az
                    else:
                        vaa = vaa + az

                if any(vaa > 360.0):
                    I = vaa > 360.0
                    vaa[I] = 360.0 - vaa[I]
                if any(vaa < 0):
                    I = vaa < 0
                    vaa[I] = 360.0 + vaa[I]

                self.vaa_granule[:,ialong,icross] = vaa

                # given vza and vaa, get lat/lon of pixel
                lon0 = self.LONGITUDE
                lat0 = self.LATITUDE
                vza = self.vza_granule[:,ialong,icross]
                vaa = self.vaa_granule[:,ialong,icross]
                lon,lat = self.pixel_loc(lon0,lat0,vza,vaa)

                self.lon_granule[:,ialong,icross] = lon
                self.lat_granule[:,ialong,icross] = lat

    # --
    def granuleSolarAngles(self):
        """
        use fortran code to get solar angles for granule
        """
        istart = self.Istyme 
        iend   = self.Istyme + self.ntyme

        self.sza_granule = np.zeros([self.ntyme,self.nalong,self.ncross])
        self.saa_granule = np.zeros([self.ntyme,self.nalong,self.ncross])

        p = Pool(self.nalong)
        args = [(self.tyme[istart:iend],self.lat_granule[istart:iend,ialong,:],self.lon_granule[istart:iend,ialong,:]) for ialong in range(self.nalong)]
        result = p.map(SolarAngles,args)
        for ialong,r in enumerate(result):
            self.sza_granule[:,ialong,:] = r[:,:,0]
            self.saa_granule[:,ialong,:] = r[:,:,1]
        p.close()

    # --
    def get_vza_cross(self):
        """ get approximate cross track swath full angle using 
            spherical earth approximation

            Note: use reference satellite altitude for calculations
        """

        # get cross track swath full angle
        swath_width = self.cross_track_km*0.5
        lambr = swath_width/self.Re
        slamb = np.sin(lambr)
        clamb = np.cos(lambr)
        srho = self.Re/(self.Re+self.hgtssRef)
        teta = srho*slamb/(1-srho*clamb)

        eta = np.degrees(np.arctan(teta))

        self.cross_track_deg = eta*2

        # get approximate cross track swath delta angle for resolution given
        pixel_width = self.cross_track_dkm*0.5
        lambr = pixel_width/self.Re
        slamb = np.sin(lambr)
        clamb = np.cos(lambr)
        srho = self.Re/(self.Re+self.hgtssRef)
        teta = srho*slamb/(1-srho*clamb)

        eta = 2.*np.degrees(np.arctan(teta))

        # get number of cross track pixels
        # round to whole number
        ncross = int(self.cross_track_deg/eta) + 1
        deta   = self.cross_track_deg/(ncross-1)

        self.ncross = ncross
        self.vza_cross = -1.*self.cross_track_deg*0.5 + np.arange(ncross)*deta
    # --
    def nadirAngles(self):
        SZA   = []
        SAA   = []
        VAA   = []
        
        for i,tyme in enumerate(self.tyme):
            CLAT = self.LATITUDE[i]
            CLON = self.LONGITUDE[i]
            
            SLAT = self.LATITUDE[i]
            SLON = self.LONGITUDE[i]
            sat_angles = LidarAngles_.satangles(tyme.year,tyme.month,tyme.day,
                                                tyme.hour,tyme.minute,tyme.second,
                                                CLAT,CLON,
                                                SLAT,SLON,
                                                0.0,
                                                self.hgtss)
            SZA.append(sat_angles[3][0])
            SAA.append(sat_angles[2][0])

            # get local azimuth angle for direction satellite is moving in
            if i == self.ntymeTotal-1:
                CLAT = self.LATITUDE[i-1]
                CLON = self.LONGITUDE[i-1]
            else:
                CLAT = self.LATITUDE[i+1]
                CLON = self.LONGITUDE[i+1]
            sat_angles = LidarAngles_.satangles(tyme.year,tyme.month,tyme.day,
                                                tyme.hour,tyme.minute,tyme.second,
                                                CLAT,CLON,
                                                SLAT,SLON,
                                                0.0,
                                                self.hgtss)

            if i == self.ntymeTotal-1:
                # add 180 to look forward
                vaa = sat_angles[0][0] + 180.0
                if vaa >=360:
                    vaa = vaa - 360.0

                VAA.append(vaa)
            else:
                VAA.append(sat_angles[0][0])


        self.SZA = np.array(SZA)
        self.SAA = np.array(SAA)
        self.VAA = np.array(VAA)  

    #---
    def readSampledGEOS(self):
        """
        Read in model sampled track
        """
        col = 'aer_Nv'
        filelist = []
        if self.inFilem1 is not None:
            filelist.append(self.inFilem1)

        filelist.append(self.inFile)
        if self.inFilep1 is not None:
            filelist.append(self.inFilep1)

        for inFile in filelist:
            if self.verbose: 
                print 'opening file',inFile.replace('%col',col)
            nc       = Dataset(inFile.replace('%col',col))

            for sds in self.SDS_AER:
                sds_ = sds
                if sds in ncALIAS:
                    sds_ = ncALIAS[sds]
                var = nc.variables[sds_][:]
                self.__dict__[sds].append(var)



#------------------------------------ M A I N ------------------------------------

if __name__ == "__main__":

    # Defaults
    DT_hours = 1

#   Parse command line options
#   --------------------------

    parser = argparse.ArgumentParser()
    parser.add_argument("iso_t1",
                        help="starting iso time")
    parser.add_argument("iso_t2",
                        help="ending iso time")

    parser.add_argument("track_pcf",
                        help="prep config file with track input file names")

    parser.add_argument("orbit_pcf",
                        help="prep config file with orbit variables")

    parser.add_argument("inst_pcf",
                        help="prep config file with instrument variables")

    parser.add_argument("-D","--DT_hours", default=DT_hours, type=int,
                        help="Timestep in hours for each file (default=%i)"%DT_hours)

    parser.add_argument("-v", "--verbose",action="store_true",
                        help="Verbose mode (default=False).")

    parser.add_argument("-r", "--dryrun",action="store_true",
                        help="do a dry run (default=False).")    

    args = parser.parse_args()

    # Parse prep config
    # -----------------
    cf             = Config(args.track_pcf,delim=' = ')
    inTemplate     = cf('inDir')     + '/' + cf('inFile')         

    cf             = Config(args.orbit_pcf,delim=' = ')
    orbitname      = cf('orbitname')
    ORBITNAME      = orbitname.upper()
    HGT            = float(cf('HGT'))

    cf             = Config(args.inst_pcf,delim=' = ')
    instname       = cf('instname')

    try:
        cross_track_km = float(cf('cross_track_km'))
        cross_track_dkm = float(cf('cross_track_dkm'))
    except:
        cross_track_km = None
        cross_track_dkm = None

    try: 
        along_track_deg = cf('along_track_deg')
    except:
        along_track_deg = None


    # Loop through dates, calculating geometry
    # ------------------------------------
    date      = isoparser(args.iso_t1)
    enddate   = isoparser(args.iso_t2)
    Dt        = timedelta(hours=args.DT_hours)

    while date < enddate:
        nymd  = str(date.date()).replace('-','')
        year  = str(date.year)
        month = str(date.month).zfill(2)
        day   = str(date.day).zfill(2)
        hour  = str(date.hour).zfill(2)    

        inFile     = inTemplate.replace('%year',year).replace('%month',month).replace('%day',day).replace('%nymd',nymd).replace('%hour',hour).replace('%orbitname',orbitname).replace('%ORBITNAME',ORBITNAME)
        outFile    = inTemplate.replace('%col',instname).replace('%year',year).replace('%month',month).replace('%day',day).replace('%nymd',nymd).replace('%hour',hour).replace('%orbitname',orbitname).replace('%ORBITNAME',ORBITNAME)

        datem1 = date - Dt
        if datem1 < datetime(2006,01,01,00):
            inFilem1 = None
        else:
            nymd  = str(datem1.date()).replace('-','')
            year  = str(datem1.year)
            month = str(datem1.month).zfill(2)
            day   = str(datem1.day).zfill(2)
            hour  = str(datem1.hour).zfill(2)    
            inFilem1   = inTemplate.replace('%year',year).replace('%month',month).replace('%day',day).replace('%nymd',nymd).replace('%hour',hour).replace('%orbitname',orbitname).replace('%ORBITNAME',ORBITNAME)

        datep1 = date + Dt
        if datep1 >= datetime(2007,01,01,00):
            inFilep1 = None
        else:
            nymd  = str(datep1.date()).replace('-','')
            year  = str(datep1.year)
            month = str(datep1.month).zfill(2)
            day   = str(datep1.day).zfill(2)
            hour  = str(datep1.hour).zfill(2)    
            inFilep1   = inTemplate.replace('%year',year).replace('%month',month).replace('%day',day).replace('%nymd',nymd).replace('%hour',hour).replace('%orbitname',orbitname).replace('%ORBITNAME',ORBITNAME)


        # Initialize SWATH class and create outfile
        # -----------------------------------------------------------
        print '++++Getting polarimeter swath with the following arguments+++'
        print '>>>inFile:    ',inFile
        print '>>>outFile:   ',outFile
        print '>>>HGT:       ',HGT
        print '>>>verbose:   ',args.verbose
        print '++++End of arguments+++'
        if not args.dryrun:
            swath = SWATH(date,inFilem1,inFile,inFilep1,outFile,HGT,along_track_deg,
                          cross_track_km=cross_track_km,
                          cross_track_dkm=cross_track_dkm)

        date += Dt
