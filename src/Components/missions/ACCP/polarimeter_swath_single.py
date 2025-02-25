#!/usr/bin/env python3

"""
    Gets geometry for the polarimeter swath

    Patricia Castellanos, Dec, 2019

"""

import os
import argparse
from   datetime        import datetime, timedelta
from   dateutil.parser import parse         as isoparser
from   pyobs.sgp4      import getTrack as getTrackTLE
from   MAPL.config     import Config
import numpy  as np
from netCDF4 import Dataset
from  py_leo_vlidort import LidarAngles_ 
from mpl_toolkits.basemap import Basemap
from   py_leo_vlidort.copyvar  import _copyVar
from multiprocessing import Pool
from scipy import signal

# Generic Lists of Varnames and Units
ncALIAS = {'LONGITUDE': 'trjLon',
           'LATITUDE': 'trjLat'}

#radius of the earth km
Re = 6378.0
hgtssRef = 500.


# --
def distance(lon0,lat0,lon,lat,Re):
    # get index of closest lon and lat to lon0, lat0

    # use haversine formula for distance
    phi0, phi, lam0, lam = list(map(np.radians,[lat0,lat,lon0,lon]))
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

def get_az_vna_vza_lambrview(args):
    srho,along_deg,vna_cross,ntymeTotal,HEAD = args
    ncross = len(vna_cross)

    # azimuth of satellite measured from North at SSP
    vna_granule = np.zeros([ntymeTotal,ncross])
    az_granule  = np.zeros([ntymeTotal,ncross])

    # view zenith angle
    vza_granule  = np.zeros([ntymeTotal,1,ncross])

    # subtended angle along view
    lambr_view   = np.zeros([ncross])

    eta_along = np.radians(abs(along_deg))
    seta_along = np.sin(eta_along)
    cepsilon_along   = seta_along/srho
    epsilon_along    = np.arccos(cepsilon_along)
    #subtended angle in along track direction
    lambr_along = np.pi*0.5 - eta_along - epsilon_along

    for icross,cross_deg in enumerate(vna_cross):
        eta_cross = np.radians(abs(cross_deg))
        seta_cross = np.sin(eta_cross)
        cepsilon_cross = seta_cross/srho
        epsilon_cross = np.arccos(cepsilon_cross)
        # subtended angle in cross track direction
        lambr_cross = np.pi*0.5 - eta_cross - epsilon_cross

        # use spherical pythagorean theorem to get
        # subtended angle along view
        clambr_view = np.cos(lambr_cross)*np.cos(lambr_along)
        lambr_view[icross] = np.arccos(clambr_view)
        slambr_view = np.sin(lambr_view[icross])

        # view nadir angle
        teta_view = srho*slambr_view/(1-srho*clambr_view)
        eta_view  = np.arctan(teta_view)
        vna_granule[:,icross] = np.degrees(eta_view)

        #view zenith angle
        seta_view = np.sin(eta_view)
        cepsilon_view = seta_view/srho
        epsilon_view = np.arccos(cepsilon_view)
        vza_granule[:,0,icross] = 90. - np.degrees(epsilon_view)

        # get azimuth angle using spherical law of cosines
        clambr_cross = np.cos(lambr_cross)
        clambr_along = np.cos(lambr_along)
        slambr_along = np.sin(lambr_along)
        if slambr_along == 0.0:
            # special case for nadir angle
            caz = 0.0
        else:
            caz = (clambr_cross - clambr_view*clambr_along)/(slambr_view*slambr_along)
    
        if caz > 1.0:
            caz = 1.0
        azr  = np.arccos(caz)
        az   = np.degrees(azr)
        if along_deg >= 0:
            # looking forwards
            aa = HEAD
            if cross_deg >= 0:
                aa = aa + az
            else:
                aa = aa - az
        else:
            # looking backwards
            aa = HEAD + 180.0
            if any(aa > 360.0):
                I = aa > 360.0
                aa[I] = aa[I] - 360.0
            if cross_deg >= 0:
                aa = aa - az
            else:
                aa = aa + az
        # make sure between 0-360
        if any(aa > 360.0):
            I = aa > 360.0
            aa[I] = aa[I] - 360.0
        if any(aa < 0):
            I = aa < 0
            aa[I] = 360.0 + aa[I]

        az_granule[:,icross] = aa

    return az_granule,vna_granule,vza_granule,lambr_view
    # --

def get_vaa(args):
    ntymeTotal,ncross,longitude,latitude,lon0,lat0,lambr_view = args
    # use law of cosines on spherical triangle formed by 
    # target, north pole and satellite subpoint
    vaa_granule = np.zeros([ntymeTotal,1,ncross])

    for icross in range(ncross):
        lon = longitude[:,icross]
        lat = latitude[:,icross]
    
        # target
        latr = np.radians(lat)
        clatp = np.cos(latr)
        slatp = np.sin(latr)

        # satellite subpoint
        latr = np.radians(lat0)
        slatss = np.sin(latr)

        # subtended angle along view
        clambr_view = np.cos(lambr_view[icross])
        slambr_view = np.sin(lambr_view[icross])

        if slambr_view == 0.0:
            # special case for nadir angle
            cvaar = np.ones(ntymeTotal)        
        else:
            cvaar = (slatss - slatp*clambr_view)/(clatp*slambr_view)
        if any(cvaar > 1):
            I = cvaar > 1
            cvaar[I] = 1.0
        if any(cvaar < -1):
            I = cvaar < -1
            cvaar[I] = -1.0
        vaa = np.degrees(np.arccos(cvaar))
        # where is satellite subpoint to the west of target
        # these VAA's are negative
        I = lon > lon0
        if any(I):
            vaa[I] = 360.0 - vaa[I]
        vaa_granule[:,0,icross] = vaa
    return vaa_granule

class SWATH(object):
    """
    get polarimeter swath
    """
    def __init__(self,starttyme,endtyme,dt_secs,trjFile,outFile,hgtss,along_track_deg,
                 verbose=True,
                 cross_track_km=None,cross_track_dkm=None,
                 no_ss=False,
                 single=False):
        self.starttyme = starttyme
        self.endtyme   = endtyme
        self.dt_secs   = dt_secs
        self.dt        = timedelta(seconds=dt_secs)
        self.trjFile   = trjFile
        self.outFile  = outFile
        self.hgtss    = hgtss
        self.verbose  = verbose
        self.cross_track_km = cross_track_km
        self.cross_track_dkm = cross_track_dkm
        self.Re       = Re
        self.hgtssRef = hgtssRef
        self.no_ss    = no_ss
        self.single   = single

        # Parse along_track_deg variable
        # vna == view nadir angle.  view angle between satellite subpoint and target at the satellite 
        if ',' in along_track_deg:
            # list
            self.vna_along = np.array(along_track_deg.split(',')).astype(float)
            self.nalong    = len(self.vna_along)
        elif ':' in along_track_deg:
            vnamin,vnamax,nalong = np.array(along_track_deg.split(':')).astype(float)

            self.vna_along = np.linspace(vnamin,vnamax,nalong)
            self.nalong    = int(nalong)

        # calculate trajectory from TLE file
        t1 = starttyme
        t2 = endtyme
        if self.single:
            t1,t2 = self.get_single_orbit(trjFile,t1,t2,dt_secs)
            print('new t1,t2',t1,t2)

        LONGITUDE, LATITUDE, tyme = getTrackTLE(trjFile, t1, t2, dt_secs)
        self.ntyme = len(LONGITUDE)

        isotime = []
        for tt in tyme:
            iso = list(tt.isoformat())
            isotime.append(iso)

        self.isotime = np.array(isotime)

        # get a few minutes prior to start
        t1 = starttyme-200*dt
        t2 = starttyme-self.dt
        lon, lat, ptyme = getTrackTLE(trjFile, t1, t2, dt_secs)
        self.Istyme = len(ptyme)
        LONGITUDE = np.append(lon,LONGITUDE)
        LATITUDE  = np.append(lat,LATITUDE)
        tyme      = np.append(np.array(ptyme),tyme)

        # get a few minutes after start
        t1 = endtyme+self.dt
        t2 = endtyme+200*dt
        lon, lat, atyme = getTrackTLE(trjFile, t1, t2, dt_secs)
        LONGITUDE = np.append(LONGITUDE,lon)
        LATITUDE  = np.append(LATITUDE,lat)
        tyme      = np.append(tyme,np.array(atyme))

        self.tyme = tyme
        self.LONGITUDE = LONGITUDE
        self.LATITUDE  = LATITUDE
        self.ntymeTotal = len(self.tyme)

        # Get Satellite Heading
        self.hgtss = hgtss
        self.nadirAngles()

        # Get cross track swath deg full angle
        if self.cross_track_km is not None:
            self.get_vna_cross()

        # Calculate Full Pixel Locations and view angles
        self.granulePixels()

        # Calculate Full Pixel solar angles
        self.granuleSolarAngles()

        # Calculate scattering angle
        istart = self.Istyme 
        iend   = self.Istyme + self.ntyme
        self.scatAngle_granule = self.getScatAngle(self.sza_granule[istart:iend,:,:],self.saa_granule[istart:iend,:,:],self.vza_granule[istart:iend,:,:],self.vaa_granule[istart:iend,:,:])

        if self.no_ss is False:
            # figure out when satellite pixels overlap with each satellite sub-point
            self.subpointView()

            # Put together arrays of angles for when 
            # satellite points at the satellite sub-point
            self.subpointViewAngles()

        # write to file
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

        # Create dimensions
        # -----------------
        nt = nc.createDimension('time',self.ntyme)
        ls = nc.createDimension('ls',19)
        x  = nc.createDimension('x',1)
        y  = nc.createDimension('y',1)
        nalong = nc.createDimension('nalong',self.nalong)
        ncross = nc.createDimension('ncross',self.ncross)

        # Coordinate variables
        # --------------------
        time = nc.createVariable('time','i4',('time',),zlib=True)
        time.long_name = 'Time'
        t0 = self.starttyme
        time.units = 'seconds since %s'%t0.isoformat(' ')
        tyme = self.tyme[self.Istyme:self.Istyme+self.ntyme]
        time[:] = np.array([(t-t0).total_seconds() for t in tyme])

        # Add fake dimensions for GrADS compatibility
        # -------------------------------------------
        x = nc.createVariable('x','f4',('x',),zlib=True)
        x.long_name = 'Fake Longitude for GrADS Compatibility'
        x.units = 'degrees_east'
        x[:] = np.zeros(1)
        y = nc.createVariable('y','f4',('y',),zlib=True)
        y.long_name = 'Fake Latitude for GrADS Compatibility'
        y.units = 'degrees_north'
        y[:] = np.zeros(1)

        # Trajectory coordinates
        # ----------------------
        lon = nc.createVariable('trjLon','f4',('time',),zlib=True)
        lon.long_name = 'Trajectory Longitude'
        lon.units = 'degrees_east'
        lon[:] = self.LONGITUDE[self.Istyme:self.Istyme+self.ntyme]
        lat = nc.createVariable('trjLat','f4',('time',),zlib=True)
        lat.long_name = 'Trajectory Latitude'
        lat.units = 'degrees_north'
        lat[:] = self.LATITUDE[self.Istyme:self.Istyme+self.ntyme]


        # Time in ISO format
        # ------------------
        isotime = nc.createVariable('isotime','S1',('time','ls'),zlib=True)
        isotime.long_name = 'Time (ISO Format)'
        isotime[:] = self.isotime

        # Create Variables
        # ------------------
        if self.no_ss is False:
            dim = ('time','nalong',)
            this = nc.createVariable('time_ss','i4',dim,zlib=True)
            this[:] = self.time
            this.long_name = 'time when polarimeter views satellite subpoint'
            this.units = 'seconds since ' + self.starttyme.isoformat()

            this = nc.createVariable('sza_ss','f4',dim,zlib=True)
            this[:] = self.sza
            this.long_name = 'SZA when polarimeter views satellite subpoint'
            this.units = 'degrees, >90 is below horizon'

            this = nc.createVariable('saa_ss','f4',dim,zlib=True)
            this[:] = self.saa
            this.long_name = 'SAA when polarimeter views satellite subpoint'
            this.units = 'degrees, clockwise from north 0-360'

            this = nc.createVariable('vza_ss','f4',dim,zlib=True)
            this[:] = self.vza
            this.long_name = 'VZA when polarimeter views satellite subpoint'
            this.units = 'degrees'

            this = nc.createVariable('vaa_ss','f4',dim,zlib=True)
            this[:] = self.vaa
            this.long_name = 'VAA when polarimeter views satellite subpoint'
            this.units = 'degrees, clockwise from north 0-360'

            this = nc.createVariable('scatAngle_ss','f4',dim,zlib=True)
            this[:] = self.scatAngle
            this.long_name = 'scattering angle when polarimeter view satellite subpoint'
            this.units = 'degrees, 0 = forward scattering'

        dim = ('time','nalong','ncross',)
        istart = self.Istyme 
        iend   = self.Istyme + self.ntyme

        this = nc.createVariable('longitude','f4',dim,zlib=True)
        this[:] = self.lon_granule[istart:iend,:,:]
        this.long_name = 'longitude of granule pixel centers'
        this.units  = 'degrees_east'

        this = nc.createVariable('latitude','f4',dim,zlib=True)
        this[:] = self.lat_granule[istart:iend,:,:]
        this.long_name = 'latitude of granule pixel centers'
        this.units = 'degrees_north'

        this = nc.createVariable('vza','f4',dim,zlib=True)
        this[:] = self.vza_granule[istart:iend,:,:]
        this.long_name = 'view zenith angle'
        this.units = 'degrees'

        this = nc.createVariable('vaa','f4',dim,zlib=True)
        this[:] = self.vaa_granule[istart:iend,:,:]
        this.long_name = 'view azimuth angle'
        this.units = 'degrees, clockwise from north 0-360'

        this = nc.createVariable('sza','f4',dim,zlib=True)
        this[:] = self.sza_granule[istart:iend,:,:]
        this.long_name = 'solar zenith angle'
        this.units = 'degrees'

        this = nc.createVariable('saa','f4',dim,zlib=True)
        this[:] = self.saa_granule[istart:iend,:,:]
        this.long_name = 'solar azimuth angle'
        this.units = 'degrees, clockwise from north 0-360'

        this = nc.createVariable('scatAngle','f4',dim,zlib=True)
        this[:] = self.scatAngle_granule
        this.long_name = 'scattering angle'
        this.units = 'degrees, 0 = forward scattering'        

        nc.close()

# ---
    def get_single_orbit(self,trjFile,t1,t2,dt_secs):
        """
        pick out one orbit closest to prime meridion
        """
        # calculate orbit for full day
        CLON, CLAT, tyme = getTrackTLE(trjFile, t1, t2, dt_secs)
        year = []
        month = []
        day = []
        hour = []
        minute = []
        second = []
        for t in tyme:
            year.append(t.year)
            month.append(t.month)
            day.append(t.day)
            hour.append(t.hour)
            minute.append(t.minute)
            second.append(t.second)


        # calculate sza time series
        saa,sza = LidarAngles_.solarangles(year,month,day,
                                                hour,minute,second,
                                                CLAT,CLON,
                                                0.0)

        # get all the peaks in the time series of sza
        # these are the orbit start indeces
        a = signal.find_peaks(sza,width=10)
        a = a[0]

        # find the orbit that is close to prime meridion
        peak = None
        for i in range(len(a)-1):
            istart = a[i]
            iend   = a[i+1]
            llon = CLON[istart:iend]
            llat = CLAT[istart:iend]
            ssza = sza[istart:iend]
            if any((llon<10) & (llon>-10) & (llat<10) & (llat>-10)& (ssza<90)):
                peak = i

        if peak is None:
            # try again
            for i in range(len(a)-1):
                istart = a[i]
                iend   = a[i+1]
                llon = lon[istart:iend]
                llat = lat[istart:iend]
                ssza = sza[istart:iend]
                if any((llon<17) & (llon>-17) & (llat<10) & (llat>-10)& (ssza<90)):
                    peak = i

        if peak is None:
            print("coulnd't find peak! Exiting")
            print(date.isoformat())
            sys.exit()

        t1_single = tyme[a[peak]]
        t2_single = tyme[a[peak+1]]

        return t1_single,t2_single
# ---
    def getScatAngle(self,sza,saa,vza,vaa):
        """
        The angle between the sun, the pixel and sensor.  A value of
        180.0 would imply pure backscatter.  A value of zero would
        imply pure forward scatter (not seen in satellite remote sensing).
        Values less than 90 imply forward scattering.
        Range is 0 to 180.0 
        """       
    
        raa = np.abs(vaa - saa)
        I = raa > 180
        raa[I] = 360.0 - raa[I]

        szar = np.radians(sza)
        vzar = np.radians(vza)
        raar = np.radians(raa)
        coscattering_angle = np.cos(szar)*np.cos(vzar) + np.sin(szar)*np.sin(vzar)*np.cos(raar)
        if np.any(coscattering_angle) > 1.0:
            I = coscattering_angle > 1.0
            coscattering_angle[I] = 1.0

        if np.any(coscattering_angle) < -1.0:
            I = coscattering_angle < -1.0
            coscattering_angle[I] = -1.0

        scattering_angle = 180. - np.degrees(np.arccos(coscattering_angle))

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
            dtyme = self.tyme[self.Itymeview[:,ialong]] - self.starttyme
            for ityme in range(self.ntyme):
                sza[ityme,ialong] = self.sza_granule[:,ialong,:][self.Itymeview[ityme,ialong],self.Icrossview[ityme,ialong]]
                saa[ityme,ialong] = self.saa_granule[:,ialong,:][self.Itymeview[ityme,ialong],self.Icrossview[ityme,ialong]]
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
        
    # --
    def min_distance(self,lon0,lat0,lon,lat):
        # get index of closest lon and lat to lon0, lat0

        # use haversine formula for distance
        phi0, phi, lam0, lam = list(map(np.radians,[lat0,lat,lon0,lon]))
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
        p.join()
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
        """
        Calculate full Granule Pixel Locations and Angles
        """

        self.lon_granule = np.zeros([self.ntymeTotal,self.nalong,self.ncross])
        self.lat_granule = np.zeros([self.ntymeTotal,self.nalong,self.ncross])
        self.vza_granule = np.zeros([self.ntymeTotal,self.nalong,self.ncross])
        self.vaa_granule = np.zeros([self.ntymeTotal,self.nalong,self.ncross])

        # get vna and az - azimuth of satellite measured from North at SSP
        srho = self.Re/(self.Re+self.hgtss)
        p = Pool(self.nalong)
        args = []
        for along_deg in self.vna_along:
            args.append([srho,along_deg,self.vna_cross,self.ntymeTotal,self.HEAD])

        result = p.map(get_az_vna_vza_lambrview,args)
        p.close()
        p.join()

        az_granule,vna_granule,vza_granule,lambr_view = [],[],[],[]
        for r in result:
            az,vna,vza,lambr = r
            az_granule.append(az)
            vna_granule.append(vna)
            vza_granule.append(vza)
            lambr_view.append(lambr)

        self.vza_granule = np.concatenate(vza_granule,axis=1)

        # given vna and az, get lat/lon of pixel
        lon0 = self.LONGITUDE
        lat0 = self.LATITUDE
        for ialong,along_deg in enumerate(self.vna_along):            
            for icross,cross_deg in enumerate(self.vna_cross):
                # given vna and az, get lat/lon of pixel
                vna  = vna_granule[ialong][:,icross]
                az   = az_granule[ialong][:,icross]
                lon,lat = self.pixel_loc(lon0,lat0,vna,az)

                self.lon_granule[:,ialong,icross] = lon
                self.lat_granule[:,ialong,icross] = lat

        # get azimuth of satellite measure from North at target
        p = Pool(self.nalong)
        args = []
        for ialong in range(self.nalong):
            args.append([self.ntymeTotal,self.ncross,self.lon_granule[:,ialong,:],self.lat_granule[:,ialong,:],lon0,lat0,lambr_view[ialong]])

        vaa_granule = p.map(get_vaa,args)
        p.close()
        p.join()

        self.vaa_granule = np.concatenate(vaa_granule,axis=1)
        
    # --
    def granuleSolarAngles(self):
        """
        use fortran code to get solar angles for granule
        """
        self.sza_granule = np.zeros([self.ntymeTotal,self.nalong,self.ncross])
        self.saa_granule = np.zeros([self.ntymeTotal,self.nalong,self.ncross])

        p = Pool(self.nalong)
        args = [(self.tyme,self.lat_granule[:,ialong,:],self.lon_granule[:,ialong,:]) for ialong in range(self.nalong)]
        result = p.map(SolarAngles,args)
        for ialong,r in enumerate(result):
            self.sza_granule[:,ialong,:] = r[:,:,0]
            self.saa_granule[:,ialong,:] = r[:,:,1]
        p.close()
        p.join()
    # --
    def get_vna_cross(self):
        """ 
        get approximate cross track swath view angle using 
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
        self.vna_cross = -1.*self.cross_track_deg*0.5 + np.arange(ncross)*deta
    # --
    def nadirAngles(self):
        """
        Get the satellite heading (viewing azimuth angle at satellite subpoint)
        """
        VAA   = []
        
        for i,tyme in enumerate(self.tyme):
            SLAT = self.LATITUDE[i]
            SLON = self.LONGITUDE[i]
        
            # get local azimuth angle for direction satellite is moving in
            if i == self.ntymeTotal-1:
                CLAT = self.LATITUDE[i-1]
                CLON = self.LONGITUDE[i-1]
            else:
                CLAT = self.LATITUDE[i+1]
                CLON = self.LONGITUDE[i+1]
            sat_angles = LidarAngles_.satangles_fromspace(tyme.year,tyme.month,tyme.day,
                                                tyme.hour,tyme.minute,tyme.second,
                                                CLAT,CLON,
                                                SLAT,SLON,
                                                0.0,
                                                self.hgtss)
            vaa = sat_angles[0][0]
            # check if going from +180 to -180
            if (SLON > 179) and (CLON < 0):
                vaa = 360. - vaa
                print('here',SLON,CLON,vaa,i)
            # check if going from -180 to +180
            if (SLON < -179) and (CLON > 0):
                vaa = 360. - vaa
                print('here again',SLON,CLON,vaa,i)
            if i == self.ntymeTotal-1:
                # add 180 to look forward
                vaa = vaa + 180.0
                if vaa >=360:
                    vaa = vaa - 360.0

                VAA.append(vaa)
            else:
                VAA.append(vaa)


        self.HEAD = np.array(VAA)  

#------------------------------------ M A I N ------------------------------------

if __name__ == "__main__":

    # Defaults
    DT_hours = 1
    dt_secs  = 1

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

    parser.add_argument("-d","--dt_secs", default=dt_secs, type=int,
                        help="Timestep in seconds for the trajectory (default=%i)"%dt_secs)

    parser.add_argument("-D","--DT_hours", default=DT_hours, type=int,
                        help="Timestep in hours for each file (default=%i)"%DT_hours)

    parser.add_argument("-v", "--verbose",action="store_true",
                        help="Verbose mode (default=False).")

    parser.add_argument("-r", "--dryrun",action="store_true",
                        help="do a dry run (default=False).")    

    parser.add_argument("--no_ss",action="store_true",
                        help="don't do satellite subpoint aggregation (default=False).")

    parser.add_argument("--single",action="store_true",
                        help="only calculate 1 orbit per day (default=False).")

    args = parser.parse_args()

    # Parse prep config
    # -----------------
    cf             = Config(args.track_pcf,delim=' = ')
    inTemplate     = cf('inDir')     + '/' + cf('inFile')         

    cf             = Config(args.orbit_pcf,delim=' = ')
    orbitname      = cf('orbitname')
    ORBITNAME      = orbitname.upper()
    HGT            = float(cf('HGT'))
    trjFile        = cf('TLEFILE')

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
    dt        = timedelta(seconds=args.dt_secs)
    if args.single:
        Dt        = timedelta(hours=24)
    else:
        Dt        = timedelta(hours=args.DT_hours)

    while date < enddate:
        edate = date + Dt - dt
        nymd  = str(date.date()).replace('-','')
        year  = str(date.year)
        month = str(date.month).zfill(2)
        day   = str(date.day).zfill(2)
        hour  = str(date.hour).zfill(2)    

        outFile    = inTemplate.replace('%col',instname).replace('%year',year).replace('%month',month).replace('%day',day).replace('%nymd',nymd).replace('%hour',hour).replace('%orbitname',orbitname).replace('%ORBITNAME',ORBITNAME)


        # Initialize SWATH class and create outfile
        # -----------------------------------------------------------
        print('++++Getting polarimeter swath with the following arguments+++')
        print('>>>starttyme: ',date.isoformat())
        print('>>>outFile:   ',outFile)
        print('>>>HGT:       ',HGT)
        print('>>>TLEFILE:   ',trjFile)
        print('>>>verbose:   ',args.verbose)
        print('>>>no_ss:     ',args.no_ss)
        print('++++End of arguments+++')
        if not args.dryrun:
            swath = SWATH(date,edate,args.dt_secs,trjFile,outFile,HGT,along_track_deg,
                          cross_track_km=cross_track_km,
                          cross_track_dkm=cross_track_dkm,
                          no_ss=args.no_ss,
                          single=args.single)
            swath = None
        date += Dt
