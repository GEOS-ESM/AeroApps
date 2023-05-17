"""
Classes for Handling D-AQ merger product.
"""

import xlrd
from numpy import array, NaN, size, ones, sort
from datetime import datetime, timedelta

MISSING = -9999
SITE_NAMES = ['Beltsville    ',
              'Padonia       ',
              'Fairhill      ',
              'Aldino        ',
              'Edgewood      ',
              'Essex         ',
              'Unknown       ',
              'Chesepeake Bay']


class xlsTable:

    def __init__(self,filename,skipRows=0,cf_meta=False):
        """
        Reads a Excell spreadsheet with named columns. The colums are
        returned as arrays. On input,
          skipRows  ---  how many rows to skip til variable names
          cf_meta   ---  whether CF metadata is available
        """
        book = xlrd.open_workbook(filename)
        self.sheet = book.sheet_by_index(0)
        if cf_meta:
            self.Long_name = dict()
            self.Std_name = dict()
            self.Units = dict()
        for j in range(self.sheet.ncols):
            var = str(self.sheet.row(skipRows)[j].value)
            if cf_meta:
                self.Long_name[var] = str(self.sheet.row(skipRows+1)[j].value)
                self.Std_name[var] = str(self.sheet.row(skipRows+2)[j].value)
                self.Units[var] = str(self.sheet.row(skipRows+3)[j].value)
            self.__dict__[var] = array(self.sheet.col_values(j,start_rowx=skipRows+3*cf_meta+1))

class P3_Profiles:

    def __init__(self, column_fn  = 'DISCOVER-AQ_vertical_column_estimates_R0.xls',
                       profile_fn = 'DISCOVER-AQ_vertical_profiles_R0.xls',
                       year=2011,Verbose=True):

        C = xlsTable(column_fn,skipRows=8)
        P = xlsTable(profile_fn,cf_meta=True)

        # Coordinate variables indexed by profiles
        # ----------------------------------------
        Tyme_start, Tyme_stop, Tyme = ([],[],[])
        for date,doy,t1,t2 in zip(C.Date,C.Jday,C.UTC_start,C.UTC_stop):
            today = datetime(year,1,1) + timedelta(seconds=(doy-1)*86400)
            Tyme_start += [today + timedelta(seconds=t1),]
            Tyme_stop += [today + timedelta(seconds=t2),]
            Tyme += [today + timedelta(seconds=int((t1+t2)/2.)),]

        # Coordinate variables
        # --------------------
        self.Profile = C.ProflSeqN
        self.Altitude = P.Nominal_alt[P.PrflSeqNum==1001]

        # Other flight information
        # ------------------------
        self.Longitude = C.Long_avg
        I = (self.Longitude>180.)
        self.Longitude[I] = self.Longitude[I] - 360. # make sure longitude in [-180,180]
        self.Latitude = C.Lat_avg
        self.Time_start = array(Tyme_start)
        self.Time = array(Tyme)
        self.Time_stop = array(Tyme_stop)
        self.Flight = C.Flight
        self.Site = C.SiteID

        self.Long_name = P.Long_name
        self.Std_name = P.Std_name
        self.Units = P.Units

        # Fundamental dimensions
        # ----------------------
        self.Nz = size(self.Altitude)
        self.Np = size(self.Profile)

        # Create 2D arrays for profile data
        # ---------------------------------
        exclude = ( 'PrflSeqNum', 'Nominal_alt', 'Units', 'sheet', 'Std_name', 'Long_name')
        self.Variables = dict()
        for var in P.__dict__:
            if var in exclude:
                continue
            if Verbose:
                print("<> Processing %s"%var)
            V = MISSING * ones((self.Np,self.Nz))
            for i in range(self.Np):
                p = self.Profile[i]
                I = (P.PrflSeqNum==p)
                V[i,:] = P.__dict__[var][I]
            V[V==MISSING] = NaN
            self.Variables[var] = V

#---
    def writeNC ( self, filename, format='NETCDF3_CLASSIC', zlib=False, Verbose=True):
        """
        Write a NetCDF file with the P3 profiles.
        """
        from netCDF4 import Dataset

        nc = Dataset(filename,'w',format=format)

        # Set global attributes
        # ---------------------
        nc.title = 'DISCOVER-AQ Vertical Profiles'
        nc.institution = 'NASA'
        nc.source = 'Instruments on P3B Aircraft'
        nc.history = 'CF compliant NetCDF file created from spreasheet DISCOVER-AQ_vertical_profiles_R0.xlsx'
        nc.references = 'n/a'
        nc.comment = 'This file contains 0.1 km gridded vertical profiles for temperature, potential temperature, specific humidity, O3, CO, CO2, NO, NO2, NOy, total dry scattering (550nm), total dry absorption (532nm), CN (>3nm), submicron volume density from UHSAS, black carbon mass concentration, and water-soluble organic carbon.  Variables are indexed by (profile,altitude).  The *altitude* variable gives the nominal altitude at the center of the bin.  The altitude for each species/parameter (e.g., O3_alt) is the average (within the bin) of the altitude where observations were made.   For each trace gas, 4 data column are provide for the average, standard deviation, number of the points, and average number density.  The number density column is not given for other parameters.'
        nc.contact = 'Original Data: Gao Chen <gao.chen@nasa.gov>, CF File: Arlindo da Silva <arlindo.dasilva@nasa.gov>'
        nc.Conventions = 'CF-1.5'
        nc.featureType = 'profile'

        # Create dimensions
        # -----------------
        p = nc.createDimension('profile',self.Np)
        z = nc.createDimension('altitude',self.Nz)
        s = nc.createDimension('site',len(SITE_NAMES))
        l = nc.createDimension('site_strlen',len(SITE_NAMES[0]))

        # Coordinate variables
        # --------------------
        profile = nc.createVariable('profile','f4',('profile',))
        profile.cf_role = 'profile id'
        profile[:] = self.Profile[:]
        
        altitude = nc.createVariable('altitude','f4',('altitude',))
        altitude.standard_name = 'altitude'
        altitude.long_name = 'Altitude above sea-level'
        altitude.units = 'km'
        altitude.positive = 'up'
        altitude.axis = 'z'
        altitude[:] = self.Altitude[:]
        
        longitude = nc.createVariable('longitude','f4',('profile',))
        longitude.standard_name = 'longitude'
        longitude.long_name = 'Average longitude of profile'
        longitude.units = 'degrees_east'
        longitude[:] = self.Longitude[:]

        latitude = nc.createVariable('latitude','f4',('profile',))
        latitude.standard_name = 'latitude'
        latitude.long_name = 'Average latitude of profile'
        latitude.units = 'degrees_north'
        latitude[:] = self.Latitude[:]

        siteName = nc.createVariable('siteName','c',('site','site_strlen',))
        siteName.long_name = 'Name of the spiral site'
        siteName.units = 'none'
        siteName[:] = SITE_NAMES[:]

        siteIndex = nc.createVariable('siteId','i4',('profile',),zlib=zlib)
        siteIndex.long_name = 'Index of the spiral site, 0 thru 7 (see siteName)'
        siteIndex.instance_dimension = 'site'
        siteIndex.units = 'none'
        siteIndex[:] = self.Site[:] - 1 # starts at 0 per CF conventions
            
        time = nc.createVariable('time','i4',('profile',))
        time.standard_name = 'time'
        time.long_name = 'Average UTC time of profile'
        time.units = 'seconds since 2011-06-01 00:00:00'
        time[:] = _secondsSince(self.Time, t_ref=datetime(2011,0o6,15))
            
        time_start = nc.createVariable('time_start','i4',('profile',))
        time_start.standard_name = 'time'
        time_start.long_name = 'Starting UTC time of profile'
        time_start.units = 'seconds since 2011-06-01 00:00:00'
        time_start[:] = _secondsSince(self.Time_start, t_ref=datetime(2011,0o6,15))

        time_stop = nc.createVariable('time_stop','i4',('profile',))
        time_stop.standard_name = 'time'
        time_stop.long_name = 'End UTC time of profile'
        time_stop.units = 'seconds since 2011-06-01 00:00:00'
        time_stop[:] = _secondsSince(self.Time_stop, t_ref=datetime(2011,0o6,15))

        year = nc.createVariable('year','i1',('profile',),zlib=zlib)
        year.long_name = 'Average UTC year of profile'
        year[:] = [ t.year for t in self.Time ]

        month = nc.createVariable('month','i1',('profile',),zlib=zlib)
        month.long_name = 'Avergage UTC month of profile'
        month[:] = [ t.month for t in self.Time ]

        day = nc.createVariable('day','i1',('profile',),zlib=zlib)
        day.long_name = 'Average UTC day of profile'
        day[:] = [ t.day for t in self.Time ]

        hour = nc.createVariable('hour','i1',('profile',),zlib=zlib)
        hour.long_name = 'Average UTC hour of profile'
        hour[:] = [ t.hour for t in self.Time ]

        minute = nc.createVariable('minute','i1',('profile',),zlib=zlib)
        minute.long_name = 'Average UTC minute of profile'
        minute[:] = [ t.minute for t in self.Time ]

        second = nc.createVariable('second','i1',('profile',),zlib=zlib)
        second.long_name = 'Average UTC second of profile'
        second[:] = [ t.second for t in self.Time ]

        flight = nc.createVariable('flight','i4',('profile',),zlib=zlib)
        flight.standard_name = 'flight'
        flight.long_name = 'Flight Id'
        flight.units = 'none'
        flight[:] = self.Flight[:]

        # Create data variables
        # ---------------------
        for var in sort(list(self.Variables.keys())):
            if Verbose:
                print("[] Writing %s"%var)
            v = nc.createVariable(var,'f4',('profile','altitude',),fill_value=NaN,zlib=zlib)
            v.standard_name = self.Std_name[var]
            v.long_name = self.Long_name[var]
            v.units = self.Units[var]
            v.coordinates = "time longitude latitude altitude"
            v.missing_value = NaN
            v[:,:] = self.Variables[var][:,:]

        # Close the file
        # --------------
        nc.close()

#........................................................................

def _secondsSince(tyme, t_ref):
    """
    Return array of with secodns sicne t_ref.
    """
    dt = tyme - t_ref
    for n in range(size(dt)):
        dt[n] = dt[n].total_seconds()

    return dt
    
#........................................................................
if __name__ == "__main__":

    p3 = P3_Profiles()
    p3.writeNC('DISCOVER-AQ_vertical_profiles_R0.nc',format='NETCDF3_CLASSIC')
    p3.writeNC('DISCOVER-AQ_vertical_profiles_R0.h5',format='NETCDF4',zlib=True)

    
        
