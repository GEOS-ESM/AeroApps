"""
Simple classes to read merge and other NetCDF from ORACLES.

"""

from datetime import datetime, timedelta
from dateutil.parser import parse as isoparser

from numpy    import ones, NaN, concatenate, array, pi, cos, sin, arccos, zeros
from glob     import glob

from netCDF4  import Dataset

MISSING = -99999.0
ULOD    = -77777.0 # Upper Limit of Detection (LOD) flag
LLOD    = -88888.0 # Lower Limit of Detection (LOD) flag

class ORACLES(object):
    """Reads ORACLES merge  files into Numpy arrays"""

    def __init__ (self,Filename,Vars=[],Verbose=False):
        """
        Opends a NetCDF merge file, reads in coordinates and 
        """

        nc = Dataset(Filename)

        # time
        # ----
        time = nc.variables['time'][:]
        tunits = nc.variables['time'].units
        onesec = timedelta(seconds=1)
        toff  = isoparser(tunits.split(' since ')[1].replace(' ','T'))
        self.tyme = array([toff + int(time[i]) * onesec for i in range(len(time))])

        # Spatial coordinates
        # -------------------
        self.lon = nc.variables['Longitude'][:]
        self.lat = nc.variables['Latitude'][:]
        self.alt = nc.variables['GPS_Altitude'][:]
        self.prs = nc.variables['Static_Pressure'][:]

        # for the other variables, just get a netCDF variable reference
        # -------------------------------------------------------------
        for name in nc.variables:
            if name in Vars:
                if Verbose:
                    print('[] reading ', name)
                self.__dict__[name] = nc.variables[name][:]
            else:
                if Verbose:
                    print('<> referencing ', name)
                self.__dict__[name] = nc.variables[name]

#--------------------------------------------------------------------------------------

if __name__ == "__main__":

    x = ORACLES('/home/adasilva/iesa/aerosol/campaigns/ORACLES/P3/MRG1/mrg1_P3_20160906_R32.nc',Verbose=True)

    
        
