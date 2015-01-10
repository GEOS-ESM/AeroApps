"""

Very minimalistic ODS reader.

"""

VARS = ['lon', 'lat', 'lev', 'time', 'obs', 'qcexcl']
EXTRA = ['kt', 'kx', 'ks', 'qchist', 'xm', 'xvec', 'omf', 'oma']


from netCDF4 import Dataset

from numpy import array
from datetime import timedelta
from dateutil.parser import parse as isoparser

class ODSreader(object):

    def __init__ (self, filename, Extra=None):
        """
        Read select attributes from an ODS file.
        """ 

        nc = Dataset(filename)

        Vars = VARS
        if Extra is not None:
             Vars += Extra

        self._getVars(nc,Vars)

        self.nobs = len(self.lon)

        # Create python time
        # ------------------
        tunits = nc.variables['time'].units
        onemin = timedelta(minutes=1)
        toff  = isoparser(tunits.split(' since ')[1].replace(' ','T'))
        t = array([toff + self.time[i] * onemin for i in range(self.nobs)])

        self.time = t # array of python datetime objects

    def _getVars(self,nc,Vars):
        """
        Read and flatten variables on file. 
        Assuems 1 synoptic time per file.
        """
        for name in Vars:
            v = nc.variables[name]
            q = v[:,:]
            try:
                q = q * v.scale_factor
            except:
                pass
            try:
                q = q + v.add_offset
            except:
                pass
                
            I = q.mask==False  # not fill values
            self.__dict__[name] = q[I].data

#-------------

if __name__ == "__main__":

    m = ODSreader('nnr_002.MOD04_L2a.land.20140119_0300z.ods')

    for i in range(m.nobs):
        print "%6.2f %6.2f %s %7.4f"%(m.lon[i],m.lat[i],m.time[i].isoformat(),m.obs[i]) 
