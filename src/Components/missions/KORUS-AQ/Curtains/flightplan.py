"""
Simple class for parsing ASCII flight path files.
"""

from datetime import datetime, timedelta
from numpy    import array, pi, cos, sin, zeros, arccos, sqrt, linspace, \
                     arcsin, arctan, NaN
import dateutil.parser
from types import *
import sys

class FlightPlan(object):

    def __init__(self,filename,ds_km=10):
        """
        Load flight plan file.
        """
        Lines = open(filename).readlines()

        self.header = Lines[0].replace('\n','').replace('\r','')
        self.aircraft = self.header.split()[0]
        
        self.hrs = []
        self.tyme = []
        self.lon = []
        self.lat = []
        self.alt = []
        gotit = False
    
        for line in Lines:
            tokens = line.split()
            if len(tokens) == 0:
                continue
            elif 'Takeoff' in tokens[0]: 
                self.Takeoff = tokens[2]
                self.takeoff = dateutil.parser.parse(self.Takeoff)
                print 'Takeoff: ',self.takeoff
            elif 'Elapsed' in tokens[0]:
                gotit = True
            elif '(hrs)' in tokens[0]:
                continue
            elif gotit:
                alt_factor = 0.3048
                self.hrs += [float(tokens[0]),]
                self.tyme += [self.takeoff + timedelta(seconds=float(tokens[0])*60*60),]
                self.lat += [float(tokens[1]),]
                self.lon += [float(tokens[2]),]
                self.alt += [alt_factor*float(tokens[3]),] # in meters

        self.hrs = array(self.hrs)
        self.lon = array(self.lon)
        self.lat = array(self.lat)
        self.alt = array(self.alt)
        self.tyme = array(self.tyme)

        # Hour of the day (UTC)
        # ---------------------
        t0 = self.takeoff
        tbeg = datetime(t0.year,t0.month,t0.day)
        self.hour = array([(t-tbeg).total_seconds()/3600. for t in self.tyme])

        # Refine trajectory if so desired
        # -------------------------------
        self.refine(ds_km=ds_km) # refine trajectory
        
    def refine(self,ds_km=10.,dt_min=1.):
        """
        Refine the trajectory so that each point is about ds_km apart.
        """
        ## Interpolate in distance
        self.dst = _getDist(self.lon,self.lat)        
        ds = ds_km * 1000.#1000 m/km
        ## Interpolate in time
        dt = zeros(self.lon.shape)
        dt[1:] = [self.hour[i+1]-self.hour[i] for i in range(len(self.hour)-1)]
        self.dt = dt # time (in hours) between points
        
        X, Y, Z, T = [],[],[], [] # Distance
        X2, Y2, Z2, T2 = [],[],[], [] # Time
        for i in range(1,len(self.lon)):
            n = int(0.5+self.dst[i]/ds)
            X += list(linspace(self.lon[i-1],self.lon[i],n,endpoint=False))
            Y += list(linspace(self.lat[i-1],self.lat[i],n,endpoint=False))
            Z += list(linspace(self.alt[i-1],self.alt[i],n,endpoint=False))
            T += list(linspace(self.hour[i-1],self.hour[i],n,endpoint=False))

            n2 = int(self.dt[i]*60.)
            X2 += list(linspace(self.lon[i-1],self.lon[i],n2,endpoint=False))
            Y2 += list(linspace(self.lat[i-1],self.lat[i],n2,endpoint=False))
            Z2 += list(linspace(self.alt[i-1],self.alt[i],n2,endpoint=False))
            T2 += list(linspace(self.hour[i-1],self.hour[i],n2,endpoint=False))
        X += [self.lon[-1],]
        Y += [self.lat[-1],]
        Z += [self.alt[-1],]
        T += [self.hour[-1],]
        X2 += [self.lon[-1],]
        Y2 += [self.lat[-1],]
        Z2 += [self.alt[-1],]
        T2 += [self.hour[-1],]

        self.Longitude,self.Latitude,self.Altitude,self.Hour = array(X),array(Y),array(Z),array(T)
        self.Longitude,self.Latitude,self.Altitude,self.Hour = array(X2),array(Y2),array(Z2),array(T2)

        today = datetime(self.takeoff.year,self.takeoff.month,self.takeoff.day)
        self.Tyme = array([today+timedelta(seconds=h*3600) for h in T])
        self.Tyme = array([today+timedelta(seconds=h*3600) for h in T2])

def _gcDist(x1,y1,z1,x2,y2,z2):
    """
    Return great circle distance.
    """
    a = (6378137.00+6356752.3142)/2. # mean earth radius
    cosa = x1*x2 + y1*y2 + z1*z2
    cosa[cosa>1.] = 1.
    return a * arccos(cosa)
    
def _chDist(x1,y1,z1,x2,y2,z2):
    """
    Return chordal distance.
    """
    a = (6378137.00+6356752.3142)/2. # mean earth radius
    return a * sqrt((x1-x2)**2 + (y1-y2)**2 + (z1-z2)**2)
    
def _getDist(lon,lat):
    """
    Return distance along trajectory.
    """
    d2r = pi / 180.
    rlon, rlat = (d2r * lon, d2r * lat)
    x, y, z = (cos(rlat)*cos(rlon),cos(rlat)*sin(rlon),sin(rlat))
    dist = zeros(lon.shape)
    dist[1:] = _gcDist(x[:-1],y[:-1],z[:-1],
                       x[1:], y[1:], z[1:])
    return dist

#--
if __name__ == "__main__":

    from grads import GrADS
    
    asm_Np = 'http://opendap.nccs.nasa.gov:9090/dods/GEOS-5/fp/0.25_deg/fcast/inst3_3d_asm_Np.latest'
    ext_Np = 'http://opendap.nccs.nasa.gov:9090/dods/GEOS-5/fp/0.25_deg/fcast/inst3_3d_ext_Np.latest'
    fname = 'SEAC4RS_130816_s3_er2_1_table_path.txt'

    asm_Np = 'http://opendap.nccs.nasa.gov:9090/dods/GEOS-5/fp/0.25_deg/fcast/inst3_3d_asm_Np/inst3_3d_asm_Np.20130814_12'
    ext_Np = 'http://opendap.nccs.nasa.gov:9090/dods/GEOS-5/fp/0.25_deg/fcast/inst3_3d_ext_Np/inst3_3d_ext_Np.20130814_12'
    
    f = FlightPlan(fname)
        
    ga = GrADS(Window=False,Echo=False)
    fh = ga.open(asm_Np)
    ga('set lev 1000 300')
    
    f.addVar(ga,('h','T'))
    
    
