"""

Implements base class to Gap Field albedo datasets used by the MODIS clouds group.

"""

from datetime import datetime, timedelta
from pyhdf.SD import SD

from glob import glob

from numpy import zeros, ones, array, int

class MCD43GF(object):

    def AlbedoOpenFile(self,time,root='/nobackup/10/MODIS/005/Level3/Albedo/data',ymin=2001,ymax=2012):
        """
        Returns albedo file name.
        """
        y = min(max(time.year,ymin),ymax)
        t = datetime(y,time.month,time.day,time.hour,time.minute,time.second)
        doy = t.toordinal() - datetime(y-1,12,31).toordinal()
        doy_ = 1 + 8 * int(0.5+(doy-1)/8.) # discrete doy
        filename = root+'/'+'%d/00-05.%03d/MCD43GF_wsa_Band4_%03d_%d_f.hdf'%(y,doy_,doy_,y)
        print "%03d %03d %s"%(doy,doy_,filename)
        return SD(filename)

    def AlbedoGetSDS(self,time,**kwopts):
        """
        Returns gridded albedo at time.
        """
        h = self.AlbedoOpenFile(time,**kwopts)
        return h.select('Albedo_Map_0.555')
 
    def AlbedoSample(self,lon,lat,time,Verbose=False,**kwopts):
        """
        Nearest-neighbor sampling of albedo file.
        """

        # Open first albedo file
        # ---------------------- 
        h = self.AlbedoOpenFile(time[0],**kwopts)

        # Scaling
        # -------
        att = h.select('Albedo_Map_0.555').attributes()
        fill, add, scale = att['_FillValue'], att['add_offset'], att['scale_factor']

        # coordinates
        # -----------
        Lon = h.select('Longitude')[:]
        Lat = - h.select('Latitude')[:] # really, negative of latitude
        dLon = Lon[1]-Lon[0]
        dLat = Lat[1]-Lat[0] 
        
        # Get unique list of 8-day DOY (i.e., days with albedo file)
        # ----------------------------------------------------------
        n = len(time)
        Doy = array([time[i].toordinal() - datetime(time[i].year-1,12,31).toordinal() for i in range(n)])
        DOY = array([1 + 8 * int(0.5+(Doy[i]-1)/8.) for i in range(n)])
        dicDOY = dict()
        for doy in DOY:
            dicDOY[doy] = 1

        # Loop over dates on a single albedo file
        # ---------------------------------------
        a = - ones(len(time))
        for doy in sorted(dicDOY.keys()):

            # coordinates for this 8-day doy
            # ------------------------------
            N = (DOY==doy)
            lon_, lat_, time_ = lon[N], -lat[N], time[N] # gather, notice -lat

            # Handle for albedo
            # -----------------
            A =  self.AlbedoGetSDS(time_[0],**kwopts)

            # Read from file
            # --------------
            n_ = len(time_)
            a_ = - ones(n_)
            for k in range(n_):
                i = int(0.5+(lon_[k]-Lon[0])/dLon)
                j = int(0.5+(lat_[k]-Lat[0])/dLat)
                a_[k] = A[j,i] # read one at a time --- can be optimized if needed

                if Verbose:
                    if a_[k] >= 0.0 and a_[k] != fill:
                        print time_[k], "%8.3f %8.3f %6.3f  ...%8.3f%%"\
                        %(lon_[k],lat_[k],scale*a_[k],100.*k/float(n_))

                
            a[N] = a_[:] # scatter
       
        # Scale and replace fill values
        # -----------------------------
        a[a==fill] = -1
        K = (a>0)
        a[K] = scale * a[K] + add

        self.albedo = a
        
        return
        
#--------------------------------------------------------------------------------------------

if __name__ == "__main__":

    from mapss import ANET
    
    m = MCD43GF()

    tdir = '/nobackup/MAPSS'
    adir = tdir+'/AERONET/AerosolOpticalDepth/2008'
    Files = sorted(glob(adir+'/mapss.ARNT.2008.07.??.*.dat'))

    a = ANET(Files,Verbose=True)

    I = (a.tau550>0)

    A = m.AlbedoSample(a.lon[I],a.lat[I],a.time[I],Verbose=True)
    
    
def later():
    t0 = datetime(2008,1,1)
    dt = timedelta(days=1)
    for n in range(31):
        t = t0 + n * dt
        a = m.AlbedoGrid(t)
    
