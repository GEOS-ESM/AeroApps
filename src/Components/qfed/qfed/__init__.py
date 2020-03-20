"""
         QFED - Quick Fire Emission Dataset

This module contains QFED support modules.

"""

from PlumeRise import PLUME_L2,PLUME_L3

class QFED(PLUME_L2):

    def print_l2b(self,filename=None):
        """
        Writes ASCII file with lon/lat/FRP/biome type and % fire area.
        """

        try:
            self.veg
        except:
            raise ValueError, '"veg" attribute not defined; use methods GetSimpleVeg() or GetDetailedVeg() to define vegetation type.'

        try:
            self.p
        except:
            raise ValueError, 'fire size attribute not defined; use method dozier() to define it.'

        if filename is None:
            f = sys.stdout
        else:
            f = open(filename,"w")

        N = self.lon.size
        jday = int(self.yyyy[N/2]*1000 + self.jjj[N/2])
        print >>f, '<fire_locations version=2.00 satellite="%s" total_number=%d jdate=%d time=#1 lon=#2 lat=#3 pow=#4 veg=#5 fire_size=#6 md5sum=n/a>'%(self.sat,N,jday)
        for i in range(N):
            print >>f, '%02d%02d %10.5f %10.5f %10.5f %2d %10.5f'%(\
                  self.hh[i],self.nn[i],self.lon[i],self.lat[i],\
                  self.pow[i],self.veg[i],self.p[i])

#.........................................................................

    def print_l2c(self,filename=None):
        """
        Writes ASCII file with lon/lat/FRP/biome type and % fire area.
        """

        try:
            self.veg
        except:
            raise ValueError, '"veg" attribute not defined; use methods GetSimpleVeg() or GetDetailedVeg() to define vegetation type.'

        try:
            self.p
        except:
            raise ValueError, 'fire size attribute not defined; use method dozier() to define it.'

        if filename is None:
            f = sys.stdout
        else:
            f = open(filename,"w")

        N = self.lon.size
        jday = int(self.yyyy[N/2]*1000 + self.jjj[N/2])
        print >>f, '<fire_locations version=2.00 satellite="%s" total_number=%d jdate=%d time=#1 lon=#2 lat=#3 pow=#4 veg=#5 r_F=#6 plume_k1=#7 plume_k2=#15 md5sum=n/a>'%(self.sat,N,jday)
        for i in range(N):
            pp = ( self.hh[i],self.nn[i],self.lon[i],self.lat[i],
                   self.pow[i],self.veg[i],self.r_F[i] )
            k1 = tuple(self.k_plume[:,i,0].astype(int))
            k2 = tuple(self.k_plume[:,i,1].astype(int))
            tt = pp + k1 + k2
            print >>f, \
            '%02d%02d %10.5f %10.5f %10.5f %2d %4.2f %2d %2d %2d %2d %2d %2d %2d %2d %2d %2d %2d %2d %2d %2d %2d %2d'%tt
             
             

#.........................................................................

    def print_l2d(self,filename=None):
        """
        Writes ASCII file with lon/lat/FRP/biome type and % fire area.
        """

        try:
            self.veg
        except:
            raise ValueError, '"veg" attribute not defined; use methods GetSimpleVeg() or GetDetailedVeg() to define vegetation type.'

        try:
            self.p
        except:
            raise ValueError, 'fire size attribute not defined; use method dozier() to define it.'

        if filename is None:
            f = sys.stdout
        else:
            f = open(filename,"w")

        N = self.lon.size
        jday = int(self.yyyy[N/2]*1000 + self.jjj[N/2])
        print >>f, '<fire_locations version=2.00 satellite="%s" total_number=%d jdate=%d time=#1 lon=#2 lat=#3 pow=#4 veg=#5 fire_size=#6 plume_p1=#7 plume_p2=#8 plume_k1=#9 plume_k2#10 md5sum=n/a>'%(self.sat,N,jday)
        for i in range(N):
            print >>f, \
                  '%02d%02d %10.5f %10.5f %10.5f %2d %10.5f %10.3f %10.3f %3d %3d'%(\
                  self.hh[i],self.nn[i],self.lon[i],self.lat[i],\
                  self.pow[i],self.veg[i],self.p[i],\
                  self.p_plume[i][0],self.p_plume[i][1],
                  self.k_plume[i][0],self.k_plume[i][1])

#.........................................................................


