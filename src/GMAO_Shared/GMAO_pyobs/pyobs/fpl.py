#--------------------------------------------------------------------------
#
#    Copyright (C) 2006-2011 by Arlindo da Silva <dasilva@opengrads.org>
#    All Rights Reserved.
#
#    This program is free software; you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation# using version 2 of the License.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY# without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with this program# if not, please consult  
#              
#              http://www.gnu.org/licenses/licenses.html
#
#    or write to the Free Software Foundation, Inc., 59 Temple Place,
#    Suite 330, Boston, MA 02111-1307 USA
#
#------------------------------------------------------------------------

"""
This module implements a simple class for parsing KML files exported
from Google Maps with driving directions.
"""

__version__ = '1.0.0'

import re
import xml.sax.handler

from numpy    import array, ones, zeros, concatenate, pi, sin, cos, arccos, cumsum, interp, savez
from datetime import datetime, timedelta

class FPL(object):
    """
    Class to read KML files exported from ForeFlights.
    """

    def __init__(self,kml_filename):
        """
        Parses a KML object with Google Maps driving directions.
        """
        kml = _xml2obj(open(kml_filename))
        self.kml= kml

        self.aircraft = kml['aircraft']['aircraft_tailnumber']

        # Lon/lat coordinates
        # -------------------
        self.waypoints = kml['waypoint_table']['waypoint']
        self.n = len(self.waypoints)
        self.identifier, self.lon, self.lat = [], [], []
        for wp in self.waypoints:
            print wp
            self.identifier += [wp['identifier'],]
            self.lon += [float(wp['lon']),]
            self.lat += [float(wp['lat']),]
        self.lon = array(self.lon)
        self.lat = array(self.lat)
            
    def getCoords(self,speed_kts=250,dt_min=1,t0=None,npzFile=None):
        """
        Return refined (lon,lat) coordinates of all waypoints, along with
        time vectors:

              lon, lat, time = self.getCoords()

        where speed_kts is the average speed in knots.  The initial time *t0* (a datetime
        object in UTC, not local time) defaults to now if not specified.

        By default, single point placemarks are ignored.
        """

        # distance
        # --------
        dst = _getDist(self.lon,self.lat)             # meters
        dst = dst / 1000.               # km
        dst = cumsum(dst)

        # Elapsed time (in secs)
        # ----------------------
        speed_kmh = speed_kts * 1.852
        dt_secs = 3600. * dst / speed_kmh

        # Estimated time
        # --------------
        if t0 is None:
            t0 = datetime(1,1,1).utcnow()

        t0 = datetime(t0.year,t0.month,t0.day,t0.hour,t0.minute,int(t0.second),0) # zero microsecond
            
        time = array([t0+timedelta(seconds=int(dt)) for dt in dt_secs])
        t_secs = array([(t-time[0]).total_seconds() for t in time ])
        
        # Next refine coordinates given dt
        # --------------------------------
        Mins = 1 + int( (time[-1]-time[0]).total_seconds()/60. )  # total number of minutes
        rTime_secs = array ( [ m*60 for m in range(Mins) ] ) # refined time in secs
        rTime_secs[-1] = t_secs[-1] # hack to make sure flight path closes
        rTime = array([t0+timedelta(seconds=int(dt)) for dt in rTime_secs])

        # Simple minded approach: interp lat/lon
        # Fix: great circle calculation
        # --------------------------------------
        rLon = interp(rTime_secs, t_secs, f.lon)
        rLat = interp(rTime_secs, t_secs, f.lat)
        print rTime_secs[-1], t_secs[-1]
        
        f.time = time

        # Write NPZ file
        # --------------
        if npzFile is not None:
            savez(npzFile,lon=rLon,lat=rLat,time=rTime)
        
        return (rLon, rLat, rTime)

def _gcDist(x1,y1,z1,x2,y2,z2):
    """
    Return great circle distance.
    """
    a = (6378137.00+6356752.3142)/2. # mean earth radius
    cosa = x1*x2 + y1*y2 + z1*z2
    cosa[cosa>1.] = 1.
    return a * arccos(cosa)
    
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

# .....................................................................
# The code below is from http://code.activestate.com/recipes/534109/
# and is distributed under the PSF License.
#
def _xml2obj(src):
    """
    A simple function to converts XML data into native Python object.
    """

    non_id_char = re.compile('[^_0-9a-zA-Z]')
    def _name_mangle(name):
        return non_id_char.sub('_', name)

    class DataNode(object):
        def __init__(self):
            self._attrs = {}    # XML attributes and child elements
            self.data = None    # child text data
        def __len__(self):
            # treat single element as a list of 1
            return 1
        def __getitem__(self, key):
            if isinstance(key, basestring):
                return self._attrs.get(key,None)
            else:
                return [self][key]
        def __contains__(self, name):
            return self._attrs.has_key(name)
        def __nonzero__(self):
            return bool(self._attrs or self.data)
        def __getattr__(self, name):
            if name.startswith('__'):
                # need to do this for Python special methods???
                raise AttributeError(name)
            return self._attrs.get(name,None)
        def _add_xml_attr(self, name, value):
            if name in self._attrs:
                # multiple attribute of the same name are represented by a list
                children = self._attrs[name]
                if not isinstance(children, list):
                    children = [children]
                    self._attrs[name] = children
                children.append(value)
            else:
                self._attrs[name] = value
        def __str__(self):
            return self.data or ''
        def __repr__(self):
            items = sorted(self._attrs.items())
            if self.data:
                items.append(('data', self.data))
            return u'{%s}' % ', '.join([u'%s:%s' % (k,repr(v)) for k,v in items])

    class TreeBuilder(xml.sax.handler.ContentHandler):
        def __init__(self):
            self.stack = []
            self.root = DataNode()
            self.current = self.root
            self.text_parts = []
        def startElement(self, name, attrs):
            self.stack.append((self.current, self.text_parts))
            self.current = DataNode()
            self.text_parts = []
            # xml attributes --> python attributes
            for k, v in attrs.items():
                self.current._add_xml_attr(_name_mangle(k), v)
        def endElement(self, name):
            text = ''.join(self.text_parts).strip()
            if text:
                self.current.data = text
            if self.current._attrs:
                obj = self.current
            else:
                # a text only node is simply represented by the string
                obj = text or ''
            self.current, self.text_parts = self.stack.pop()
            self.current._add_xml_attr(_name_mangle(name), obj)
        def characters(self, content):
            self.text_parts.append(content)

    builder = TreeBuilder()
    if isinstance(src,basestring):
        xml.sax.parseString(src, builder)
    else:
        xml.sax.parse(src, builder)
    return builder.root._attrs.values()[0]

#...........................................................................................

if __name__ == "__main__":

    f = FPL('wecan_2018-08-23.fpl')

    print 'Waypoints: ', f.identifier[:]
    print '      lon: ', f.lon[:]
    print '      lat: ', f.lat[:]
    
    lon,lat,time = f.getCoords(t0=datetime(2018,8,23,18,30,0),
                               npzFile='wecan.plan.2018-08-23.npz')

    
