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
#import xml.sax.handler

from numpy    import array, ones, zeros, concatenate, pi, sin, cos, arccos, cumsum
from datetime import datetime, timedelta

class SimpleKML(object):
    """
    Class to read KML files exported from Google Maps with directions.
    """

    def __init__(self,kml_filename):
        """
        Parses a KML object with Google Maps driving directions.
        """
        kml = _xml2obj(open(kml_filename))
        self.Document = kml.Document

    def getCoords(self,speed=90.,t0=None,metric=True,noSinglePoint=True):
        """
        Return (lon,lat) coordinates of all placemarks, along with
        distance and time vectors:

              lon, lat, dst, time = self.getCoords()

        where speed is the average speed, in km/h if *metric* is
        True, or miles/h otherwise. The initial time *t0* (a datetime
        object in UTC, not local time) defaults to now if not specified.

        On output, the distance *dst* is given in either km or miles,
        depending on *metric*.
        
        By default, single point placemarks are ignored.
        """

        # Lon/lat coordinates
        # -------------------
        lon, lat = ([], [])
        for n in range(len(self.Document.Placemark)):
            lon_, lat_ = self.getLonLat(n)
            if len(lon_)==1 and noSinglePoint:
                continue 
            lon.append(lon_)
            lat.append(lat_)
        lon = concatenate(lon)
        lat = concatenate(lat)

        # distance
        # --------
        dst = _getDist(lon,lat)             # meters
        if metric:
            dst = dst / 1000.               # km
        else:
            dst = (0.621371192/1000.) * dst # miles

        dst = cumsum(dst)

        # Elapsed time (in secs)
        # ----------------------
        DT = 3600. * dst / speed

        # Estimated time
        # --------------
        if t0 is None:
            t0 = datetime(1,1,1).utcnow()

        time = array([t0+timedelta(seconds=int(dt)) for dt in DT])
            
        return (lon,lat,dst,time)

    def nPlacemarks(self):
        """
        Return number of placemarks.
        """
        return len(self.Document.Placemark)
    
    def getLonLat(self,n):
        """
        Return coordinates of Placemark number "n"
        """
        if self.Document.Placemark[n].LineString is not None:
            coords = array([ array(s.split(',')).astype('float') for s in \
                             self.Document.Placemark[n].LineString.coordinates.split('\n')])
        elif self.Document.Placemark[n].Point is not None:
            coords = ones((1,3))
            coords[0,:] = array(self.Document.Placemark[n].Point.coordinates.split(',')).astype('float')
        else:
            raise ValueError, 'Placemak does not have Point nor LineString'
        
        return (coords[:,0], coords[:,1])

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

    class TreeBuilder(): #xml.sax.handler.ContentHandler):
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

    kml = SimpleKML("Directions_to_Acadia.kml")

    print "KML Document Name: ", kml.Document.name
    
    lon, lat, dst, time = kml.getCoords()
    
