"""
class to read in AMES BRDF data on tiles
"""

from netCDF4 import Dataset
import numpy as np
import xarray as xr
from glob import glob

# number of tiles
nH = 60
nV = 20
nlon = 600
nlat = 600


SDS =  ['Kg','Ki','Kv','mask']


# dictionary of tiles lon/lon extents
textent = {}
for H in range(nH):
    for V in range(nV):
        latmax = 60.0 - (V*6 + 0.0*0.01) - 0.005
        latmin = 60.0 - (V*6 + (nlat-1)*0.01) - 0.005

        lonmin = -180.0 + (H*6 + 0.0*0.01) + 0.005
        lonmax = -180.0 + (H*6 + (nlon-1)*0.01) + 0.005
        textent['H{}V{}'.format(H,V)] = [lonmin,lonmax,latmin,latmax]


class AMES_BRDF(object):
    """
    class for ames brdf data that is on tiles
    must provide either H & V that defines the tile, 
    or a numpy boolean array (grid) of shape nH,nV with the indeces of the tiles desired set to True
    desired tiles must be continuous

    return_ds = return an xarray dataset rather than individual numpy arrays. Right now only works is grid
    """

    def __init__(self,Path,doy,H=None,V=None,grid=None,SDS=SDS,return_ds=False):
        self.SDS = SDS
        self.return_ds = return_ds

        if H is not None:
            self.readTile(Path,doy,H,V)

        elif grid is not None:
            self.readGrid(Path,doy,grid)

        else:
            raise SyntaxError("Either grid or H & V must be not None")

    
    def readTile(self,Path,doy,H,V):
        """
        reads a single tile
        """
        fname = '{}/h{:02d}v{:02d}/sbg_mdl_h{:02d}v{:02d}_2019{}_brdf.nc4'.format(Path,H,V,H,V,doy)
        nc = Dataset(fname)
        
        # only read wavelength once
        if not hasattr(self, 'wavelength'):
            self.wavelength = nc.groups['ancillary'].variables['wavelength'][:]
            self.nwav       = len(self.wavelength)

        # read is sds        
        for sds in self.SDS:
            self.__dict__[sds] = nc.groups['gridded'].variables[sds][:]
        
        self.fname = fname

    def readGrid(self,Path,doy,grid):
        """
        reads several tiles defined in a grid
        """
        hextent,vextent = np.where(grid)
        flist = []
        DS = []
        for H,V in zip(hextent,vextent):
            fname = '{}/h{:02d}v{:02d}/sbg_mdl_h{:02d}v{:02d}_2019{}_brdf.nc4'.format(Path,H,V,H,V,doy)
            ds = xr.open_dataset(fname,group='gridded')
            tlons,tlats = self.getCoords(H,V)
            
            flist.append(fname)
            DS.append(ds.assign_coords({"lon":tlons,"lat":tlats}))

            if not hasattr(self, 'wavelength'):
                nc = Dataset(fname)
                self.wavelength = nc.groups['ancillary'].variables['wavelength'][:]
                self.nwav       = len(self.wavelength)
                nc.close()

        # stich together the tiles
        combined = xr.combine_by_coords(DS)

        if self.return_ds:
            self.ds = combined
        else:
            for sds in self.SDS:
                self.__dict__[sds] = combined.variables[sds].values

        self.flist = flist
    def getCoords(self,H,V):
        tlats = 60.0 - (V*6 + np.arange(nlat)*0.01) - 0.005
        tlons = -180.0 + (H*6 + np.arange(nlon)*0.01) + 0.005    

        return tlons,tlats

def tiles(box):
    """
    uses shapely polygons to find if a tile intersects with a given bounding box
    you can use shapely.plotting.plot_polygon to plot the polygons
    """
    from shapely.geometry import Polygon
    from shapely import intersects

    blonmin,blonmax,blatmin,blatmax = box
    coords = ((blonmin,blatmin),(blonmin,blatmax),(blonmax,blatmax),(blonmax,blatmin))
    boxpoly = Polygon(coords)

    grid = np.zeros([nH,nV]).astype(bool)
    for H in range(nH):
        for V in range(nV):
            lonmin,lonmax,latmin,latmax = textent['H{}V{}'.format(H,V)]
            coords = ((lonmin,latmin),(lonmin,latmax),(lonmax,latmax),(lonmax,latmin))
            tilepoly = Polygon(coords)
            grid[H,V] = intersects(boxpoly,tilepoly)
               

    return grid 


        



