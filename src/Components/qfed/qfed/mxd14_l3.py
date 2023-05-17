import os
import re

from types      import *
from glob       import glob
from datetime   import date, timedelta

from numpy      import (array, zeros, zeros_like, ones_like, any, 
                        linspace, arange, sum, sin, cos, sqrt, 
                        logical_and, logical_or, bitwise_and)

from pyhdf.SD   import *

from gfio       import GFIO
from pyobs      import IGBP_
from pyobs.binObs_ import binareas, binareasnr

__VERSION__ = 2.5
__CVSTAG__  = '@CVSTAG'

# MxD14 collection 6 fire mask pixel classes
WATER  = 3     # non-fire water pixel
CLOUD  = 4     # cloud (land or water)
NOFIRE = 5     # non-fire land pixel

# MxD14 collection 6 algorithm quality assessment bits: land/water state (bits 0-1)
QA_WATER = 0
QA_COAST = 1
QA_LAND  = 2


MODIS_TERRA = 'MODIS_TERRA'
MODIS_AQUA  = 'MODIS_AQUA'

class MxD14_L3(object):
    """
    This class implements the MODIS Level 2 fire products, usually
    referred to as MOD14 (TERRA satellite) and MYD14 (AQUA satellite).
    This version handles pixels with no fires.
    """

    def __init__ (self,Path,GeoDir,IgbpDir,refine=16,res=None,Verb=0):
       """
       Reads individual granules or a full day of Level 2 MOD14/MYD14 files
       present on a given *Path* and returns a single object with
       all data concatenated. On input, 

         Path -- can be a single file, a single directory, of a list
                 of files and directories.  Directories are
                 transversed recursively. If a non MOD14/MYD14 Level 2
                 file is encountered, it is simply ignored.
                 
         GeoDir -- directory where to find the geolocation MOD03/MYD03
                 files.

         IgbpDir -- directory for the IGBP database

         refine  -- refinement level for a base 4x5 GEOS-5 grid
                       refine=1  produces a      4x5       grid
                       refine=2  produces a      2x2.50    grid
                       refine=4  produces a      1x1.25    grid
                       refine=8  produces a   0.50x0.625   grid
                       refine=16 produces a   0.25x0.3125  grid
                       refine=32 produces a  0.125x0.15625 grid

       Alternatively, one can specify the grid resolution with a
       single letter:

         res     -- single letter denoting GEOS-5 resolution,
                       res='a'  produces a      4x5       grid
                       res='b'  produces a      2x2.50    grid
                       res='c'  produces a      1x1.25    grid
                       res='d'  produces a   0.50x0.625   grid
                       res='e'  produces a   0.25x0.3125  grid
                       res='f'  produces a  0.125x0.15625 grid
                    or,
                       res='0.1x0.1' produces   0.1x0.1   grid

                   NOTE: *res*, if specified, supersedes *refine*.

         Verb -- Verbose level:
                 0 - really quiet (default)
                 1 - Warns if invalid file is found
                 2 - Prints out non-zero number of fires in each file.

       """

#      Output grid resolution
#      ----------------------
       cubed = re.compile('c[0-9]+')
       if cubed.match(res):
          refine = None
          isCubed = True
       else:
          isCubed = False
          if res is not None:
              if res=='a': refine = 1 
              if res=='b': refine = 2
              if res=='c': refine = 4
              if res=='d': refine = 8
              if res=='e': refine = 16
              if res=='f': refine = 32
              if res=='0.1x0.1': refine = None

#      Initially are lists of numpy arrays for each granule
#      ----------------------------------------------------
       self.verb = Verb
       self.sat  = None # Satellite name
       self.GeoDir = GeoDir
       self.IgbpDir = IgbpDir

#      Place holder for date
#      ---------------------
       self.doy  = None
       self.date = None 
       self.col  = None # collection number
       
#      Lat lon grid
#      ------------
       self.grid_type = None

       if refine is not None:
           dx = 5. / refine
           dy = 4. / refine
           im = int(360. / dx)
           jm = int(180. / dy + 1)

           self.grid_type = 'GEOS-5 A-Grid' 

           self.im = im
           self.jm = jm
           self.glon = linspace(-180.,180.,self.im,endpoint=False)
           self.glat = linspace(-90.,90.,self.jm)
       else:
           if isCubed:
                 self.grid_type = 'GEOS-5 A-Grid'
                 self.im = int(res[1:len(res)])
                 self.jm = self.im*6
                 self.glon = linspace(1,self.im,self.im)
                 self.glat = linspace(1,self.jm,self.jm)
           else:
              if res=='0.1x0.1':
                 self.grid_type = 'DE_x_PE'

                 self.im = 3600
                 self.jm = 1800

                 d_lon = 360.0 / self.im
                 d_lat = 180.0 / self.jm

                 self.glon = linspace(-180+d_lon/2, 180-d_lon/2, self.im)
                 self.glat = linspace( -90+d_lat/2,  90-d_lat/2, self.jm)

       assert (self.grid_type != None)

#      Gridded accumulators
#      --------------------
       self.gLand  = zeros((self.im,self.jm))
       self.gWater = zeros((self.im,self.jm))
       self.gCloud = zeros((self.im,self.jm))
       self.gFRP   = zeros((4,self.im,self.jm))

#      Land/water state
#      ----------------
       self.lws = None

       if type(Path) is not ListType:
           Path = [Path, ]

#      Fire Products
#      -------------
       self._readList(Path)


#---
    def write(self,filename=None,dir={'ana':'.', 'bkg':'.'},expid='qfed2',tag=None, col=None, qc=True, bootstrap=False):
       """
       Writes gridded Areas and FRP to file.
       """

       meta = dict(vname  = ('land', 'water', 'cloud', 'frp_tf', 'frp_xf', 'frp_sv', 'frp_gl', 
                             'fb_tf', 'fb_xf', 'fb_sv', 'fb_gl'),
                   vtitle = ('Observed Clear Land Area',
                             'Water Area',
                             'Obscured by Clouds Area',
                             'Fire Radiative Power (Tropical Forests)',
                             'Fire Radiative Power (Extra-tropical Forests)',
                             'Fire Radiative Power (Savanna)',
                             'Fire Radiative Power (Grasslands)',
                             'Background FRP Density (Tropical Forests)',
                             'Background FRP Density (Extra-tropical Forests)',
                             'Background FRP Density (Savanna)',
                             'Background FRP Density (Grasslands)'),
                  vunits  = ('km2', 'km2', 'km2', 'MW', 'MW', 'MW', 'MW', 
                             'MW km-2', 'MW km-2', 'MW km-2', 'MW km-2'),
                  title   = 'QFED Level3a v%3.1f (%s) Gridded FRP Estimates'%(__VERSION__, _getTagName(tag)),
                  source  = 'NASA/GSFC/GMAO GEOS Aerosol Group',
                  contact = ('%s; %s') % ('arlindo.dasilva@nasa.gov', 'anton.darmenov@nasa.gov'))


       if self.date is None:
           print("[x] did not find matching files, skipped writing an output file")
           return

       if qc == True:
           self.qualityControl()
       else:
           print('[!] skipping QC procedures')

       self._write_ana(filename=filename,dir=dir['ana'],expid=expid,file_meta=meta,bootstrap=bootstrap)
    

    def _write_ana(self,filename=None,dir='.',expid='qfed2',file_meta=None,bootstrap=False):
       """
       Writes gridded Areas and FRP to file.
       """

       nymd = 10000*self.date.year + 100*self.date.month + self.date.day
       nhms = 120000

       if filename is None:
           self.filename = os.path.join(dir, '%s.frp.%s.%d.nc4'%(expid, self.col, nymd))
       else:
           self.filename = filename


       if bootstrap:
           print() 
           print('[i] Bootstrapping FRP forecast!')
           print()

           f = GFIO()

           f.create(self.filename, file_meta['vname'], nymd, nhms, 
                    lon=self.glon, lat=self.glat,
                    vtitle=file_meta['vtitle'], vunits=file_meta['vunits'], 
                    title=file_meta['title'], source=file_meta['source'], 
                    contact=file_meta['contact'])
       else:
           f = GFIO(self.filename, 'w')

       f.write(file_meta['vname'][0], nymd, nhms, self.gLand)
       f.write(file_meta['vname'][1], nymd, nhms, self.gWater)
       f.write(file_meta['vname'][2], nymd, nhms, self.gCloud)
       f.write(file_meta['vname'][3], nymd, nhms, self.gFRP[0,:,:])
       f.write(file_meta['vname'][4], nymd, nhms, self.gFRP[1,:,:])
       f.write(file_meta['vname'][5], nymd, nhms, self.gFRP[2,:,:])
       f.write(file_meta['vname'][6], nymd, nhms, self.gFRP[3,:,:])

       if bootstrap:
           f.write(file_meta['vname'][ 7], nymd, nhms, zeros_like(self.gFRP[0,:,:]))
           f.write(file_meta['vname'][ 8], nymd, nhms, zeros_like(self.gFRP[0,:,:]))
           f.write(file_meta['vname'][ 9], nymd, nhms, zeros_like(self.gFRP[0,:,:]))
           f.write(file_meta['vname'][10], nymd, nhms, zeros_like(self.gFRP[0,:,:]))

       try:
           f.close()
       except:
           pass

       if self.verb >=1:
           print("[w] Wrote file %s" % self.filename)


    def _write_bkg(self,filename=None,dir='.',expid='qfed2',file_meta=None, FillValue=1.0e20):
       """
       Creates L3a file valid at day + 1. Data fields are filled with FillValue.
       """

       _date = self.date + timedelta(days=1)

       nymd = 10000*_date.year + 100*_date.month + _date.day
       nhms = 120000

       if filename is None:
           _filename = os.path.join(dir, '%s.frp.%s.%d.nc4'%(expid, self.col, nymd))
       else:
          _dir, _file = os.path.split(filename)
          _filename = os.path.join(_dir, '_%s'%_file)

       f = GFIO()

       f.create(_filename, file_meta['vname'], nymd, nhms, 
                lon=self.glon, lat=self.glat,
                vtitle=file_meta['vtitle'], vunits=file_meta['vunits'], 
                title=file_meta['title'], source=file_meta['source'], 
                contact=file_meta['contact'])

       missing = zeros_like(self.gFRP[0,:,:])
       missing[:,:] = FillValue

       for v in file_meta['vname']:
           f.write(v, nymd, nhms, missing)
 
       try:
           f.close()
       except:
           pass

       if self.verb >=1:
           print("[w] Wrote file %s" % _filename)


#     ......................................................................

#---
    def _readList(self,List):
        """
        Recursively, look for files in list; list items can
        be files or directories.
        """
        for item in List:
            if os.path.isdir(item):      self._readDir(item)
            elif os.path.isfile(item):   self._readGranule(item)
            else:
                print("%s is not a valid file or directory, ignoring it"%item)
#---
    def _readDir(self,dir):
        """Recursively, look for files in directory."""
        for item in os.listdir(dir):
            path = dir + os.sep + item
            if os.path.isdir(path):      self._readDir(path)
            elif os.path.isfile(path):   self._readGranule(path)
            else:
                print("%s is not a valid file or directory, ignoring it"%item)

#---
    def _readGranule(self,filename):
        """Reads one MOD14/MYD14 granule with Level 2 fire data."""

#       Don't fuss if the file cannot be opened
#       ---------------------------------------
        try:
           mxd14 = SD(filename)
        except HDF4Error:
            if self.verb >= 1:
                print("- %s: not recognized as an HDF file"%filename)
            return 

#       Figure out MxD03 pathname
#       -------------------------
        str = mxd14.attributes()['MOD03 input file']
        try:
            i = str.index('MOD03.A')
            MxD03 = 'MOD03'
            self.sat = MODIS_TERRA
        except:
            i = str.index('MYD03.A')
            MxD03 = 'MYD03'
            self.sat = MODIS_AQUA

        MxD03 = str[i:i+5]
        base  = str[i:i+23]      # e.g., MOD03.A2003001.0000.005
        year  = str[i+7:i+11]    # e.g., 2003
        doy   = str[i+11:i+14]   # e.g., 001

        try:
            gfilename = glob(self.GeoDir+'/'+MxD03+'/'+year+'/'+doy+'/'+base+'*.hdf')[0]
        except:
            print("[x] cannot find geo-location file for <%s>, ignoring granule"%base)
            return 

#       Record date and collection
#       --------------------------
        if self.doy is None:
            self.doy  = int(doy)
            self.date = date(int(year),1,1) + timedelta(days=self.doy - 1)

        if self.col is None:
            # adopt the collection version of the MxD14 files 
            str = os.path.basename(filename)
            self.col = str.split('.')[3]
            
#       Read FRP retrieval and bin it
#       -----------------------------
        if mxd14.select('FP_longitude').checkempty():
            if self.verb >= 2:
                print("[x] no fires in granule <%s>, ignoring it <"%base)
            n_fires = 0

#       There are fires in granule
#       --------------------------
        else:

#           Open geolocation file
#           ---------------------
            try:
                mxd03 = SD(gfilename)
            except HDF4Error:
                if self.verb >= 1:
                    print("[x] cannot open geo-location file <%s>, ignoring granule"%gfilename)
                return 

#           Get no-fire areas
#           -----------------
            rc = self._readAreas(mxd14,mxd03)
            if rc:
                print("[x] problems with geo-location file <%s>, ignoring granule"%gfilename)
                return # inconsistent MOD03, skip granule

#           Read fire properties
#           --------------------
            lon = mxd14.select('FP_longitude').get()
            lat = mxd14.select('FP_latitude').get()
            frp = mxd14.select('FP_power').get()

            fp_line   = mxd14.select('FP_line').get()
            fp_sample = mxd14.select('FP_sample').get()

            area = _pixar(fp_sample)


#           Determine if there are fires from water pixels (likely offshore gas flaring) and exclude them
#           ---------------------------------------------------------------------------------------------
            n_fires_initial = frp.size

            i = [n for n in range(n_fires_initial) if self.lws[fp_line[n],fp_sample[n]] == QA_WATER]
            if len(i) > 0:
                if self.verb > 1:
                    print("      --> found %d FIRE pixel(s) over water" % len(i))

                self.gWater += _binareas(lon[i],lat[i],area[i],self.im,self.jm,grid_type=self.grid_type)

            i = [n for n in range(n_fires_initial) if self.lws[fp_line[n],fp_sample[n]] in (QA_COAST,QA_LAND)]
            if len(i) > 0:
                lon = lon[i]
                lat = lat[i]
                frp = frp[i]
                fp_line = fp_line[i]
                fp_sample = fp_sample[i]
                area = area[i] 
            else:
                if self.verb > 1:
                    print("      --> no FIRE pixels over land/coast")

                return
 
            n_fires = frp.size

            if n_fires_initial != n_fires:
                if self.verb > 1:
                    print("      --> reduced the number of FIRE pixels from %d to %d" % (n_fires_initial, n_fires))


#           Bin area of burning pixels
#           --------------------------
            self.gLand += _binareas(lon,lat,area,self.im,self.jm,grid_type=self.grid_type)

#           Bin FRP for each biome
#           ----------------------
            veg = _getSimpleVeg(lon,lat,self.IgbpDir)
            for b in (1,2,3,4):
                i = (veg==b)
                if any(i):
                    blon = lon[i]
                    blat = lat[i]
                    bfrp = frp[i]
                    self.gFRP[b-1,:,:] += _binareas(blon,blat,bfrp,self.im,self.jm,grid_type=self.grid_type)

        if n_fires>0 and self.verb>=1:
            print('[ ] Processed <'+os.path.basename(filename)[0:19]+'> with %4d fires'%n_fires)

#---
    def _readAreas(self,mxd14,mxd03):
        """Process an already open MOD03 file for NO-FIRE areas."""
    
#       Read coordinates and fire mask
#       -----------------------------
        Lon = mxd03.select('Longitude').get()
        Lat = mxd03.select('Latitude').get()
        fmask = mxd14.select('fire mask').get()

        nLines,  nSamples = Lon.shape
        nLines2, nSamples2 = fmask.shape

        mismatch = (nLines!=nLines2) | (nSamples!=nSamples2)

#       Find indices corresponding to no fires
#       --------------------------------------
        if mismatch:
            # nLines = nLines2
            # Lon = [0:nLines,0:nSample]
            # Lat = [0:nLines,0:nSample]
            return 1
      
#       Find indices with coordinates within the valid ranges
#       -----------------------------------------------------
        try:
            minLon, maxLon = mxd03.select('Longitude').getrange()
        except:
            minLon, maxLon = (-180.0, 180.0)

        try:
            minLat, maxLat = mxd03.select('Latitude').getrange()
        except:
            minLat, maxLat = (-90.0, 90.0)

        i_lon = logical_and(Lon >= minLon, Lon <= maxLon)
        i_lat = logical_and(Lat >= minLat, Lat <= maxLat)
        valid = logical_and(i_lon, i_lat)


#       algorithm QA 
#       ------------
        qa = mxd14.select('algorithm QA').get()
        self.lws = bitwise_and(qa, 3)  # land/water state is stored in bits 0-1

        i_water = logical_and(self.lws==QA_WATER, valid)
        i_land  = logical_and(logical_or(self.lws==QA_COAST, self.lws==QA_LAND), valid)


#       Calculate pixel area
#       --------------------
        Area = zeros((nLines,nSamples))
        Area[:] = _pixar(1+arange(nSamples))

        # non-fire land pixel
        i = logical_and(fmask==NOFIRE, valid)

        if any(i):

#           Condensed 1D arrays of clear-land not burning pixels
#           ----------------------------------------------------
            lon = Lon[i].ravel()
            lat = Lat[i].ravel()
            area = Area[i].ravel()

#           Bin areas of no-fires and sum
#           -----------------------------
            self.gLand += _binareas(lon,lat,area,self.im,self.jm,grid_type=self.grid_type)

        else:
            if self.verb > 1:
                print("      --> no NOFIRE pixel for granule")

        # non-fire water or cloud over water
        i = logical_or(logical_and(fmask==WATER, valid), logical_and(logical_and(fmask==CLOUD, valid), i_water))

        if any(i):

#           Condensed 1D arrays of water pixels
#           -----------------------------------
            lon = Lon[i].ravel()
            lat = Lat[i].ravel()
            area = Area[i].ravel()

#           Bin areas of water and sum
#           --------------------------
            self.gWater += _binareas(lon,lat,area,self.im,self.jm,grid_type=self.grid_type)

        else:
            if self.verb > 1:
                print("      --> no WATER pixel for granule")

        # cloud over land only
        i = logical_and(logical_and(fmask==CLOUD, valid), i_land)

        if any(i):

#           Condensed 1D arrays of cloud pixels
#           -----------------------------------
            lon = Lon[i].ravel()
            lat = Lat[i].ravel()
            area = Area[i].ravel()

#           Bin areas of cloud and sum
#           --------------------------
            self.gCloud += _binareas(lon,lat,area,self.im,self.jm,grid_type=self.grid_type)

        else:
            if self.verb > 1:
                print("      --> no CLOUD pixel for granule")

        return 0
    
###        if self.verb>=1:
###            print '[2] Processed NO-FIRE AREAS for <'+os.path.basename(filename)[0:19]+'>'

#---
    def qualityControl(self, pom_oc_ratio=1.4, oc_mass_ext_coeff=4.0, max_aod=10.0):
        """
        Cap FRP to prevent large AOD values.

        Assumes:
               *) AOD > 10 should not be allowed
               *) AOD_OC = 0.9 AOD
               *) OC mass extinction coefficient = 4 m2 g-1
               *) OC to POM conversion factor = 1.4

        For more details see emissions.py.
        """     

        # Biome-dependent CO emission factors (Andreae & Merlet 2001)
        # Units: g(species) / kg(dry mater)
        # --------------------------------
        B_f = (5.20, 9.14, 3.40, 3.40)

        # Combustion rate constant
        Alpha = 1.37e-6 # kg(dry mater)/J
        A_f   = Alpha * array((2.5, 4.5,  1.8, 1.8))
       
        # Satellite Fudge Factor
        # ----------------------
        S_f = {}
        S_f[MODIS_TERRA] = 1.385
        S_f[MODIS_AQUA ] = 0.473

        # factor needed to convert B_f from [g/kg] to [kg/kg]
        units_factor = 1.0e-3

        n_biomes = 4

        E = zeros((n_biomes, self.im, self.jm))

        A_l = self.gLand
        A_w = self.gWater
        A_c = self.gCloud

        A_o = A_l + A_w

        # apply the 'sequential-b0' method to compute the emissions
        for b in range(n_biomes):
            E[b,:,:] += units_factor * A_f[b] * S_f[self.sat] * B_f[b] * self.gFRP[b,:,:]

        i = (A_l > 0)
        for b in range(n_biomes):
            E_b = E[b,:,:]
            E_b[i] = E_b[i] / (A_o[i] + A_c[i]) * ((A_l[i] + 2*A_c[i]) / (A_l[i] + A_c[i]))
            E[b,:,:] = E_b
        
        E_total = sum(E[:,:,:], axis=0)
        
        # column density of OC emitted for 24 hours, g m-2
        M = (1e3 * E_total) * (24 * 3600) 

        # cap FRP if the emissions are too strong
        aod_oc = oc_mass_ext_coeff * (pom_oc_ratio * M)

        f_phys = 6  # empirical factor - accounts for absence of removal processes
        max_aod_oc = f_phys * max_aod

        i_cap = aod_oc > max_aod_oc

        if (self.verb > 1) and any(i_cap):
            print('Maximum estimated AOD(OC) : %.2f' % aod_oc.max())
            print('Maximum FRP(TF)           : %.2f' % self.gFRP[0,:,:].max())
            print('Maximum FRP(XF)           : %.2f' % self.gFRP[1,:,:].max())
            print('Maximum FRP(SV)           : %.2f' % self.gFRP[2,:,:].max())
            print('Maximum FRP(GL)           : %.2f' % self.gFRP[3,:,:].max())

        q = ones_like(aod_oc)
        q[i_cap] = max_aod_oc / aod_oc[i_cap]

        for b in range(n_biomes):
            FRP = q*self.gFRP[b,:,:]
            self.gFRP[b,:,:] = FRP
        
        n_cap = sum(i_cap)
        if n_cap > 0:
            p_cap = 100.0 * n_cap / (sum(E_total.ravel() > 0))
            print('[!] FRPs in %d grid cells (%.1f%% of grid cells with fires) were capped.' % (n_cap, p_cap))

        if (self.verb > 1) and any(i_cap):
            print('After capping the emissions')
            print('Maximum FRP(TF)           : %.2f' % self.gFRP[0,:,:].max())
            print('Maximum FRP(XF)           : %.2f' % self.gFRP[1,:,:].max())
            print('Maximum FRP(SV)           : %.2f' % self.gFRP[2,:,:].max())
            print('Maximum FRP(GL)           : %.2f' % self.gFRP[3,:,:].max())

#..............................................................

def _binareas(lon, lat, area, im, jm, grid_type='GEOS-5 A-Grid'):

    if grid_type == 'GEOS-5 A-Grid':
        result = binareas(lon,lat,area,im,jm)
    elif grid_type == 'DE_x_PE':
        result = binareasnr(lon,lat,area,im,jm)
    else:
        result = None        
   
    return result


def _pixar(sample):
    """
    Compute pixel area given the sample number. Description is
    given in 'MODIS Collection 5 Active Fire Product User's Guide Version 2.4',
    Giglio, L., p. 44, 2010.
    """

    # parameters
    s  = 0.0014184397
    Re = 6378.137                # Earth radius, km
    h  = 705.0                   # satellite altitude, km

    sa = s * (sample - 676.5)    # scan angle, radians

    r  = Re + h
    q  = Re / r

    cos_sa   = cos(sa)
    sin_sa   = sin(sa)
    sqrt_trm = sqrt(q*q - sin_sa*sin_sa)
    
    dS = s*Re * (cos_sa / sqrt_trm - 1)
    dT = s*r  * (cos_sa - sqrt_trm)

    area = dS * dT

    return area


def _getSimpleVeg(lon,lat,Path,nonVeg=None):
    """
        Given the pathname for the IGBP datasets, returns
    information about the *aggregated* vegetation type:  

         1  Tropical Forest         IGBP  1, 30S < lat < 30N
         2  Extra-tropical Forests  IGBP  1, 2(lat <=30S or lat >= 30N),
                                          3, 4, 5
         3  cerrado/woody savanna   IGBP  6 thru  9
         4  Grassland/cropland      IGBP 10 thru 17

     the new attribute name is "veg". Notice that this module also
     defines the following constants:

         TROPICAL = 1
         EXTRA_TROPICAL = 2
         SAVANNA = 3
         GRASSLAND = 4

     corresponding to each aggregated vegetation type.
    """

    nobs = lon.shape[0]
    veg = IGBP_.getsimpleveg(lon,lat,Path+'/IGBP',nobs) # Fortran

    # substitute non vegetation (water, snow/ice) data with 
    # another type, e.g. GRASSLAND by default
    if nonVeg != None:
        i = logical_or(veg == -15, veg == -17)
        veg[i] = nonVeg # could be one of the biomes, i.e., 1, 2, 3 or 4

    return veg

def _getTagName(tag):
    if tag != None:
        tag_name = tag
    else:    
        if __CVSTAG__ not in (None, ''):
            tag_name = __CVSTAG__
        else:
            tag_name = 'unknown'

    return tag_name
    

# .............................................................

if __name__ == "__main__":

    Path = ['/nobackup/2/MODIS/Level2/MOD14/2003/001',
            '/nobackup/2/MODIS/Level2/MYD14/2003/001',
           ]
    GeoDir = '/nobackup/2/MODIS/Level1'
    IgbpDir = '/nobackup/Emissions/Vegetation/GL_IGBP_INPE'

    fires = MxD14_L3(Path,GeoDir,IgbpDir,res=1,Verb=1)
