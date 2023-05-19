"""
  Calculate emissions from Fire Radiative Flux (FRP/Area).
"""

__VERSION__ = 2.5
__CVSTAG__  = '@CVSTAG'

import warnings
warnings.simplefilter('ignore',DeprecationWarning)

import os
from datetime import date, timedelta

from numpy    import array, zeros, zeros_like, linspace, sum, exp,    any, all, logical_or

from gfio     import GFIO

from qfed.mxd14_l3 import __VERSION__ as L3A_VERSION
from qfed.mxd14_l3 import __CVSTAG__  as L3A_CVSTAG


#                  -------------------
#                  Internal Parameters
#                  -------------------

# Fire emission flux density [kg/s/m^2]
#
#     f_s = S_f(sat) * (FRP/A) * Alpha * B_f(s)
#

# Biome-dependent Emission Factors
# (From Andreae & Merlet 2001)
# Units: g(species) / kg(dry mater)
# --------------------------------
B_f  = {}   
eB_f = {} # errors

#                Tropical  Extratrop.    
#                 Forests     Forests   Savanna  Grasslands
#                --------  ----------   -------  ----------
B_f['CO2']  = (  1580.00,    1569.00,  1613.00,     1613.00  )
B_f['CO']   = (   104.00,     107.00,    65.00,       65.00  )
B_f['SO2']  = (     0.57,       1.00,     0.35,        0.35  )
B_f['OC']   = (     5.20,       9.14,     3.40,        3.40  )
B_f['BC']   = (     0.66,       0.56,     0.48,        0.48  )
B_f['NH3']  = (     1.30,       1.40,     1.05,        1.05  )
B_f['PM25'] = (     9.10,      13.00,     5.40,        5.40  )
B_f['TPM']  = (     8.50,      17.60,     8.30,        8.30  ) # note that TPM < PM2.5 for Tropical Forests
B_f['NO']   = (     1.60,       3.00,     3.90,        3.90  ) # NOx as NO
B_f['MEK']  = (     0.43,       0.45,     0.26,        0.26  ) # Methyl Ethyl Ketone
B_f['C3H6'] = (     0.55,       0.59,     0.26,        0.26  ) # Propene/Propylene
B_f['C2H6'] = (     1.20,       0.60,     0.32,        0.32  ) # Ethane
B_f['C3H8'] = (     0.15,       0.25,     0.09,        0.09  ) # Propane
B_f['ALK4'] = (     0.056,      0.091,    0.025,       0.025 ) # C4,5 alkanes (C4H10): n-butane + i-butane
B_f['ALD2'] = (     0.65,       0.50,     0.50,        0.50  ) # Acetaldehyde (C2H4O)
B_f['CH2O'] = (     1.40,       2.20,     0.26,        0.26  ) # Formaldehyde (HCHO)
B_f['ACET'] = (     0.62,       0.56,     0.43,        0.43  ) # Acetone (C3H6O)
B_f['CH4']  = (     6.80,       4.70,     2.30,        2.30  ) # Methene (CH4)


#                Tropical  Extratrop.    
#                 Forests     Forests   Savanna  Grasslands
#                --------  ----------   -------  ----------
eB_f['CO2']  = (   90.00,     131.00,    95.00,       95.00  )
eB_f['CO']   = (   20.00,      37.00,    20.00,       20.00  )
eB_f['SO2']  = (    0.23,       0.23,     0.16,        0.16  )
eB_f['OC']   = (    1.50,       0.55,     1.40,        1.40  )
eB_f['BC']   = (    0.31,       0.19,     0.18,        0.18  )
eB_f['NH3']  = (    0.80,       0.80,     0.45,        0.45  )
eB_f['PM25'] = (    1.50,       7.00,     1.50,        1.50  )
eB_f['TPM']  = (    2.00,       6.40,     3.20,        3.20  )
eB_f['NO']   = (    0.70,       1.40,     2.40,        2.40  )
eB_f['MEK']  = (    0.22,       0.28,     0.13,        0.13  )
eB_f['C3H6'] = (    0.25,       0.16,     0.14,        0.14  )
eB_f['C2H6'] = (    0.70,       0.15,     0.16,        0.16  )
eB_f['C3H8'] = (    0.10,       0.11,     0.03,        0.03  )
eB_f['ALK4'] = (    0.03,       0.05,     0.09,        0.09  )
eB_f['ALD2'] = (    0.32,       0.02,     0.39,        0.39  )
eB_f['CH2O'] = (    0.70,       0.50,     0.44,        0.44  )
eB_f['ACET'] = (    0.31,       0.04,     0.18,        0.18  )
eB_f['CH4']  = (    2.00,       1.90,     0.90,        0.90  )

# Scaling of C6 based on C5 (based on OC tuning)
# ----------------------------------------------
alpha = array([0.96450253,1.09728882,1.12014982,1.22951496,1.21702972])
for s in list(B_f.keys()):
    B_f[s] = list(array(B_f[s]) * alpha[1:])
    
# Combustion rate constant
# (ECMWF Tech Memo 596)
# It could be biome-dependent in case we want to tinker
# with the A-M emission factors
# -----------------------------------------------------
Alpha = 1.37e-6 # kg(dry mater)/J
A_f = {}
#                           Tropical  Extratrop.    
#                           Forests     Forests    Savanna  Grasslands
#                           --------  ----------   -------  ----------
A_f['CO2']  = Alpha * array(( 1.000,     1.000,     1.000,       1.000 ))
A_f['CO']   = Alpha * array(( 1.000,     1.000,     1.000,       1.000 ))
A_f['SO2']  = Alpha * array(( 2.500,     4.500,     1.800,       1.800 ))
A_f['OC']   = Alpha * array(( 2.500,     4.500,     1.800,       1.800 ))
A_f['BC']   = Alpha * array(( 2.500,     4.500,     1.800,       1.800 ))
A_f['NH3']  = Alpha * array(( 2.500,     4.500,     1.800,       1.800 ))
A_f['PM25'] = Alpha * array(( 2.500,     4.500,     1.800,       1.800 ))
A_f['TPM']  = Alpha * array(( 2.500,     4.500,     1.800,       1.800 ))
A_f['NO']   = Alpha * array(( 1.000,     1.000,     1.000,       1.000 ))
A_f['MEK']  = Alpha * array(( 1.000,     1.000,     1.000,       1.000 ))
A_f['C3H6'] = Alpha * array(( 1.000,     1.000,     1.000,       1.000 ))
A_f['C2H6'] = Alpha * array(( 1.000,     1.000,     1.000,       1.000 ))
A_f['C3H8'] = Alpha * array(( 1.000,     1.000,     1.000,       1.000 ))
A_f['ALK4'] = Alpha * array(( 1.000,     1.000,     1.000,       1.000 ))
A_f['ALD2'] = Alpha * array(( 1.000,     1.000,     1.000,       1.000 ))
A_f['CH2O'] = Alpha * array(( 1.000,     1.000,     1.000,       1.000 ))
A_f['ACET'] = Alpha * array(( 1.000,     1.000,     1.000,       1.000 ))
A_f['CH4']  = Alpha * array(( 1.000,     1.000,     1.000,       1.000 ))


# Satellite Fudge Factor
# ----------------------
S_f = {}
S_f['MODIS_TERRA'] = 1.385 * alpha[0] # C6 scaling based on C5 above
S_f['MODIS_AQUA' ] = 0.473

#                     Tropical  Extratrop.    
#                      Forests     Forests   Savanna  Grasslands
#                     --------  ----------   -------  ----------
#
#S_f['MODIS_TERRA'] = ( 1.000,      1.000,    1.000,       1.000 )
#S_f['MODIS_AQUA' ] = ( 1.000,      1.000,    1.000,       1.000 )



class Emissions(object):
    """
    Class for computing emissions from pre-gridded FRP
    estimates.
    """

#---
    def __init__(self, date, FRP, F, Land, Water, Cloud, biomes='default', Verb=0):
        """
        Initializes an Emission object. On input,

          date  ---   date object

          FRP   ---   Dictionary keyed by satellite name
                      with each element containing a
                      tuple with gridded Fire Radiative Power
                      (MW); each tuple element corresponds to a
                      different biome type, e.g.,
                      
                      FRP['MODIS_TERRA'] = (frp_tf,frp_xf,frp_sv,frp_gl)

          F     ---   Dictionary keyed by satellite name
                      with each element containing a
                      tuple with gridded forecast of FRP density 
                      (MW km-2); each tuple element corresponds to a
                      different biome type, e.g.,
                      
                      F['MODIS_TERRA'] = (f_tf,f_xf,f_sv,f_gl)
        
          Land  ---   Dictionary keyed by satellite name
                      with each element containing a
                      observed clear-land area [km2] for each gridbox
          
          Water ---   Dictionary keyed by satellite name
                      with each element containing a
                      water area [km2] for each gridbox

          Cloud ---   Dictionary keyed by satellite name
                      with each element containing a
                      cloud area [km2] for each gridbox            
        """

#       Save relevant information
#       -------------------------
        self.Land  = Land
        self.Water = Water
        self.Cloud = Cloud
        self.FRP   = FRP
        self.F     = F
        self.Sat   = list(Land.keys())
        self.date  = date
        self.verb  = Verb


#       Filter missing data
#       -------------------
        eps = 1.0e-2
        FillValue=1.0e20
        
        missing = []
        for sat in self.Sat:
            m = Land[sat][:,:]  > (1 - eps)*FillValue
            m = logical_or(m, Water[sat][:,:] > (1 - eps)*FillValue)
            m = logical_or(m, Cloud[sat][:,:] > (1 - eps)*FillValue)

            n_biomes = len(FRP[sat])
            for b in range(n_biomes):
                m = logical_or(m, FRP[sat][b][:,:] > (1 - eps)*FillValue)

            if any(m): 
                print('[w] Detected missing area or FRP values in %s QFED/L3A file on %s' % (sat, self.date))
           
            Land[sat][m]  = 0.0
            Water[sat][m] = 0.0 
            Cloud[sat][m] = 0.0
            for b in range(n_biomes):
                FRP[sat][b][m] = 0.0

            if all(m):
                missing.append(True)
            else:
                missing.append(False)

        assert not all(missing), '[x] No valid L3A input data. Please persist emissions from the last known good date.'
                

#       Biomes
#       -----------
        if biomes == 'default':
            self.biomes = ('Tropical Forest', 'Extratropical Forests', 'Savanna', 'Grassland')
        else:
            self.biomes = biomes[:]


#       Record grid
#       -----------
        self.im, self.jm = Land[self.Sat[0]].shape
        if (5*self.im - 8*(self.jm - 1)) == 0:
            self.lon  = linspace(-180.,180.,self.im,endpoint=False)
            self.lat  = linspace(-90.,90.,self.jm)
        elif (self.im*6 == self.jm):
            self.lon = linspace(1,self.im,self.im)
            self.lat = linspace(1,self.jm,self.jm)
        else:
            d_lon = 360.0 / self.im
            d_lat = 180.0 / self.jm
            self.lon = linspace(-180+d_lon/2, 180-d_lon/2, self.im)
            self.lat = linspace( -90+d_lat/2,  90-d_lat/2, self.jm)

#---
    def calculate(self, Species='all', method='default'):
    
        """
        Calculate emissions for each species using built-in
        emission coefficients and fudge factors.

        The default list of species is: 
        Species = ('CO', 'CO2','SO2','OC','BC','NH3','PM25','NO','MEK', 
                   'C3H6','C2H6','C3H8','ALK4','ALD2','CH2O','ACET','CH4')

        The default method for computing the emissions is:
            method = 'sequential-zero' 
        """


        if (Species == 'all') or (Species == 'default'):
            species = ('CO'  , 'CO2' , 'SO2' , 'OC'  , 'BC'  , 'NH3' , 
                       'PM25', 'NO'  , 'MEK' , 'C3H6', 'C2H6', 'C3H8', 
                       'ALK4', 'ALD2', 'CH2O', 'ACET', 'CH4')
        else:
            species = Species[:]


        # factor needed to convert B_f from [g/kg] to [kg/kg]
        units_factor = 1e-3

        n_biomes = len(self.biomes)

        A_l = zeros((self.im, self.jm))
        A_w = zeros((self.im, self.jm))
        A_c = zeros((self.im, self.jm))

        for sat in self.Sat:
            A_l += self.Land[sat]
            A_w += self.Water[sat]
            A_c += self.Cloud[sat]

        A_o = A_l + A_w

        i = (A_l > 0)
        j = ((A_l + A_c) > 0)

        E = {}
        E_= {}
        for s in species:
            E[s]  = zeros((n_biomes, self.im, self.jm))
            E_[s] = zeros((n_biomes, self.im, self.jm))

            for sat in self.Sat:
                FRP = self.FRP[sat]
                F   = self.F[sat]
                A_  = self.Cloud[sat]
                
                for b in range(n_biomes):
                    E[s][b,:,:]  += units_factor * A_f[s][b] * S_f[sat] * B_f[s][b] * FRP[b]
                    E_[s][b,:,:] += units_factor * A_f[s][b] * S_f[sat] * B_f[s][b] * F[b] * A_
           
            for b in range(n_biomes):
                E_b  = E[s][b,:,:]
                E_b_ = E_[s][b,:,:]

                if (method == 'default') or (method == 'sequential'):
                    E_b[j] = ( (E_b[j]  / (A_o[j] + A_c[j])) * (1 + A_c[j] / (A_l[j] + A_c[j])) +
                               (E_b_[j] / (A_o[j] + A_c[j])) * (    A_c[j] / (A_l[j] + A_c[j])) )

                if method == 'sequential-zero':
                    E_b[i] = E_b[i] / (A_o[i] + A_c[i]) * (1.0 + A_c[i] / (A_l[i] + A_c[i]))

                if method == 'nofires':
                    E_b[i] = E_b[i] / (A_o[i] + A_c[i])

                if method == 'similarity':
                    E_b[i] = E_b[i] / A_l[i] * ((A_l[i] + A_c[i]) / (A_o[i] + A_c[i]))

                if method == 'similarity_qfed-2.2':
                    E_b[i] = E_b[i] / A_o[i]
                   

#       Save Forecast dictionary
#       ------------------------
        dt  = 1.0    # days
        tau = 3.0    # days

        for sat in self.Sat:
            for b in range(n_biomes):
                s = species[0]

                self.F[sat][b][:,:] = (E[s][b,:,:] / (units_factor * A_f[s][b] * S_f[sat] * B_f[s][b])) * exp(-dt/tau)
                self.F[sat][b][j] = self.F[sat][b][j] * ((A_o[j] + A_c[j]) / (A_l[j] + A_c[j]))

#       Save Emission dictionary
#       ------------------------
        self.Species = species
        self.Emissions = E


#---
    def total(self, specie):
        """
        Calculates the emissions from all biomes.
        """

        return sum(self.Emissions[specie][:,:,:], axis=0)


#---
    def _write_ana(self, filename=None, dir='.', expid='qfed2', col='sfc', tag=None):
       """
       Writes gridded emissions. You must call method
       calculate() first. Optional input parameters:

       filename  ---  file name; if not specified each
                      species will be written to a separate
                      file, e.g.,
                         qfed2.emis_co.sfc.20030205.nc4
       dir       ---  optional directory name, only used
                      when *filename* is omitted
       expid     ---  optional experiment id, only used
                      when *filename* is omitted
       col       ---  collection
       tag       ---  tag name, by default it will be set to
                      the QFED CVS tag name as part of the 
                      installation procedure
       
       """
       
       title = 'QFED Level3b v%3.1f (%s) Gridded Emission Estimates' % (__VERSION__, _getTagName(tag))
       source = 'NASA/GSFC/GMAO GEOS Aerosol Group'
       contact = ('%s; %s') % ('arlindo.dasilva@nasa.gov', 'anton.darmenov@nasa.gov')


#      Create directory for output file
#      --------------------------------
       dir = os.path.join(dir, 'Y%04d'%self.date.year, 'M%02d'%self.date.month)
       rc = os.system("/bin/mkdir -p %s"%dir)
       if rc:
           raise IOError('cannot create output directory')


       nymd = 10000*self.date.year + 100*self.date.month + self.date.day
       nhms = 120000

#      Loop over species
#      -----------------
       vname_ = ( 'biomass', 'biomass_tf', 'biomass_xf', 'biomass_sv', 'biomass_gl' )

       filename_ = filename
       self.filename = {}
       for s in self.Species:

           if filename_ is None:
               filename = dir+'/%s.emis_%s.%s.%d.nc4'\
                          %(expid,s.lower(),col,nymd)
               vname = vname_ 
           else:
               vname = [ '%s_%s' % (s.lower(), name) for name in vname_]

           vtitle = [ '%s Biomass Emissions' % s, 
                      '%s Biomass Emissions from Tropical Forests' % s, 
                      '%s Biomass Emissions from Extratropical Forests' % s,
                      '%s Biomass Emissions from Savanna' % s,
                      '%s Biomass Emissions from Grasslands' % s ]
                      
           vunits = [ 'kg s-1 m-2', 'kg s-1 m-2', 'kg s-1 m-2', 'kg s-1 m-2', 'kg s-1 m-2' ]

           self.filename[s] = filename

#          Open file if it exists, otherwise create it
#          -------------------------------------------
           if os.path.exists(filename):
               f = GFIO(filename,'w')
           else:
               f = GFIO()
               f.create(filename, vname, nymd, nhms,
                        lon=self.lon, lat=self.lat,
                        vtitle=vtitle, vunits=vunits, timinc=240000,
                        title=title, source=source, contact=contact)

           f.write(vname[0], nymd, nhms, self.total(s))

           for b in range(len(self.biomes)):
               name = vname[b + 1]
               f.write(name, nymd, nhms, self.Emissions[s][b,:,:])

           try:
               f.close()
           except:
               pass

           if self.verb >=1:
               print("[w] Wrote file "+filename)


#---
    def _write_fcs(self, forecast, FillValue=1.0e20):
       """
       Writes gridded emissions. You must call method
       calculate() first. Input parameter(s):

       forecast        ---  L3a file names
       forecast_fields ---  Variable names of FRP density forecast
       """
      
       vname  = ('land', 'water', 'cloud', 
                 'frp_tf', 'frp_xf', 'frp_sv', 'frp_gl', 
                 'fb_tf', 'fb_xf', 'fb_sv', 'fb_gl')
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
                 'Background FRP Density (Grasslands)')
       vunits  = ('km2', 'km2', 'km2', 'MW', 'MW', 'MW', 'MW', 
                  'MW km-2', 'MW km-2', 'MW km-2', 'MW km-2')
       title   = 'QFED Level3a v%3.1f (%s) Gridded FRP Estimates'%(L3A_VERSION, _getTagName(L3A_CVSTAG))
       source  = 'NASA/GSFC/GMAO GEOS Aerosol Group'
       contact = ('%s; %s') % ('arlindo.dasilva@nasa.gov', 'anton.darmenov@nasa.gov')

       d = self.date + timedelta(days=1)
       nymd = 10000*d.year + 100*d.month + d.day
       nhms = 120000

       for sat in self.Sat:
           f = GFIO()
           
           f.create(forecast[sat], vname, nymd, nhms, lon=self.lon, lat=self.lat,
                    vtitle=vtitle, vunits=vunits, title=title, source=source, contact=contact)
          
           missing = zeros_like(self.F[sat][0])
           missing[:,:] = FillValue

           f.write('land',   nymd, nhms, missing)
           f.write('water',  nymd, nhms, missing)
           f.write('cloud',  nymd, nhms, missing)
           f.write('frp_tf', nymd, nhms, missing)
           f.write('frp_xf', nymd, nhms, missing)
           f.write('frp_sv', nymd, nhms, missing)
           f.write('frp_gl', nymd, nhms, missing)
         
           f.write('fb_tf', nymd, nhms, self.F[sat][0])
           f.write('fb_xf', nymd, nhms, self.F[sat][1])
           f.write('fb_sv', nymd, nhms, self.F[sat][2])
           f.write('fb_gl', nymd, nhms, self.F[sat][3])

           try:
               f.close()
           except:
               pass

#---
    def write(self, filename=None, dir='.', forecast=None, expid='qfed2', col='sfc', 
                    tag=None, ndays=1, uncompressed=False):
       """
       Writes gridded emissions that can persist for a number of days. You must 
       call method calculate() first. Optional input parameters:

       filename      ---  file name; if not specified each
                          species will be written to a separate
                          file, e.g.,
                          qfed2.emis_co.sfc.20030205.nc4
       dir           ---  optional directory name, only used
                          when *filename* is omitted
       expid         ---  optional experiment id, only used
                          when *filename* is omitted
       col           ---  collection
       tag           ---  tag name, by default it will be set to
                          the QFED CVS tag name as part of the 
                          installation procedure
       ndays         ---  persist emissions for a number of days
       uncompressed  ---  use n4zip to compress gridded output file
       
       """

#      Write out the emission files
#      ----------------------------
       self._write_fcs(forecast)

       for n in range(ndays):

#          Write out the emission files
#          ----------------------------
           self._write_ana(filename=filename,dir=dir,expid=expid,col=col,tag=tag)

#          Compress the files by default
#          -----------------------------
           if not uncompressed:
               for s in list(self.filename.keys()):
                   rc = os.system("n4zip %s"%self.filename[s])
                   if rc:
                       warnings.warn('cannot compress output file <%s>'%self.filename[s])

#          Increment date by one day
#          ---------------------------------
           self.date = self.date + timedelta(days=1)

#..............................................................

def _getTagName(tag):
    if tag != None:
        tag_name = tag
    else:    
        if __CVSTAG__ not in (None, ''):
            tag_name = __CVSTAG__
        else:
            tag_name = 'CVSTAG_UNKNOWN'

    return tag_name

