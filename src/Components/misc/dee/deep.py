"""
  Implements direct emission estimates based on MODIS Deep Blue C6
  retrievals.

"""

import os

from numpy    import linspace, array, savez
from glob     import glob
from datetime import datetime, timedelta

from pyobs.mxd04 import MxD04_L2, SDS, MISSING, ALIAS, CHANNELS
from binObs_     import binobs2d, binobs3d

from grads.gacore import gat2dt, dt2gat

SDS['META'] = ('Longitude', 'Latitude', 'Scan_Start_Time',)

SDS['DEEP'] = ( 'Deep_Blue_Spectral_Single_Scattering_Albedo_Land',
                'Deep_Blue_Angstrom_Exponent_Land',
                'Deep_Blue_Aerosol_Optical_Depth_550_Land_Best_Estimate',
                'Deep_Blue_Aerosol_Optical_Depth_550_Land_QA_Flag',
                'Deep_Blue_Algorithm_Flag_Land',)

c412, c660 = ( 0, 2 )

class DEEP(MxD04_L2):
    """
    Extends the MODIS MxD04 class with 
    """
    
    def __init__(self,Path,Verb=0):
        """
        Initializes base class MODIS MxD04_L2 class for the Deep Blue
        algorithm with dust specific filtering based on the Angstrom
        exponent. 
        """

        # Initialize base class
        # ---------------------
        MxD04_L2.__init__(self,Path,'DEEP',Verb=Verb,only_good=True)

        self.nch = len(CHANNELS)

        # Screening for Dust
        #     angstrom in [-5,1]
        #     ssa(670)>ssa(412)
        # ----------------------
        r = self.ssa[:,c660]/self.ssa[:,c412]
        I = (r>1.) & (self.angstrom<1.) & (self.atype==0) & (self.aod550>-0.01)
        self.reduce(I)  # keep only dust obs

#---
    def write(self,filename=None,dirn='.',expid=None,Verb=1):
        """
        Writes the un-gridded reduced object to a numpy npz file. 
        """

        # Stop here is no good obs available
        # ----------------------------------
        if self.nobs == 0:
            return # no data to work with
        if any(self.iGood) == False:
            return # no good data to work with

        if expid == None:
            expid = self.ident

        if filename is None:
            filename = '%s/%s.dudb.%d.npz'%(dirn,expid,self.nymd/100)

        version = 1 # File format version
        meta = [self.nymd,self.nhms,self.nobs,self.nch,version]
        savez(filename,
                            meta = meta,
                             lon = self.lon,
                             lat = self.lat,
                            time = self.Time,
                        channels = self.channels,
                         qa_flag = self.qa_flag,
                          aod550 = self.aod550,
                        angstrom = self.angstrom,
                             ssa = self.ssa)

        if Verb >=1:
            print("[w] Wrote file "+filename)

#---
    def writeg(self,filename=None,dirn='.',expid=None,refine=8,res=None,
               channels=(412,470,660),Verb=1):
       """
        Writes gridded MODIS measurements to file.

         refine  -- refinement level for a base 4x5 GEOS-5 grid
                       refine=1  produces a   4  x  5    grid
                       refine=2  produces a   2  x2.50   grid
                       refine=4  produces a   1  x1,25   grid
                       refine=8  produces a  0.50x0.625  grid
                       refine=16 produces a  0.25x0.3125 grid
        Alternatively, one can specify the grid resolution with a
        single letter:

         res     -- single letter denoting GEOS-5 resolution,
                       res='a'  produces a   4  x  5    grid
                       res='b'  produces a   2  x2.50   grid
                       res='c'  produces a   1  x1,25   grid
                       res='d'  produces a  0.50x0.625  grid
                       res='e'  produces a  0.25x0.3125 grid

                   NOTE: *res*, if specified, supersedes *refine*.

         Verb -- Verbose level:
                 0 - really quiet (default)
                 1 - Warns if invalid file is found
                 2 - Prints out non-zero number of fires in each file.


       """
       from gfio import GFIO
       
       # Stop here is no good obs available
       # ----------------------------------
       if self.nobs == 0:
           return # no data to work with

       if expid == None:
           expid = self.ident

#      Output grid resolution
#      ----------------------
       if res is not None:
           if res=='a': refine = 1 
           if res=='b': refine = 2
           if res=='c': refine = 4
           if res=='d': refine = 8
           if res=='e': refine = 16

#      Lat lon grid
#      ------------
       dx = 5. / refine
       dy = 4. / refine
       im = int(360. / dx)
       jm = int(180. / dy + 1)

       glon = linspace(-180.,180.,im,endpoint=False)
       glat = linspace(-90.,90.,jm)

       if channels is None:
           channels = self.channels
 
       levs = array(channels)

       nch = len(channels)
       nymd = self.nymd
       nhms = self.nhms

       vtitle = [ 'Aerosol Optical Depth [550 nm]',\
                  'Angstrom Exponent',\
                  'Single Scattering Albedo']

       vname  = ['aod', 'ang','ssa']
       vunits = [ '1','1', '1'] 
       kmvar  = [  0, 0, nch ]

       title = 'Gridded MODIS Aerosol Retrievals'
       source = 'NASA/GSFC/GMAO GEOS-5 Aerosol Group'
       contact = 'arlindo.dasilva@nasa.gov'

       if filename is None:
           filename = '%s/%s.dudb.%d_%02dz.nc4'%(dirn,expid,self.nymd,self.nhms/10000)

       # Create the file
       # ---------------
       f = GFIO()
       f.create(filename, vname, nymd, nhms,
                lon=glon, lat=glat, levs=levs, levunits='nm',
                vtitle=vtitle, vunits=vunits,kmvar=kmvar,amiss=MISSING,
                title=title, source=source, contact=contact)

       # Grid variable and write to file
       # -------------------------------
       f.write('aod', nymd, nhms, 
               binobs2d(self.lon,self.lat,self.aod550,im,jm,MISSING) )
       f.write('ang', nymd, nhms, 
               binobs2d(self.lon,self.lat,self.angstrom,im,jm,MISSING) )
       f.write('ssa', nymd, nhms, 
               binobs3d(self.lon,self.lat,self.ssa,im,jm,MISSING) )
           
       try:
           f.close()
       except:
           pass
     
       os.system('n4zip -v '+filename)

       if Verb >=1:
           print("[w] Wrote file "+filename)

#....................................................................
def GranulesSyn ( path, prod, syn_time, coll='006', nsyn=24):
    """
    Returns a list of MxD04 granules for a given product at given synoptic time.
    On input,

    path      ---  mounting point for the MxD04 Level 2 files
    prod      ---  either MOD04 or MYD04
    syn_time  ---  synoptic time (timedate format)

    coll      ---  collection: 005, 051 (optional)
    nsyn      ---  number of synoptic times per day (optional)

    """

    # Determine synoptic time range
    # -----------------------------
    dt = timedelta(seconds = 12. * 60. * 60. / nsyn)
    t1, t2 = (syn_time-dt,syn_time+dt)

    return Granules_(path,prod,coll,t1,t2)

#....................................................................
def mrange(year,month,gat=False):
    """
    Returns first and last time for a given month.
    Assumes hourly time steps. Returns GrADS time if gat=True
    """
    day   = timedelta(seconds=24*60*60)
    hour  = timedelta(seconds=60*60)
    hhour = timedelta(seconds=30*60)

    t1 = datetime(year,month,1)  # first of the month
    if month==12:
        t2 = datetime(year+1,1,1) - day
    else:
        t2 = datetime(year,month+1,1) - day
    t1 += hhour              #  0:30Z
    t2 += 23 * hour + hhour  # 23:30Z

    if gat:
        return (dt2gat(t1), dt2gat(t2))
    else:
        return (t1, t2)

def GranulesMonth ( path, prod, year, month, coll='006'):
    """
    Returns a list of MxD04 granules for a given product at given synoptic time.
    On input,

    path      ---  mounting point for the MxD04 Level 2 files
    prod      ---  either MOD04 or MYD04
    year      ---  year, e.g., 2003
    month     ---  month, 1..12
    coll      ---  collection: 006, 005, 051 (optional)

    """

    t1 = datetime(year,month,1) # first of the month
    if month==12:
        t2 = datetime(year+1,1,1) 
    else:
        t2 = datetime(year,month+1,1) 

    return Granules_(path,prod,coll,t1,t2)

def Granules_(path,prod,coll,t1,t2):
    """
    Find Granules in [t1,t2)
    """

    # Find MODIS granules in synoptic time range
    # ------------------------------------------
    dt = timedelta(minutes=5)
    t = t1
    Granules = []
    while t < t2:
        if t >= t1:
            doy = t.timetuple()[7]
            # print t, doy
            basen = "%s/%s/%04d/%03d/%s_L2.A%04d%03d.%02d%02d.%s.*.hdf"\
                     %(path,prod,t.year,doy,prod,t.year,doy,t.hour,t.minute,coll)
            try:
                filen = glob(basen)[0]
                Granules += [filen,]
#               print " [x] Found "+filen
            except:
                pass
        t += dt

    if len(Granules) == 0:
        print("WARNING: no %s collection %s granules found for"%(prod,coll), t1, t2)

    return Granules

#....................................................................

def select_dust(path,year1,year2,dirn='./DEE'):

    for year in range(year1,year2+1):
        for month in range(1,13):
            for prod in ('MOD04', 'MYD04'):

                Files = GranulesMonth(path,prod,year,month)

                d = DEEP(Files,Verb=0)
   
                d.nymd = year*10000 + month*100 + 00
                d.nhms = -1
                
                d.write(dirn=dirn,expid='dee_%s'%d.ident)

#....................................................................

def grid_dust(path,year1,year2,dirn='./DEE'):
    """
    Write hourly gridded files.
    """

    hour = timedelta(seconds=60*60)
    hhour = timedelta(seconds=30*60)

    for year in range(year1,year2+1):
        for month in range(1,13):
            for prod in ('MOD04', 'MYD04'):

                dirm = dirn+'/Level3/%s/Y%d/M%02d'%(prod,year,month)
                os.system('/bin/mkdir -p '+dirm)

                t1, t2 = mrange(year,month) # 0:30Z to 23:30Z 

                t  = t1
                while t <= t2:

                    Files = GranulesSyn(path,prod,t,coll='006',nsyn=24)

                    if len(Files)>0:

                        d = DEEP(Files,Verb=0)
   
                        d.nymd = t.year*10000 + t.month*100 + t.day
                        d.nhms = t.hour*10000 + t.minute*100 + t.second
         
                        filename = '%s/dee_%s.dudb.%d_%02d:%02dz.nc4'%\
                            (dirm,prod,d.nymd,t.hour,t.minute)

                        d.writeg(filename=filename,Verb=False)

                    t += hour

#....................................................................

if __name__ == "__main__":

    path = '/nobackup/MODIS/006/Level2'

    #select_dust(path,year1,year2,dirn='./DEE')
    grid_dust(path,2008,2008,dirn='./DEE')
    
