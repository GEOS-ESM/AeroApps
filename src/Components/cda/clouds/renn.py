"""
    Neural Net estimattion of cloud effective radius.
"""

from datetime import datetime
from numpy    import savez, linspace, ones, where, log
from numpy    import c_ as cat

from pyobs.mxd06 import MxD06_L2, MISSING, granules
from pyobs.npz   import NPZ

from Aero_       import getaeropbl, getaerolev, revisectp, getrh

from gfio        import GFIOctl

LIQUID, ICE = 2, 3

class RENN(MxD06_L2):

    def __init__(self,Path,kbeg=37,eRE=10.,**kwopts):
        """
        Instantiate a Re (Cloud Effective Radius) Neural Net
        object. 

        On input, *kbeg* is the top GEOS-5 vertical layer (Fortran
        1-offset); kbeg=37 corresponts to ~ 79 hPa. This device is
        used to control reading of 3D GEOS-5 fields.

        Multi-layer clouds and data with effective radius grater 
        than *eRE* are filtered out.

        """

        MxD06_L2.__init__(self,Path,**kwopts)

        self.iVal = (self.MLAYER==1) & (self.RE>0) & (self.eRE<=eRE)
        self.iLiq = self.iVal & (self.PHASE==LIQUID)
        self.iIce = self.iVal & (self.PHASE==ICE)

        self.kbeg = kbeg             # top layer, fortran 1-offset
        self.kount = 72 - kbeg + 1

    def getVcoords(self,nav_Ne,flx_Nx,asm_Nx):
        """
        Get GEOS-5 midlevel height and pressure coordinates, as well as
        PBL height. On input, GEOS-5 standard collections.
        """

        lon, lat, time = (self.lon.ravel(), self.lat.ravel(), self.Time.ravel())

        # PBL Height
        # ----------
        f = GFIOctl(flx_Nx)
        self.PBLH = f.sample('PBLH',lon,lat,time,Verbose=self.verb)

        # T, QV, PS
        # ---------
        f = GFIOctl(asm_Nx)
        self.PS   = f.sample('PS',  lon,lat,time,Verbose=self.verb)
        self.T2M  = f.sample('T2M', lon,lat,time,Verbose=self.verb)
        self.QV2M = f.sample('QV2M',lon,lat,time,Verbose=self.verb)

        # Reshape - 2D
        # ------------
        shp2d = self.lon.shape
        self.PS   = self.PS.reshape(shp2d)
        self.T2M  = self.T2M.reshape(shp2d)
        self.QV2M = self.QV2M.reshape(shp2d)
        self.PBLH = self.PBLH.reshape(shp2d)

        # T at LCL
        # --------
        self.TLCL = _getTLCL(self.T2M,self.QV2M,self.PS)

        # Height & Pressure
        # -----------------
        f = GFIOctl(nav_Ne)
        self.HE = f.sample('ZLE',lon,lat,time,kbeg=self.kbeg,kount=self.kount+1,
                           Verbose=self.verb)
        self.PE = f.sample('PLE',lon,lat,time,kbeg=self.kbeg,kount=self.kount+1,
                           Verbose=self.verb)

        # Reshape - 3D
        # ------------
        shp3d = (self.kount+1,) + self.lon.shape
        self.HE   = self.HE.reshape(shp3d)
        self.PE   = self.PE.reshape(shp3d)

        # Revised CTP
        # -----------
        rCTP = revisectp(self.CTP*100,self.CTT,self.PE,self.HE,self.PBLH)
        I = (self.CTP<0) | (self.CTT<0)
        self.rCTP = where(I,MISSING,rCTP/100.) # units: hPa

    def getAero(self, aer_Nv, npzFile=None,
                Vars=('DU001','SS001','BCPHILIC','OCPHILIC','SO4')):
        """
        Get aerosol predictors.
        """
        
        shp3d = (self.kount,) + self.lon.shape
        lon, lat, time = self.lon.ravel(), self.lat.ravel(), self.Time.ravel()
        f = GFIOctl(aer_Nv)
        
        Features = dict(rCTP = self.rCTP, TLCL=self.TLCL)

        for var in Vars:
           if self.verb:
               print(" <> Getting aerosol features for <%s> "%var)

           Q = f.sample(var,lon,lat,time,kbeg=self.kbeg,kount=self.kount,
                        Verbose=self.verb)
           Q = Q.reshape(shp3d)

           QABL, QPBL = getaeropbl(Q,self.PE,self.HE,self.PBLH)
           QLEV = getaerolev(Q,self.PE,self.HE,self.rCTP*100)

           if var[0:2] == 'SO':
              prefix = 'SU'
           else:
              prefix = var[0:2]

           self.__dict__[prefix+'LEV'] = QLEV
           self.__dict__[prefix+'ABL'] = QABL
           self.__dict__[prefix+'PBL'] = QPBL

           if npzFile is not None:
               Features[prefix+'LEV'] = QLEV
               Features[prefix+'ABL'] = QABL
               Features[prefix+'PBL'] = QPBL

        # Write out npzFile if desidered
        # ------------------------------
        if npzFile is not None:
            savez(npzFile,**Features)

    def loadAero(self, npzFile=None):
        """
        Load aerosol features from an NPZ file.
        """
        f = NPZ(npzFile)
        for var in f.__dict__:
            self.__dict__[var] = f.__dict__[var]

    def getInputs(self,I,Input):
        """
        Given a set of indices *I*, returns the corresponding
        inputs for a neural net evaluation.
        Returns: inputs
        """
        if self.verb:
            print(" ")
            print("       Feature          Min      Max")
            print("  ------------------  -------  -------")
        inputs = self.__dict__[Input[0]][I]
        if self.verb:
            print("%20s %8.4f %8.4f"%(Input[0],inputs.min(),inputs.max()))
        for var in Input[1:]:
            q = self.__dict__[var][I]
            inputs = cat[inputs,q]
            if self.verb:
                print("%20s %8.4f %8.4f"%(var,q.min(),q.max()))
        if self.verb:
            print("  ------------------  -------  -------")
            print("")
        return inputs
    
    def getTargets(self,I,Target=('RE',)):
        """
        Given a set of indices *I*, return the corresponding
        targets for a neural net evaluation:
        Returns: tagets
        """
        targets = self.__dict__[Target[0]][I]
        if self.verb:
            print(" ")
            print("       Target           Min      Max")
            print("  ------------------  -------  -------")
            print("%20s %8.4f %8.4f"%(Target[0],targets.min(),targets.max()))
        for var in Target[1:]:
            targets = cat[targets,self.__dict__[var][I]]
            if self.verb:
                print("%20s %8.4f %8.4f"%(var,targets.min(),targets.max()))
        return targets

#---
    def writeAero(self,filename, syn_time, iFilter=None,
                  refine=8,res=None, Verb=1):
       """
        Writes gridded aerosol parameters that have been co-located
        with MODIS coud retrievals. On input,

         phase  --  cloud phase; one of (None,'Liquid','Ice')
                    None means both phases.

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
       
       # Optional filter
       # ---------------
       self.iFilter = iFilter

       # Output grid resolution
       # ----------------------
       if res is not None:
           if res=='a': refine = 1 
           if res=='b': refine = 2
           if res=='c': refine = 4
           if res=='d': refine = 8
           if res=='e': refine = 16

       # Lat lon grid
       # ------------
       dx = 5. / refine
       dy = 4. / refine
       im = int(360. / dx)
       jm = int(180. / dy + 1)

       glon = linspace(-180.,180.,im,endpoint=False)
       glat = linspace(-90.,90.,jm)

       nch = 1
       levs = ones(1) # for now

       t = syn_time
       Y, M, D, h, m, s = (t.year,t.month,t.day,t.hour,t.minute,t.second)

       nymd = Y*10000+M*100+D
       nhms = h*10000+m*100+s

       vtitle = [ 'Dust Concentration at Cloud Top',
                  'Sea Salt Concentration at Cloud Top',
                  'Black Carbon Concentration at Cloud Top',
                  'Organic Carbon Concentration at Cloud Top',
                  'Sulfate Concentration at Cloud Top',
                  'Dust Average Concentration Above Boundary Layer',
                  'Sea Salt Average Concentration Above Boundary Layer',
                  'Black Carbon Average Concentration Above Boundary Layer',
                  'Organic Carbon Average Concentration Above Boundary Layer',
                  'Sulfate Average Concentration Above Boundary Layer',
                  'Dust Average Concentration Within Boundary Layer',
                  'Sea Salt Average Concentration Within Boundary Layer',
                  'Black Carbon Average Concentration Within Boundary Layer',
                  'Organic Carbon Average Concentration Within Boundary Layer',
                  'Sulfate Average Concentration Within Boundary Layer',
                 ]

       vname  = [ 'DULEV', 'SSLEV', 'BCLEV', 'OCLEV', 'SULEV',
                  'DUABL', 'SSABL', 'BCABL', 'OCABL', 'SUABL',
                  'DUPBL', 'SSPBL', 'BCPBL', 'OCPBL', 'SUPBL',
                ]
       vunits = 15 * [ 'g cm-3', ]
       kmvar  = 15 * [  0, ]

       title = 'Gridded GEOS-5 Aerosol Co-located with MODIS Cloud Retrievals'
       source = 'NASA/GSFC/GMAO'
       contact = 'arlindo.dasilva@nasa.gov'

       if filename is None:
           filename = '%s/%s.aer.%d_%02dz.nc4'%(dir,expid,self.nymd,self.nhms/10000)

       # Create the file
       # ---------------
       f = GFIO()
       f.create(filename, vname, nymd, nhms,
                lon=glon, lat=glat, levs=levs, levunits='nm',
                vtitle=vtitle, vunits=vunits,kmvar=kmvar,amiss=MISSING,
                title=title, source=source, contact=contact)

       # Grid variables and write to file
       # -------------------------------
       for name in vname:
           f.write (name, nymd, nhms, self._binobs2d(self.__dict__[name],im,jm))

       f.close()

       # Reset filter
       # ------------
       self.iFilter = None

       if Verb >=1:
           print("[w] Wrote file "+filename)

#...............................................................................

def _doNN(m):

    #from pyobs import sknet as nn
    import ffnet as nn

    Input_pbl  = ( 'CTT', 
                   'DUPBL', 'SSPBL', 'BCPBL', 'OCPBL', 'SUPBL', 
                 )
    Input_abl  = ( 'CTT', 
                   'DUPBL', 'SSPBL', 'BCPBL', 'OCPBL', 'SUPBL', 
                   'DUABL', 'SSABL', 'BCABL', 'OCABL', 'SUABL', 
                 )
    Input_lev  = ( 'CTT', 
                   'DULEV', 'SSLEV', 'BCLEV', 'OCLEV', 'SULEV', 
                 )
    Target = ('RE',)

    XXLEV = m.DULEV+m.SSLEV+m.BCLEV+m.OCLEV+m.SULEV

    # Water clouds
    # ------------
    for Input in ( Input_lev, ):

        I = m.iIce & ( m.eRE<10 ) 

        X = m.getInputs(I,Input)
        y = m.getTargets(I,Target)

        y = 1. / y

        topology = nn.mlgraph( (len(Input),len(Input),1) )

        # Instantiate Neural Net
        # ----------------------
        #net = nn.SKNET(topology)
        net = nn.ffnet(topology)

        print(" <> Starting training with %s inputs and %s targets"\
                  %(str(X.shape),str(y.shape)))

        net.train_tnc(X,y,maxfun=2500)

        y_, reg = net.test(X,y)

    return (net,X,y)

def _getTLCL(T, QV, PS):
    """
    Returns temperature at LCL.
    T --- T at 2m in Kelvin
    QV -- QV at 2m [0-1]
    PS -- surface pressure
    """

    RH = getrh (QV,T,PS)
    Tc = T - 273.16 # to Celsius
    Td = T - (14.55+0.114*Tc)*(1-RH) + \
             ((2.5+0.007*Tc)*(1-RH))**3 + \
             (15.9+0.117*Tc)*(1-RH)**14

    Tlcl = 1 / ( 1/(Td-56) + log(T/Td)/800 ) + 56 

    return Tlcl
    

#...................................................................................

if __name__ == "__main__":

      from glob import glob

      nav_Ne = '/nobackup/fp/opendap/nav_Ne.ctl'
      flx_Nx = '/nobackup/fp/opendap/flx_Nx.ctl'
      asm_Nx = '/nobackup/fp/opendap/asm_Nx.ctl'
      aer_Nv = '/nobackup/fp/opendap/aer_Nv.ctl'

      syn_time = datetime(2012,8,15,12,0,0)
      Files = granules('/nobackup/MODIS/051/Level2','MOD06',syn_time,coll='051')
      Files = sorted(glob('/nobackup/MODIS/051/Level2/MOD06/2012/228/MOD06_L2.A2012228.*.hdf'))
      
      m = RENN(Files,Verb=1)

      m.getVcoords(nav_Ne,flx_Nx,asm_Nx)
      m.getAero(aer_Nv, npzFile='aeroFeatures.npz')

      #m.loadAero('aeroFeatures.npz')

def hold():

      m.write('liq-clouds.nc4',syn_time,iFilter=m.iLiq)
      m.write('ice-clouds.nc4',syn_time,iFilter=m.iIce)
      m.writeAero('liq-aerosols.nc4',syn_time,iFilter=m.iLiq)
      m.writeAero('ice-aerosols.nc4',syn_time,iFilter=m.iIce)

      #net, X, y = _doNN(m)



      
