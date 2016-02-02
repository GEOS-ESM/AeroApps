"""
  Implements direct emission estimates based on MODIS Deep Blue C6
  retrievals.

"""

import os
import sys

from numpy    import linspace, array, savez
from datetime import datetime, timedelta
from string   import Template

from grads        import GrADS
from grads.gacore import gat2dt, dt2gat

Force = False

class myGrADS(GrADS):

    def Cmd(self,tmpl,d):
        """
        Wrapper for cmd() 
        """
        self.__call__(Template(tmpl).substitute(d))

Sat = dict ( MYD04='Aqua', MOD04='Terra')

#....................................................................
def plot_xy(path,expr,var,Var,clevs,cmin,scale,year1=2003,year2=2014):
    """
    Plot 12 month climatology.
    """
    dirn = path+'/Images/'+var
    d = dict(path=path, var=var, Var=Var, clevs=clevs,
             dirn = dirn, cmin=cmin, scale=scale, expr=expr,
             year1=str(year1), year2=str(year2))

    os.system('/bin/mkdir -p '+dirn)

    ga = myGrADS(Window=False)

    #for prod in ('MOD04', 'MYD04', 'MxD04'):
    for prod in ('MxD04', ):

        d['prod'] = prod
        d['sat'] = Sat[prod]

        ga.Cmd("""
               reinit
               open $path/opendap/$prod.$var
               set mpdset mres
               set gxout grfill
               set lon -180 180
               """,d)

        for month in range(1,13):
            M3 = datetime(year1,month,1).ctime().split()[1]
            d['m2'] = '%02d'%month
            d['M3'] = M3
            d['t1'] = M3+str(year1)
            d['t2'] = M3+str(year2)
            ga.Cmd("""
                   c
                   set grads off
                   set clevs $clevs
                   xxx = $scale*ave($expr,time=$t1,time=$t2,12)
                   d if(xxx,<,$cmin,-u,xxx)
                   cbarn
                   draw title $Var - $sat $M3 Clim. ($year1-$year2)
                   gxyat -x 1600 -y 1200 $dirn/dee_$prod.$var.$m2.png
                   """,d)


#....................................................................
def plot_DEE_xy(path,var,Var,clevs,cmin,scale,
                fmin=0.5,year1=2003,year2=2014,cdiv=1,crm=1):
    """
    Plot 12 month climatology.
    """
    dirn = path+'/Images/'+var
    d = dict(path=path, var=var, Var=Var, clevs=clevs,
             dirn = dirn, cmin=cmin, scale=scale, fmin=fmin,
             cdiv=cdiv, crm=crm,
             year1=str(year1), year2=str(year2))

    os.system('/bin/mkdir -p '+dirn)

    ga = myGrADS(Window=False)

    for prod in ('MOD04', 'MYD04'):

        d['prod'] = prod
        d['ndvi'] = prod[0:3]+'43.ndvi'
        d['sat'] = Sat[prod]

        ga.Cmd("""
               reinit
               open $path/opendap/$prod.uqvq
               open $path/opendap/$prod.durm
               open $path/opendap/$prod.foo
               open $path/opendap/$ndvi
               nbare = 0.4 * 10000 
               set mpdset mres
               set gxout grfill
               set lon -180 180
               """,d)

        for month in range(1,13):
            M3 = datetime(year1,month,1).ctime().split()[1]
            d['m2'] = '%02d'%month
            d['M3'] = M3
            d['t1'] = M3+str(year1)
            d['t2'] = M3+str(year2)
            ga.Cmd("""
                   c
                   set grads off
                   set clevs $clevs

                   dee = $scale*ave(if(ndvi.4,<,0,-u,if(ndvi.4,>,nbare,-u,1))*if(foo.3,>,$fmin/100,1,-u)*($cdiv*hdivg(dufluxu,dufluxv)+$crm*durm.2),time=$t1,time=$t2,12)
 
                   d if(dee,<,$cmin,-u,dee)
                   cbarn
                   draw title $Var - $sat $M3 Clim. ($year1-$year2)
                   gxyat -x 1600 -y 1200 $dirn/dee_$prod.$var.$m2.png
                   """,d)


#....................................................................

if __name__ == "__main__":

    path = '/nobackup/5/GAAS/DEE/'

    """
    # FoO
    # ---
    expr, var, Var = 'foo', 'foo', 'FoO [%]'
    clevs = '0.5 0.75 1.0 1.25 1.5 1.75 2 2.25 2.5 2.75 3 3.5 4 4.5 5'
    cmin = 0.5
    scale = 100.
    plot_xy(path,expr,var,Var,clevs,cmin,scale) 

    # MERRA-2 Dust Emissions
    # ----------------------
    expr, var, Var = 'duem', 'duem', 'MERRA-2 Emissions [ug/m2/s]'
    clevs = '0.1 0.14 0.2 0.28 0.4 0.57 0.8 1.13 1.6 2.3 3.2 4.5 6.4 9 12.8'
    cmin = 0.001
    scale = 1e9
    plot_xy(path,expr,var,Var,clevs,cmin,scale) 

    # Wx. Modulation Factor
    # ---------------------
    expr, var, Var = 'duwx', 'duwx', 'Wx. Modulation'
    clevs = ' 1 2 4 6 8 32 64 128 512 1024 2048 4096 8192'
    cmin = 0.00
    scale = 1
    plot_xy(path,expr,var,Var,clevs,cmin,scale) 

    # Deep Blue DOD
    # -------------
    expr, var, Var = 'aod', 'aod_o', 'Deep Blue DOD [550 nm]'
    clevs = '0.01 0.02 0.04 0.08 0.16 0.32 0.64 1.28 2.56'
    cmin = 0.00
    scale = 1
    plot_xy(path,expr,var,Var,clevs,cmin,scale) 

    # MERRA-2 DOD
    # -------------
    expr, var, Var = 'duexttau','aod_m', 'MERRA-2 DOD [550 nm]'
    clevs = '0.01 0.02 0.04 0.08 0.16 0.32 0.64 1.28 2.56'
    cmin = 0.00
    scale = 1
    plot_xy(path,expr,var,Var,clevs,cmin,scale) 

    # Direct Emission Estimate
    # ------------------------
    var, Var = 'dee', 'D.E.E. [ug/m2/s]'
    clevs = '0.1 0.14 0.2 0.28 0.4 0.57 0.8 1.13 1.6 2.3 3.2 4.5 6.4 9 12.8'
    cmin = 0.
    scale = 1e9
    fmin = 0.1 # foo threshold to set to undef
    plot_DEE_xy(path,var,Var,clevs,cmin,scale,fmin=fmin) 

    # Direct Emission Estimate: 
    # Divergence Contribution
    # ------------------------
    var, Var = 'dee_d', 'D.E.E. Div [ug/m2/s]'
    clevs = '0.1 0.14 0.2 0.28 0.4 0.57 0.8 1.13 1.6 2.3 3.2 4.5 6.4 9 12.8'
    cmin = 0.
    scale = 1e9
    fmin = 0.1 # foo threshold to set to undef
    plot_DEE_xy(path,var,Var,clevs,cmin,scale,fmin=fmin,crm=0) 

    # Direct Emission Estimate:
    # Removal contribution
    # -------------------------
    var, Var = 'dee_r', 'D.E.E. Rm [ug/m2/s]'
    clevs = '0.1 0.14 0.2 0.28 0.4 0.57 0.8 1.13 1.6 2.3 3.2 4.5 6.4 9 12.8'
    cmin = 0.
    scale = 1e9
    fmin = 0.1 # foo threshold to set to undef
    plot_DEE_xy(path,var,Var,clevs,cmin,scale,fmin=fmin,cdiv=0) 

    """
