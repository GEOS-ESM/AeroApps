#!/usr/bin/env python

"""
  Uses grads to compute merra-2 anomalies.

"""

import os
import sys

from numpy    import linspace, array, savez
from glob     import glob
from datetime import datetime, timedelta
from string   import Template

from grads import GrADS

class myGrADS(GrADS):

    def Cmd(self,tmpl,d):
        """
        Wrapper for cmd() 
        """
        # print Template(tmpl).substitute(d)
        self.__call__(Template(tmpl).substitute(d))

#....................................................................

def clim_anom(ga,coll,var,year1,year2,dirn='./Anom'):

    cFile = dirn+'/MERRA2_clm.%s.%d-%d.nc4'%(var,year1,year2)
    aFile = dirn+'/MERRA2_ano.%s.%d-%d.nc4'%(var,year1,year2)

    ga("set time jan%d dec%d"%(year1,year2))
    qh = ga.query("dims")

    d = dict ( var=var, y1=year1, y2=year2, t1=qh.t[0],t2=qh.t[1],
               aFile=aFile, cFile=cFile)

    ga.Cmd("""
           set time jan$y1 dec$y1
           define ${var}clm = ave($var,t+0,t=$t2,12)
           modify ${var}clm seasonal

           set time jan$y1 dec$y2
           define ${var}ano = $var - ${var}clm

           set sdfwrite -flt -nc4 -zip $cFile
           sdfwrite ${var}clm

           set sdfwrite -flt -nc4 -zip $aFile
           sdfwrite ${var}ano
 
           undefine ${var}clm
           undefine ${var}ano

           """,d)            

#....................................................................

if __name__ == "__main__":

    inPath = '/home/adasilva/opendap/m2/opendap'

    # Parse command line
    # ------------------
    if len(sys.argv) == 3:
        year1 = sys.argv[1]
        year2 = sys.argv[2]
    elif len(sys.argv) == 2 :
        year1 = sys.argv[1]
        year2 = year1
    else:
        year1 = 2000
        year2 = 2015

    year1, year2 = int(year1), int(year2)    

    # Instantiate GrADS
    # -----------------
    ga = myGrADS(Window=False,Echo=False)
       

    #for varcoll in ( 'tavgM_2d_slv_Nx:ts,slp,u10m,v10m',  \
    #                 'tavgM_2d_aer_Nx:ducmass,dufluxu,dufluxv' ):
    ###for varcoll in ( 'tavgM_2d_aer_Nx:ducmass,so4cmass,bccmass,occmass,so2cmass,dusmass,so4smass,bcsmass,ocsmass,so2smass,duexttau,suexttau,bcexttau,ocexttau', ):
    for varcoll in ( 'tavgM_2d_aer_Nx:sscmass,sssmass,ssexttau',):

        coll,Vars = varcoll.split(':')
        Vars = Vars.split(',')
    
        fh = ga.open("%s/%s"%(inPath,coll)) # open collection
        ga("set x 1 %d"%fh.nx)

        for var in Vars:

            print " [] Working on <%s>"%var

            clim_anom(ga,coll,var,year1,year2)

