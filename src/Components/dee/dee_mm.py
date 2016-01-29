"""
  Implements direct emission estimates based on MODIS Deep Blue C6
  retrievals.

"""

import os
import sys

from numpy    import linspace, array, savez
from glob     import glob
from datetime import datetime, timedelta
from string   import Template

from grads        import GrADS
from grads.gacore import gat2dt, dt2gat

from deep import mrange

Force = True

#....................................................................
def mflux_o(path,filename,prod,gatime):
    """
    Given a time range in gatime, use lats4d to compute mass flux month
    mean.

    path       top path for input control files, see below
    filename   output file name

    """

    if not Force:
        if os.path.exists(filename):
            print '<> File exists, skipping <%s>'%filename
            return

    d = dict(path=path, prod=prod, filename=filename, 
             t1=gatime[0], t2=gatime[1])

    tmpl = """
           lats4d.sh -v -gzip 2 \
                -i    $path/opendap/du_cm \
                -j    $path/opendap/$prod \
                -o    $filename \
                -vars dufluxu dufluxv \
                -time $t1 $t1 \
                -func 'ave(@*aod.2*(ducmass.1/duexttau.1),time=$t1,time=$t2)'
          """

    cmd = Template(tmpl).substitute(d)

    rc = os.system(cmd)

    if rc:
        raise RuntimeError, 'error on return from %s'%cmd

def mflux_m(path,filename,prod,gatime):
    """
    Given a time range in gatime, use lats4d to compute mass flux month
    mean.

    path       top path for input control files, see below
    filename   output file name

    """

    if not Force:
        if os.path.exists(filename):
            print '<> File exists, skipping <%s>'%filename
            return

    d = dict(path=path, prod=prod, filename=filename, 
             t1=gatime[0], t2=gatime[1])

    tmpl = """
           lats4d.sh -v -gzip 2 \
                -i    $path/opendap/du_cm \
                -j    $path/opendap/$prod \
                -o    $filename \
                -vars dufluxu dufluxv \
                -time $t1 $t1 \
                -func 'ave(@+0*aod.2,time=$t1,time=$t2)'
 `         """

    cmd = Template(tmpl).substitute(d)

    if os.system(cmd):
        raise RuntimeError, 'error on return from %s'%cmd

#....................................................................
def removal(path,filename,prod,gatime):
    """
    Given a time range in gatime, use lats4d to compute data constrained
    removal month mean.

    path       top path for input control files, see below
    filename   output file name

    """

    if not Force:
        if os.path.exists(filename):
            print '<> File exists, skipping <%s>'%filename
            return

    d = dict(path=path, prod=prod, filename=filename, 
             t1=gatime[0], t2=gatime[1])

    tmpl = """
           lats4d.sh -v -gzip 2 \
                -i    $path/opendap/du_rm \
                -j    $path/opendap/du_cm \
                -k    $path/opendap/$prod \
                -o    $filename \
                -vars durm \
                -time $t1 $t1 \
                -func 'ave(durm.1*aod.3/duexttau.2,time=$t1,time=$t2)'
          """

    cmd = Template(tmpl).substitute(d)

    rc = os.system(cmd)

    if rc:
        raise RuntimeError, 'error on return from %s'%cmd

#....................................................................
def aod_o(path,filename,prod,gatime):
    """
    Given a time range in gatime, use lats4d to compute

    path       top path for input control files, see below
    filename   output file name

    """

    if not Force:
        if os.path.exists(filename):
            print '<> File exists, skipping <%s>'%filename
            return

    d = dict(path=path, prod=prod, filename=filename, 
             t1=gatime[0], t2=gatime[1])

    if prod=='MxD04':
        tmpl = """
               lats4d.sh -v -gzip 2 \
                -i    $path/opendap/MOD04 \
                -j    $path/opendap/MYD04 \
                -o    $filename \
                -vars aod \
                -time $t1 $t1 \
                -func 'ave(const(@.1,0,-u)+const(@.2,0,1)/(if(@.1,==,-u,0,1)+if(@.2,==,-u,0,1)),time=$t1,time=$t2)'
               """
    else:
        tmpl = """
               lats4d.sh -v -gzip 2 \
                -i    $path/opendap/$prod \
                -o    $filename \
                -vars aod \
                -time $t1 $t1 \
                -func 'ave(aod,time=$t1,time=$t2)'
               """

    cmd = Template(tmpl).substitute(d)

    rc = os.system(cmd)

    if rc:
        raise RuntimeError, 'error on return from %s'%cmd

def xxx_m(path,filename,prod,gatime,inFile,var):
    """
    Given a time range in gatime, use lats4d to compute sampled mothly mean.

    path       top path for input control files, see below
    filename   output file name

    """

    if not Force:
        if os.path.exists(filename):
            print '<> File exists, skipping <%s>'%filename
            return

    d = dict(path=path, prod=prod, filename=filename, 
             inFile=inFile, var=var,
             t1=gatime[0], t2=gatime[1])

    tmpl = """
           lats4d.sh -v -gzip 2 \
                -i    $path/opendap/$inFile \
                -j    $path/opendap/$prod \
                -o    $filename \
                -vars $var \
                -time $t1 $t1 \
                -func 'ave(@+0*aod.2,time=$t1,time=$t2)'
          """

    cmd = Template(tmpl).substitute(d)

    rc = os.system(cmd)

    if rc:
        raise RuntimeError, 'error on return from %s'%cmd

#....................................................................
def wx_modulation(path,filename,prod,gatime):
    """
    Given a time range in gatime, use lats4d to compute the emission
    weather modulation factor monthly means. Sampled as the data.

    path       top path for input control files, see below
    filename   output file name

    """

    if not Force:
        if os.path.exists(filename):
            print '<> File exists, skipping <%s>'%filename
            return

    d = dict(path=path, prod=prod, filename=filename, 
             t1=gatime[0], t2=gatime[1])

    tmpl = """
           lats4d.sh -v -gzip 2 -udxt ginoux.udxt \
                -i    $path/opendap/w10m  \
                -j    $path/opendap/gwet  \
                -k    $path/opendap/$prod \
                -o    $filename \
                -vars w10m \
                -time $t1 $t1 \
                -func 'ave(ginoux5(w10m.1,gwettop.2)+0*aod.3,time=$t1,time=$t2)'
          """

    cmd = Template(tmpl).substitute(d)
    if os.system(cmd):
        raise RuntimeError, 'error on return from <%s>'%cmd

    cmd = "ncrename %s -v w10m,duwx"%filename
    if os.system(cmd):
        raise RuntimeError, 'error on return from <%s>'%cmd

#....................................................................

def nobs(path,filename,prod,gatime):
    """
    Given a time range in gatime, use lats4d to compute the MERRA-2
    dust emissions monthly means. Sampled as the data.

    path       top path for input control files, see below
    filename   output file name

    """

    if not Force:
        if os.path.exists(filename):
            print '<> File exists, skipping <%s>'%filename
            return

    d = dict(path=path, prod=prod, filename=filename, 
             t1=gatime[0], t2=gatime[1])

    tmpl = """
           lats4d.sh -v -gzip 2 \
                -i    $path/opendap/$prod  \
                -o    $filename \
                -vars aod \
                -time $t1 $t1 \
                -func 'ave(if(@,==,-u,0,1),time=$t1,time=$t2)'
          """

    cmd = Template(tmpl).substitute(d)
    if os.system(cmd):
        raise RuntimeError, 'error on return from <%s>'%cmd

    cmd = "ncrename %s -v aod,nobs"%filename
    if os.system(cmd):
        raise RuntimeError, 'error on return from <%s>'%cmd

#....................................................................

def foo(path,filename,prod,gatime):
    """
    Given a time range in gatime, use lats4d to compute the MERRA-2
    dust emissions monthly means. Sampled as the data.

    path       top path for input control files, see below
    filename   output file name

    """

    if not Force:
        if os.path.exists(filename):
            print '<> File exists, skipping <%s>'%filename
            return

    d = dict(path=path, prod=prod, filename=filename, 
             t1=gatime[0], t2=gatime[1])

    if prod=='MxD04':
        tmpl = """
               lats4d.sh -v -gzip 2 \
                -i    $path/opendap/MOD04  \
                -j    $path/opendap/MYD04  \
                -o    $filename \
                -vars aod \
                -time $t1 $t1 \
                -func 'ave(if(const(@.1,0,-u),>,0.2,1,0)+if(const(@.2,0,-u),>,0.2,1,0),time=$t1,time=$t2)'               """

    else:
        tmpl = """
               lats4d.sh -v -gzip 2 \
                -i    $path/opendap/$prod  \
                -o    $filename \
                -vars aod \
                -time $t1 $t1 \
                -func 'ave(if(const(@,0,-u),>,0.2,1,0),time=$t1,time=$t2)'
               """

    cmd = Template(tmpl).substitute(d)
    if os.system(cmd):
        raise RuntimeError, 'error on return from <%s>'%cmd

    cmd = "ncrename %s -v aod,foo"%filename
    if os.system(cmd):
        raise RuntimeError, 'error on return from <%s>'%cmd

#....................................................................

if __name__ == "__main__":

    path = '/nobackup/5/GAAS/DEE/'

    # Parse command line
    # ------------------
    if len(sys.argv) == 3:
        year1 = sys.argv[1]
        year2 = sys.argv[2]
    elif len(sys.argv) == 2 :
        year1 = sys.argv[1]
        year2 = year1
    else:
        print "   Usage:   %s  year1 [year2]"%sys.argv[0]
        print "Examples:   %s   2003  2005"%sys.argv[0]
        print "            %s   2003"%sys.argv[0]
        raise RuntimeError, 'not enough parameters'

    year1, year2 = int(year1), int(year2)    

    # Loop over time and products
    # ---------------------------
    for year in range(year1,year2+1):
        for month in range(1,13):
            #for prod in ('MOD04', 'MYD04'):
            for prod in ('MxD04',):

                gatime = mrange(year,month,gat=True)

                dirn = path+'/Level3/%s/Y%d/M%02d'%(prod,year,month)

                """
                # Deep Blue Mass flux
                # -------------------
                filename = '%s/dee_%s.uqvq.%d%02d.nc4'%(dirn,prod,year,month)
                mflux_o(path,filename,prod,gatime)

                # Deep Blue Removal
                # -----------------
                filename = '%s/dee_%s.durm.%d%02d.nc4'%(dirn,prod,year,month)
                removal(path,filename,prod,gatime)

                # Frequency of Occurence (based on gridbox mean)
                # ----------------------------------------------
                filename = '%s/dee_%s.foo.%d%02d.nc4'%(dirn,prod,year,month)
                foo(path,filename,prod,gatime)

                # Number of "obs"
                # ---------------
                filename = '%s/dee_%s.nobs.%d%02d.nc4'%(dirn,prod,year,month)
                nobs(path,filename,prod,gatime)

                """

                # MERRA-2 Mass flux
                # ------------------
                filename = '%s/dee_%s.uqvq_m.%d%02d.nc4'%(dirn,prod,year,month)
                mflux_m(path,filename,prod,gatime)

                # MERRA-2 AOD
                # -----------
                filename = '%s/dee_%s.aod_m.%d%02d.nc4'%(dirn,prod,year,month)
                inFile, var = 'du_cm', 'duexttau'
                xxx_m(path,filename,prod,gatime,inFile,var)

                # MERRA-2 CMASS
                # -------------
                filename = '%s/dee_%s.ducm_m.%d%02d.nc4'%(dirn,prod,year,month)
                inFile, var = 'du_cm', 'ducmass'
                xxx_m(path,filename,prod,gatime,inFile,var)

                # MERRA-2 Removal
                # ---------------
                filename = '%s/dee_%s.durm_m.%d%02d.nc4'%(dirn,prod,year,month)
                inFile, var = 'du_rm', 'durm'
                xxx_m(path,filename,prod,gatime,inFile,var)

                # MERRA-2 Emissions
                # -----------------
                filename = '%s/dee_%s.duem_m.%d%02d.nc4'%(dirn,prod,year,month)
                inFile, var = 'du_em', 'duem'
                xxx_m(path,filename,prod,gatime,inFile,var)

                # Weather Modulation
                # ------------------
                filename = '%s/dee_%s.duwx.%d%02d.nc4'%(dirn,prod,year,month)
                wx_modulation(path,filename,prod,gatime)

                # Observed AOD
                # ------------
                filename = '%s/dee_%s.aod_o.%d%02d.nc4'%(dirn,prod,year,month)
                aod_o(path,filename,prod,gatime)

