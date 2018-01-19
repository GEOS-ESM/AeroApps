#!/usr/bin/env python
# -W ignore::DeprecationWarning
""" 
Runscript for running geo_vlidort.x on NCCS


May 2017
patricia.castellanos@nasa.gov
"""

from   datetime           import datetime, timedelta 
from   dateutil.parser    import parse
import os
import shutil
import sys
import subprocess
from   distutils.dir_util import mkpath
import numpy              as np
import math
import time
import glob
import shutil
from   netCDF4           import Dataset
from   optparse          import OptionParser   # Command-line args
from   geo_vlidort_lc2   import JOBS, WORKSPACE

jobsmax   = 150
dt = timedelta(hours=1)
archive = '/archive/u/rgovinda/osse2/'


class CLD_WORKSPACE(WORKSPACE):
    """ Create Working Directories and RC files """
    def __init__(self,startdate,enddate,options):

        self.cloud = True
        for oo in options.__dict__:
            if (type(options.__dict__[oo]) is str) and (options.__dict__[oo].lower() == 'none'):
                self.__dict__[oo] = None
            else:
                self.__dict__[oo] = options.__dict__[oo]

        if self.nodemax is not None: self.nodemax = int(self.nodemax)
        self.nccs    = self.nccs + '/' + self.instname.upper() + '/CLD_DATA/' 
        self.prefix  = self.nccs + 'workdir/'

        self.indir   = self.nccs 
        self.outdir  = self.nccs + 'LevelC2'
        if (self.additional_output):
            self.addoutdir         = self.nccs + 'LevelC2'
        else:
            self.addoutdir         = None        

        # check for layout keyword. 
        # figure out number of tiles        
        if self.layout is None:
            self.ntiles = 1
        else:
            self.ntiles = int(self.layout[0])*int(self.layout[1])

        if type(self.channels) is str:
            if ',' in self.channels:
                self.channels = self.channels.replace(' ','').split(',')
            else:
                self.channels = self.channels.split()

        if type(self.c_band) is str:
            self.c_band.replace(',',' ')

        # Runmode code
        self.code = self.runmode + '.'
        if (self.interp.lower() == 'interpolate'):
            self.code += 'i'                 

        self.code += self.surface


        self.startdate = startdate
        self.enddate   = enddate

        ##################################################
        ####
        #    Loop through dates 
        #    Create working directories, SLURM scripts, and RC-files
        ####
        ##################################################

        # save run directory
        if self.verbose:
            print '++Saving run directory',os.getcwd()
        self.cwd     = os.getcwd()

        #initialize arrays to hold directory names
        self.dirstring    = []
        self.outdirstring = []
        self.nodemax_list = []
        if (self.additional_output):
            self.addoutdirstring = []

        # Loop over dates
        while (startdate <= enddate):
            # loop through tiles
            for tile in np.arange(self.ntiles):
                if self.layout is not None:
                    laycode = self.layout + str(tile)
                else:
                    laycode = None

                # check to see if there is any work to do
                # only simulating land pixels with limited SZAs and Cloud fractions
                numpixels = self.prefilter(startdate,layout=laycode) 
                
                if numpixels>0:

                    if (numpixels <= 1000 and self.nodemax is not None):
                        nodemax = 1
                    elif self.nodemax is not None:
                        nodemax = int(self.nodemax)
                    else:
                        nodemax = None

                    for i, ch in enumerate(self.channels):
                        # Create working directories for intermediate outputs
                        # create output directories
                        # save directory names
                        dirlist = self.make_workspace(startdate,ch,nodemax=nodemax,layout=laycode)

                        if (self.additional_output):
                            workdir, outdir, addoutdir_ = dirlist
                        else:
                            workdir, outdir = dirlist

                        self.dirstring.append(workdir)
                        self.outdirstring.append(outdir)
                        self.nodemax_list.append(nodemax)
                        if (self.additional_output):
                            self.addoutdirstring.append(addoutdir_)


                        # Create rcfiles - different for different surface types
                        if (self.surface.upper() == 'MAIACRTLS'):
                            self.make_maiac_rcfile(workdir,startdate,ch,nodemax=nodemax,
                                                   layout=laycode)
                        elif ('MCD43' in self.surface.upper()):
                            self.make_mcd43_rcfile(workdir,startdate,ch,nodemax=nodemax,
                                                   layout=laycode)                            
                        else:
                            self.make_ler_rcfile(workdir,startdate,ch,nodemax=nodemax,
                                                 layout=laycode)


            self.nodemax_list = np.array(self.nodemax_list)
            startdate = startdate + dt

    def prefilter(self,date,layout=None):
        if self.verbose:
            print '++Checking for good pixels in prefilter'
        g5dir = self.indir + '/LevelB/'+ 'Y'+ str(date.year) + '/M' + str(date.month).zfill(2) + '/D' + str(date.day).zfill(2) 
        nymd  = str(date.year) + str(date.month).zfill(2) + str(date.day).zfill(2)
        hour  = str(date.hour).zfill(2)

        if layout is None:
            geom  = g5dir + '/' + self.angname.lower() + '.cloud.lb2.angles.' + nymd + '_' + hour + 'z.nc4'
            land  = g5dir + '/' + self.instname.lower() + '-g5nr-icacl-TOTWPDF-GCOP-SKEWT.' + nymd + '_' + hour + 'z.nc4' 
            aer   = g5dir + '/' + self.instname.lower() + '-g5nr.lb2.aer_Nv.' + nymd + '_' + hour + 'z.nc4'
        else:
            geom  = g5dir + '/' + self.angname.lower() + '.cloud.lb2.angles.' + nymd + '_' + hour + 'z.' + layout +'.nc4'
            land  = g5dir + '/' + self.instname.lower() + '-g5nr-icacl-TOTWPDF-GCOP-SKEWT.' + nymd + '_' + hour + 'z.' + layout +'.nc4'
            aer   = g5dir + '/' + self.instname.lower() + '-g5nr.lb2.aer_Nv.' + nymd + '_' + hour + 'z.nc4'
 

        if self.verbose:
            print '++Opening geometry file ',geom
            print '++Opening land file', land

        if not os.path.exists(geom):
            self.get_from_archive(geom,date)
        ncGeom = Dataset(geom)
        SZA    = np.squeeze(ncGeom.variables[u'solar_zenith'][:])
        VZA    = np.squeeze(ncGeom.variables[u'sensor_zenith'][:])
        ncGeom.close()

        if not os.path.exists(land):
            self.get_from_archive(land,date)        
        ncLand = Dataset(land)
        FRLAND = np.squeeze(ncLand.variables[u'FRLAND'][:])
        ncLand.close()

        def clean_up(self,geom,land,aer):
            self.put_in_archive(geom)
            self.put_in_archive(land)
            self.put_in_archive(aer)

        f   = VZA < 80
        if np.sum(f) == 0:
            clean_up(self,geom,land,aer)
            return 0

        SZA    = SZA[f]
        FRLAND = FRLAND[f]
        f      = SZA < 80
        if np.sum(f) == 0:
            clean_up(self,geom,land,aer)
            return 0

        FRLAND = FRLAND[f]
        f      = FRLAND >= 0.99
        if np.sum(f) == 0:
            clean_up(self,geom,land,aer)
            return 0

        return np.sum(f)


    def destroy_workspace(self,i,jobid):
        # put LevelB files in archive or remove
        # --------------------

        #parse date from distring
        date = parse(os.path.basename(self.dirstring[i]).split('.')[1])

        g5dir = self.indir + '/LevelB/'+ 'Y'+ str(date.year) + '/M' + str(date.month).zfill(2) + '/D' + str(date.day).zfill(2) 
        nymd  = str(date.year) + str(date.month).zfill(2) + str(date.day).zfill(2)
        hour  = str(date.hour).zfill(2)

        if layout is None:
            geom  = g5dir + '/' + self.angname.lower() + '.cloud.lb2.angles.' + nymd + '_' + hour + 'z.nc4'
            land  = g5dir + '/' + self.instname.lower() + '-g5nr-icacl-TOTWPDF-GCOP-SKEWT.' + nymd + '_' + hour + 'z.nc4' 
            aer   = g5dir + '/' + self.instname.lower() + '-g5nr.lb2.aer_Nv.' + nymd + '_' + hour + 'z.nc4'
        else:
            geom  = g5dir + '/' + self.angname.lower() + '.cloud.lb2.angles.' + nymd + '_' + hour + 'z.' + layout +'.nc4'
            land  = g5dir + '/' + self.instname.lower() + '-g5nr-icacl-TOTWPDF-GCOP-SKEWT.' + nymd + '_' + hour + 'z.' + layout +'.nc4'
            aer   = g5dir + '/' + self.instname.lower() + '-g5nr.lb2.aer_Nv.' + nymd + '_' + hour + 'z.nc4'


        if self.archive_lb:
            self.put_in_archive(aer)
            self.put_in_archive(geom)
            self.put_in_archive(land)
            
        if not self.keep_lb:
            os.remove(aer)
            os.remove(geom)
            os.remove(land)


        os.chdir(self.dirstring[i])

        os.remove('Aod_EOS.rc')
        os.remove('Chem_MieRegistry.rc')
        os.remove(os.path.basename(self.execFile))
        os.remove('clean_mem.sh')
        os.remove('ExtData')
        os.remove('ExtDataCloud')
        os.remove('geo_vlidort.rc')
        os.remove(self.runfile)

        if self.nodemax is None:
            nodemax = None
        else:
            nodemax = self.nodemax_list[i]

        if self.additional_output:
            addoutdir = self.addoutdirstring[i]
        else:
            addoutdir = None

        outdir = self.outdirstring[i]

        if self.profile is False:
            if nodemax is not None and nodemax > 1:
                for a in np.arange(nodemax):
                    a = a + 1
                    errfile = 'slurm_' +jobid + '_' + str(a) + '.err'                    
                    outfile = 'slurm_' +jobid + '_' + str(a) + '.out'
                    if self.verbose:
                        print '++cleaning up errfile', errfile
                        print '++cleaning up outfile', outfile
                    os.remove(errfile)
                    os.remove(outfile)
                os.remove('slurm_%A_%a.out')

            else:
                errfile = 'slurm_' +jobid + '.err'
                os.remove(errfile)        
                outfile = 'slurm_' +jobid + '.out'
                os.remove(outfile)        

        def move_file(filelist,dest):
            for movefile in filelist:
                if os.path.exists(dest+'/'+movefile):
                    os.remove(dest+'/'+movefile)
                shutil.move(movefile,dest)

        #runscript to combine files
        if nodemax is not None and nodemax > 1:
            outfilelist = glob.glob('*.lc2.*.nc4')
            self.combine_files(outfilelist)
            outfilelist = glob.glob('*.lc2.*.nc4')
            move_file(outfilelist,outdir)
        else:
            outfilelist = glob.glob('*.lc2.*.nc4')
            move_file(outfilelist,outdir)


        if (addoutdir is not None):
            if nodemax is not None and nodemax > 1:
                outfilelist = glob.glob('*.add.*.nc4')
                self.combine_files(outfilelist)
                outfilelist = glob.glob('*.add.*.nc4')
                move_file(outfilelist,addoutdir)
            else:
                outfilelist = glob.glob('*.add.*.nc4')
                move_file(outfilelist,addoutdir)

        os.chdir(self.cwd)
        if self.profile is False:
            os.rmdir(self.dirstring[i])


    def make_maiac_rcfile(self,dirname,date,ch,nodemax=None,layout=None):
        os.chdir(dirname)

        rcfile = open('geo_vlidort.rc','w')
        rcfile.write('INDIR: '+self.indir+'\n')
        rcfile.write('OUTDIR: .\n')
        rcfile.write('DATE: '+str(date.year)+str(date.month).zfill(2)+str(date.day).zfill(2)+'\n')
        rcfile.write('TIME: '+str(date.hour).zfill(2)+'\n')
        rcfile.write('INSTNAME: ' + self.instname.lower() + '\n')
        rcfile.write('ANGNAME:' +self.angname+'\n')  
        rcfile.write('SURFNAME: MAIACRTLS\n')
        rcfile.write('SURFMODEL: RTLS\n')

        # Set Surface Configuration
        #figure out correct MODIS doy
        if (str(startdate.date()) == '2005-12-31'):
            rcfile.write('SURFDATE: 2006008\n')
        elif ((str(startdate.year) == '2007') and (str(date.month).zfill(2) == '03')):
            rcfile.write('SURFDATE: 2007072\n')
        elif ((str(startdate.year) == '2007') and (str(date.month).zfill(2) == '04')):
            rcfile.write('SURFDATE: 2007120\n')
        elif ((str(startdate.year) == '2006') and (str(date.month).zfill(2) == '07')):
            rcfile.write('SURFDATE: 2006216\n')
        elif ((str(startdate.year) == '2006') and (str(date.month).zfill(2) == '08')):
            rcfile.write('SURFDATE: 2006216\n')        
        else:
            doy = date.toordinal() - datetime(date.year-1,12,31).toordinal()
            DOY = 8*(int(doy/8) + 1)
            if (DOY > 365):
                DOY = 8
                rcfile.write('SURFDATE: '+str(date.year+1)+str(DOY).zfill(3)+'\n')
            else:
                rcfile.write('SURFDATE: '+str(date.year)+str(DOY).zfill(3)+'\n')

        if self.c_band is None:
            self.c_band = "645 858 469 555 1240 1640 2130 412"

        rcfile.write('SURFBANDM: '+ str(len(self.c_band.split(' '))) + '\n')

        rcfile.write('SURFBAND: ' + self.interp +'\n')
        if (self.interp.upper() == 'INTERPOLATE'):
            rcfile.write('SURFBAND_C: ' + self.c_band + '\n')
        else:
            #figure out index
            c_band = np.array(self.c_band.split(' ')).astype('float')
            if ch >= c_band.max():
                i_band = np.argmax(c_band)
            elif ch <= c_band.min():
                i_band = np.argmin(c_band)
            else:
                i_band = np.argmin(np.abs((ch-c_band) ))
            rcfile.write('SURFBAND_I: '+ str(i_band)  + '\n')

        if (self.code == 'scalar'):
            rcfile.write('SCALAR: true\n')
        else:
            rcfile.write('SCALAR: false\n')


        rcfile.write('CHANNELS: '+ch+'\n')
        if nodemax is not None:
            rcfile.write('NODEMAX: '+ str(nodemax) + '\n')

        if self.additional_output:
            rcfile.write('ADDITIONAL_OUTPUT: true\n')
        else:
            rcfile.write('ADDITIONAL_OUTPUT: false\n')

        if self.version is not None:
            rcfile.write('VERSION: '+self.version+'\n')

        if self.surf_version is not None:
            rcfile.write('SURF_VERSION: '+self.surf_version+'\n')

        if layout is not None:
            rcfile.write('LAYOUT: '+layout+'\n')

        # Set up cloud configuration
        rcfile.write('ICLDTABLE:'+self.icldtable+'\n') 
        rcfile.write('LCLDTABLE:'+self.lcldtable+'\n') 
        cld_band = np.array(self.cld_band.split(',')).astype('float')
        if ch >= cld_band.max():
            idxcld = np.argmax(cld_band)
        elif ch <= cld_band.min():
            idxcld = np.argmin(cld_band)
        else:
            idxcld = np.argmin(np.abs((ch-cld_band) ))

        rcfile.write('IDXCLD:'+ str(idxcld) +'\n')   
            
        rcfile.close()

        os.chdir(self.cwd)


    def make_ler_rcfile(self,dirname,date,ch,nodemax=None,layout=None):
        os.chdir(dirname)

        rcfile = open('geo_vlidort.rc','w')
        rcfile.write('INDIR: ' + self.indir + '\n')
        rcfile.write('OUTDIR: .\n')
        rcfile.write('DATE: ' + str(date.year) + str(date.month).zfill(2) + str(date.day).zfill(2) + '\n')
        rcfile.write('TIME: ' + str(date.hour).zfill(2) + '\n')
        rcfile.write('INSTNAME: ' + self.instname + '\n')
        rcfile.write('ANGNAME:' +self.angname+'\n')  
        rcfile.write('SURFNAME: LER\n')
        rcfile.write('SURFMODEL: LER\n')

        # Set Surface Configuration
        rcfile.write('SURFDATE: ' + str(date.month).zfill(2) +'\n')

        if self.c_band is None:
            self.c_band = "354 388"

        rcfile.write('SURFBANDM: '+ str(len(self.c_band.split(' '))) + '\n')

        rcfile.write('SURFBAND: ' + self.interp + '\n')
        if (self.interp.upper() == 'INTERPOLATE'):
            rcfile.write('SURFBAND_C: ' + self.c_band + '\n')
        else:
            #figure out index
            c_band = np.array(self.c_band.split(' ')).astype('float')
            if ch >= c_band.max():
                i_band = np.argmax(c_band)
            elif ch <= c_band.min():
                i_band = np.argmin(c_band)
            else:
                i_band = np.argmin(np.abs((ch-c_band) ))
            rcfile.write('SURFBAND_I: '+ str(i_band) + '\n')

        if (self.code == 'scalar'):
            rcfile.write('SCALAR: true\n')
        else:
            rcfile.write('SCALAR: false\n')


        rcfile.write('CHANNELS: ' + ch + '\n')
        if nodemax is not None:
            rcfile.write('NODEMAX: '+ str(nodemax) + '\n')

        if self.additional_output:
            rcfile.write('ADDITIONAL_OUTPUT: true\n')
        else:
            rcfile.write('ADDITIONAL_OUTPUT: false\n')

        if self.version is not None:
            rcfile.write('VERSION: '+self.version+'\n')

        if self.surf_version is not None:
            rcfile.write('SURF_VERSION: '+self.surf_version+'\n')

        if layout is not None:
            rcfile.write('LAYOUT: '+layout+'\n')

        rcfile.write('ICLDTABLE:'+self.icldtable+'\n') 
        rcfile.write('LCLDTABLE:'+self.lcldtable+'\n') 
        cld_band = np.array(self.cld_band.split(',')).astype('float')
        if ch >= cld_band.max():
            idxcld = np.argmax(cld_band)
        elif ch <= cld_band.min():
            idxcld = np.argmin(cld_band)
        else:
            idxcld = np.argmin(np.abs((ch-cld_band) ))

        rcfile.write('IDXCLD:'+ str(idxcld) +'\n')   
        rcfile.close()

        os.chdir(self.cwd)    


    def make_mcd43_rcfile(self,dirname,date,ch,nodemax=None,layout=None):
        os.chdir(dirname)

        rcfile = open('geo_vlidort.rc','w')
        rcfile.write('INDIR: '+self.indir+'\n')
        rcfile.write('OUTDIR: .\n')
        rcfile.write('DATE: '+str(date.year)+str(date.month).zfill(2)+str(date.day).zfill(2)+'\n')
        rcfile.write('TIME: '+str(date.hour).zfill(2)+'\n')
        rcfile.write('INSTNAME: ' + self.instname.lower() + '\n')
        rcfile.write('ANGNAME:' +self.angname+'\n')  
        rcfile.write('SURFNAME: '+self.surface+'\n')
        rcfile.write('SURFMODEL: RTLS\n')

        # Set Surface Configuration
        doy = date.toordinal() - datetime(date.year-1,12,31).toordinal()
        rcfile.write('SURFDATE: '+str(date.year)+str(doy).zfill(3)+'\n')

        if self.c_band is None:
            self.c_band = "645 858 469 555 1240 1640 2130"

        rcfile.write('SURFBANDM: '+ str(len(self.c_band.split(' '))) + '\n')

        rcfile.write('SURFBAND: ' + self.interp +'\n')
        if (self.interp.upper() == 'INTERPOLATE'):
            rcfile.write('SURFBAND_C: ' + self.c_band + '\n')
        else:
            #figure out index
            c_band = np.array(self.c_band.split(' ')).astype('float')
            if ch >= c_band.max():
                i_band = np.argmax(c_band)
            elif ch <= c_band.min():
                i_band = np.argmin(c_band)
            else:
                i_band = np.argmin(np.abs((ch-c_band) ))
            rcfile.write('SURFBAND_I: '+ str(i_band) + '\n')

        if (self.code == 'scalar'):
            rcfile.write('SCALAR: true\n')
        else:
            rcfile.write('SCALAR: false\n')


        rcfile.write('CHANNELS: '+ch+'\n')
        if nodemax is not None:
            rcfile.write('NODEMAX: '+ str(nodemax) + '\n')

        if self.additional_output:
            rcfile.write('ADDITIONAL_OUTPUT: true\n')
        else:
            rcfile.write('ADDITIONAL_OUTPUT: false\n')

        if self.version is not None:
            rcfile.write('VERSION: '+self.version+'\n')

        if self.surf_version is not None:
            rcfile.write('SURF_VERSION: '+self.surf_version+'\n')

        if layout is not None:
            rcfile.write('LAYOUT: '+layout+'\n')

        rcfile.write('ICLDTABLE:'+self.icldtable+'\n') 
        rcfile.write('LCLDTABLE:'+self.lcldtable+'\n') 
        cld_band = np.array(self.cld_band.split(',')).astype('float')
        if ch >= cld_band.max():
            idxcld = np.argmax(cld_band)
        elif ch <= cld_band.min():
            idxcld = np.argmin(cld_band)
        else:
            idxcld = np.argmin(np.abs((ch-cld_band) ))

        rcfile.write('IDXCLD:'+ str(idxcld) +'\n')               
        rcfile.close()

        os.chdir(self.cwd)



        
#########################################################

if __name__ == "__main__":
    # Defaults
    # ------------
    instname          = 'tempo'
    version           = '1.0'    
    channels          = '388'
    surface           = 'lambertian' #'lambertian' or 'MAIACRTLS'
    interp            = 'exact'  #'interpolate' or 'exact'
    
    nodemax           = 10
    surf_version      = '1.0'
    layout            = None
    CLDMAX            = '0.01'
    nccs              = '/discover/nobackup/projects/gmao/osse2/pub/c1440_NR/OBS/'
    
    runmode           = 'vector'
    runfile           = 'geo_vlidort_lc2.j'
    execFile          = '/discover/nobackup/pcastell/workspace/GAAS/src/Components/geo_vlidort/geo_vlidort_cloud.x'
    
    icldtable         = 'ExtDataCloud/IceLegendreCoeffs.nc4'
    lcldtable         = 'ExtDataCloud/WaterLegendreCoeffs.nc4'
    cld_band          = '650, 860, 470, 550, 1240, 1630, 2130, ' + \
                        '410, 440, 910, 936, 940, 3700, 3900, ' + \
                        '380, 6200, 7300, 8500, 11000, 12000, ' + \
                        '13200, 13400, 13800, 14200'
    cld_band          = cld_band.replace(' ','')
    #Flags
    # verbose           = False
    # additional_output = False
    # profile           = False
    # keep_lb           = False
    # archive_lb        = False

    # Parse command line options
    # ------------------------------
    parser = OptionParser(usage="Usage: %prog [options] startdate enddate")


    parser.add_option("-I", "--instname", dest="instname", default=instname,
                      help="Instrument name (default=%s)"\
                      %instname )

    parser.add_option("--angname", dest="angname", default=instname,
                      help="Instrument name for angles (default=%s)"\
                      %instname )    

    parser.add_option("-V", "--version_string", dest="version", default=version,
                      help="Version name (default=%s)"\
                      %version )        

    parser.add_option("-c", "--channels", dest="channels", default=channels,
                      help="Channels (default=%s)"\
                      %channels )  

    parser.add_option("--icldtable", dest="icldtable", default=icldtable,
                      help="Ice cloud optics table (default=%s)"\
                      %icldtable )     

    parser.add_option("--lcldtable", dest="lcldtable", default=lcldtable,
                      help="Liquid cloud optics table (default=%s)"\
                      %lcldtable ) 

    parser.add_option("--cld_band", dest="cld_band", default=cld_band,
                      help="Bands in cloud optics table(default=%s)"\
                      %cld_band )                                                 

    parser.add_option("-s", "--surface", dest="surface", default=surface,
                      help="Surface Reflectance Dataset.  Choose from 'lambertian' or 'MAIACRTLS' or 'MCD43X' "\
                      "(default=%s)"\
                      %surface )      

    parser.add_option("-S", "--surfversion", dest="surf_version", 
                      help="Surface Reflectance Dataset version.  Required if surface==MAIACRTLS' " )          

    parser.add_option("-i", "--interp", dest="interp", default=interp,
                      help="Surface reflectance channel interpolation. Choose from 'interpolate' or 'exact' "\
                      "(default=%s)"\
                      %interp )      

    parser.add_option("-b", "--c_band", dest="c_band",
                      help="Surface reflectance bands." )  

    parser.add_option("-a", "--additional",
                      action="store_true", dest="additional_output",default=False,
                      help="Turn on writing additional_output. (default=False)")                           

    parser.add_option("-v", "--verbose",action="store_true",
                      dest="verbose", default=False,
                      help="Verbose (default=False)" )    

    parser.add_option("--keep_lb",action="store_true",
                      dest="keep_lb", default=False,
                      help="keep LevelB files - do not remove after geo_vlidort is finished (default=False)" )    

    parser.add_option("--archive_lb",action="store_true",
                      dest="archive_lb", default=False,
                      help="archive LevelB files - copy to archive dir after geo_vlidort is finished (default=False)" )    

    parser.add_option("-n", "--nodemax", dest="nodemax", default=nodemax,
                      help="Max number of nodes to use. "\
                      "(default=%s)"\
                      %nodemax )    

    parser.add_option("-l", "--layout", dest="layout", default=layout,
                      help="Layout of domain. Used for breaking up high-res domains (e.g. GOES-R)"\
                      "(default=%s)"\
                      %layout )  

    parser.add_option("-d", "--dir", dest="nccs", default=nccs,
                      help="Root level directory for inputs/outputs "\
                      "(default=%s)"\
                      %nccs )                                 

    parser.add_option("-p", "--profile", 
                      action="store_true",dest="profile", default=False,
                      help="Leave slurm output files in working directory. (default=False)")   

    parser.add_option("-r", "--runmode", dest="runmode", default=runmode,
                      help="VLIDORT run mode. Either 'scalar' or 'vector' "\
                      "(default=%s)"\
                      %runmode )       

    parser.add_option("-f", "--runfile", dest="runfile", default=runfile,
                      help="slurm script template "\
                      "(default=%s)"\
                      %runfile )       

    parser.add_option("-e", "--execfile", dest="execFile", default=execFile,
                      help="geo_vlidort executable "\
                      "(default=%s)"\
                      %execFile )     

    parser.add_option("-A", "--archive", dest="archive", default=archive,
                      help="where to look for missing data on archive"\
                      "(default=%s)"\
                      %archive )                                                                            


    ################
    ###
    #    End of uper inputs
    ###
    ################    
    (options, args) = parser.parse_args()

    if len(args) == 2:
        startdate, enddate = args
    else:
        parser.error("must have 2 arguments: startdate and enddate")    


    # Setup workspace and runscript
    # ------------------------------
    startdate = parse(startdate)
    enddate   = parse(enddate)
    workspace = CLD_WORKSPACE(startdate,enddate,options)

    # Submit and Handle Jobs
    # -----------------------------
    if (workspace.dirstring) > 0:
        workspace.handle_jobs()
    else:
        print 'No model hours to run'

    


