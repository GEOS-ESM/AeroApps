#!/bin/python
# -*- coding: utf-8 -*-
""" Runscript for geo_vlidort episodes"""
from   datetime           import datetime, timedelta 
from   dateutil.parser    import parse
import os
import subprocess
from   distutils.dir_util import mkpath
import numpy              as np
import math
import time
import glob
import shutil
from   netCDF4           import Dataset
from   optparse        import OptionParser   # Command-line args

jobsmax   = 150
dt = timedelta(hours=1)

class JOBS(object):
    def handle_jobs(self):
        # Figure out how many jobs you need to submit
        runlen  = len(self.dirstring)   

        if self.nodemax is not None:
            numjobs = sum(self.nodemax_list)
        else:
            numjobs = runlen
        
        devnull = open(os.devnull, 'w')
        if numjobs <= jobsmax:   
            countRun   = runlen    
            node_tally   = numjobs  
        else:
            if self.nodemax is not None:
                keep_adding = True
                node_tally = 0
                countRun = 0
                while(keep_adding):
                    if (node_tally + self.nodemax_list[countRun] <= jobsmax):
                        node_tally = node_tally + self.nodemax_list[countRun]
                        countRun = countRun + 1
                    else:
                        keep_adding = False                
            else:
                countRun = jobsmax
                node_tally = jobsmax

        workingJobs = np.arange(countRun)


        # Submit first JOBSMAX jobs
        jobid = np.empty(0)
        for i in workingJobs:
            s = self.dirstring[i]
            os.chdir(s)
            jobid = np.append(jobid,subprocess.check_output(['qsub',self.runfile]))
            
        os.chdir(self.cwd)


        # Monitor jobs 1-by-1 
        # Add a new job when one finishes 
        # Until they are all done
        stat = subprocess.call(['qstat -u pcastell'], shell=True, stdout=devnull)
        while (stat == 0):
            stat = subprocess.call(['qstat -u pcastell'], shell=True, stdout=devnull)
            finishedJobs = np.empty(0,dtype='int')
            for ii,i in enumerate(workingJobs):
                s = jobid[i]
                s = s.strip('\n')
                finished = False

                # Check to see if this job is finished
                if self.nodemax is not None and self.nodemax_list[i] > 1:
                    # Loop through the nodes working on this job
                    finishedCNT = 0
                    for a in np.arange(nodemax_list[i]):
                        a = a + 1

                        try:
                            result = subprocess.check_output(['squeue','-j',s+'_'+ str(a)],stderr=subprocess.STDOUT)
                            if (s+'_'+ str(a) not in result):
                                #print ' Job Not Found! '+s+'_'+ str(a)
                                finishedCNT = finishedCNT + 1
                        except subprocess.CalledProcessError as e:
                            #print ' ERROR Job Not Found! '+s+'_'+ str(a)
                            finishedCNT = finishedCNT + 1

                    if (finishedCNT == self.nodemax_list[i]): 
                        finished = True
                else:
                    result = subprocess.call(['qstat',s], stdout=devnull)
                    if (result != 0):
                        finished = True   

                # if the job is finished clean up the workspace
                if finished:
                    #print 'Job finished, cleaning up', s, i 
                    finishedJobs = np.append(finishedJobs,ii)
                    errcheck = self.check_for_errors(i,s)               
                    if (errcheck is False):
                        self.destroy_workspace(i,s)
                    else:
                        print 'Jobid ',s,' exited with errors'

            # finished checking up on all the jobs
            # Remove finished jobs from the currently working list
            if len(finishedJobs) != 0:
                print 'deleting finishedJobs',finishedJobs,jobid[workingJobs[finishedJobs]]
                if self.nodemax is not None:
                    node_tally  = node_tally - sum(self.nodemax_list[workingJobs[finishedJobs]])
                else:
                    node_tally  = node_tally - len(finishedJobs)

                workingJobs = np.delete(workingJobs,finishedJobs)

            # Add more jobs if needed
            # reinitialize stat variable
            if (runlen > countRun) and (node_tally < jobsmax):
                #print 'adding new jobs'
                if self.nodemax is not None:
                    keep_adding = True
                    newRun      = 0
                    while (keep_adding):
                        if (node_tally + self.nodemax_list[countRun + newRun] <= jobsmax):
                            node_tally = node_tally + self.nodemax_list[countRun + newRun]
                            newRun = newRun + 1
                            if (countRun + newRun == runlen):
                                keep_adding = False
                        else:
                            keep_adding = False
                else:
                    newRun = jobsmax - node_tally
                    node_tally = jobsmax

                newjobs  = countRun + np.arange(newRun)
                workingJobs = np.append(workingJobs, newjobs)
                for i in newjobs:
                    s = dirstring[i]
                    os.chdir(s)
                    jobid = np.append(jobid,subprocess.check_output(['qsub',runfile]))

                os.chdir(cwd)
                countRun = countRun + newRun
                stat = subprocess.call(['qstat -u pcastell'], shell=True, stdout=devnull)


            print 'Waiting 1 minutes'
            time.sleep(60)
            

        # Exited while loop
        print 'All jobs done'


    def check_for_errors(self,i,jobid):
        os.chdir(self.dirstring[i])  

        error = False  
        if self.nodemax is None:
            nodemax = None
        else:
            nodemax = self.nodemax_list[i]


        if self.nodemax is not None and nodemax > 1:
            for a in np.arange(nodemax):
                a = a + 1
                errfile = 'slurm_' +jobid + '_' + str(a) + '.err'
                statinfo = os.stat(errfile)
                if (statinfo.st_size != 0):
                    error = True
        else:
            errfile = 'slurm_' +jobid + '.err'
            statinfo = os.stat(errfile)
            if (statinfo.st_size != 0):
                error = True

        os.chdir(self.cwd)
        return error


class WORKSPACE(JOBS):
    """ Create Working Directories and RC files """
    def __init__(self,startdate,enddate,options):


        for oo in options.__dict__:
            self.__dict__[oo] = options.__dict__[oo]

        self.runfile = 'geo_vlidort_run_array.j'
        self.nccs    = self.nccs + self.instname.upper() + '/DATA/' 
        self.prefix  = self.nccs + 'workdir/'

        self.indir   = nccs 
        self.outdir  = nccs + 'LevelC2'
        if (self.additional_output):
            self.addoutdir         = self.nccs + 'LevelC2'
        else:
            self.addoutdir         = None        

        # check for layout keyword. 
        # figure out number of tiles        
        if self.layout is None:
            self.ntiles = 1
        else:
            self.ntiles = int(layout[0])*int(layout[1])

        if type(self.channels) is str:
            self.channels = self.channels.split()

        if type(self.i_band) is str:
            self.i_band = self.i_band.split()

        # Runmode code
        self.code = self.runmode + '.'
        if (self.interp.lower() == 'interpolate'):
            self.code += 'i'                 

        self.code += surface


        self.startdate = startdate
        self.enddate   = enddate

        ##################################################
        ####
        #    Loop through dates 
        #    Create working directories, SLURM scripts, and RC-files
        ####
        ##################################################

        # save run directory
        self.rundir     = os.getcwd()

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
                    else:
                        nodemax = self.nodemax

                    for i, ch in enumerate(self.channels):
                        try:
                            i_band = self.i_band[i]
                        except:
                            i_band = None

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
                        self.nodemax_list.append(int(nodemax))
                        if (self.additional_output):
                            self.addoutdirstring.append(addoutdir)


                        # Create rcfiles - different for different surface types
                        if (self.surface.upper() == 'MAIACRTLS'):
                            self.make_maiac_rcfile(workdir,startdate,ch,nodemax=nodemax,i_band=i_band,
                                                   layout=laycode)
                        else:
                            self.make_ler_rcfile(workdir,startdate,ch,nodemax=nodemax,i_band=i_band,
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
            met   = g5dir + '/' + self.instname.lower() + '-g5nr.lb2.met_Nv.' + nymd + '_' + hour + 'z.nc4'
            geom  = g5dir + '/' + self.instname.lower() + '.lb2.angles.' + nymd + '_' + hour + 'z.nc4'
            land  = self.indir + '/LevelB/invariant/' + self.instname.lower() + '-g5nr.lb2.asm_Nx.nc4'  
        else:
            met   = g5dir + '/' + self.instname.lower() + '-g5nr.lb2.met_Nv.' + nymd + '_' + hour + 'z_' + layout +'.nc4'
            geom  = g5dir + '/' + self.instname.lower() + '.lb2.angles.' + nymd + '_' + hour + 'z_' + layout +'.nc4'
            land  = self.indir + '/LevelB/invariant/' + self.instname.lower() + '-g5nr.lb2.asm_Nx_' + layout + '.nc4'  

        if self.verbose:
            print '++Opening metfile ',met
        ncMet = Dataset(met)
        Cld   = np.squeeze(ncMet.variables[u'CLDTOT'][:])
        f     = np.where(Cld <= float(self.CLDMAX))
        ncMet.close()
        if len(f[0]) == 0:
            return False, 0

        ncGeom = Dataset(geom)
        SZA    = np.squeeze(ncGeom.variables[u'solar_zenith'][:])
        VZA    = np.squeeze(ncGeom.variables[u'sensor_zenith'][:])
        ncGeom.close()
        ncLand = Dataset(land)
        FRLAND = np.squeeze(ncLand.variables[u'FRLAND'][:])
        ncLand.close()

        SZA = SZA[f]
        VZA = VZA[f]
        FRLAND = FRLAND[f]

        f   = np.where(VZA < 80)
        if len(f[0]) == 0:
            return 0

        SZA    = SZA[f]
        FRLAND = FRLAND[f]
        f      = np.where(SZA < 80)
        if len(f[0]) == 0:
            return 0

        FRLAND = FRLAND[f]
        f      = np.where(FRLAND >= 0.99)
        if len(f[0]) == 0:
            return 0

        return len(f[0])


    def make_workspace(date,ch,nodemax=None,layout=None):

        # creating working directory
        dirname = '{}/{}.{}T{}.{}'.format(self.prefix,self.instname.lower(),date.date(),str(date.hour).zfill(2),ch,self.code)
        jobname = '{}.{}T{}.{}'.format(instname.lower(),date.date(),str(date.hour).zfill(2),ch)

        if layout is not None:
            dirname = dirname + '.' + layout
            jobname = jobname + '.' + layout

        bindir = dirname
        if not os.path.exists(dirname):
            os.makedirs(dirname)            

        # outdirname
        outdir = '{}/Y{}/M{}/D{}'.format(self.outdir,date.year,str(date.month).zfill(2),str(date.day).zfill(2))
        if self.addoutdir is not None:
            addoutdir = '{}/Y{}/M{}/D{}'.format(self.addoutdir,date.year,str(date.month).zfill(2),str(date.day).zfill(2))

        # create some symlinks to shared files
        os.chdir(dirname)
        if not os.path.isfile('geo_vlidort.x'):
            os.symlink(self.cwd+'/geo_vlidort.x','geo_vlidort.x')
        if not os.path.exists('ExtData'):
            os.symlink(self.cwd+'/ExtData','ExtData')
        if not os.path.isfile('Chem_MieRegistry.rc'):
            os.symlink(self.cwd+'/Chem_MieRegistry.rc','Chem_MieRegistry.rc')
        if not os.path.isfile('clean_mem.sh'):
            os.symlink(self.cwd+'/clean_mem.sh','clean_mem.sh')
            
        # Create slurm runfile - using self.runfile as template           
        source = open(self.cwd+'/'+self.runfile,'r')
        destination = open(self.runfile,'w')
        for line in source:
            if (line[0:18] == '#SBATCH --job-name'):            
                destination.write('#SBATCH --job-name='+jobname+'\n')

            elif (line[0:16] == '#SBATCH --output'):
                if (nodemax is not None and nodemax > 1):
                    destination.write('#SBATCH --output='+'slurm'+"_%A_%a.out"+'\n')
                else:
                    destination.write('#SBATCH --output='+'slurm'+"_%j.out"+'\n')

            elif (line[0:15] == '#SBATCH --error'):
                if (nodemax is not None and nodemax > 1):
                    destination.write('#SBATCH --error='+'slurm'+"_%A_%a.err"+'\n')      
                else:
                    destination.write('#SBATCH --error='+'slurm'+"_%j.err"+'\n')

            elif (line[0:15] == '#SBATCH --array'):
                if (nodemax is not None and nodemax > 1):
                    destination.write('#SBATCH --array=1-'+str(nodemax)+'\n') 

            elif (line[0:13] == 'setenv GEOBIN'):
                destination.write('setenv GEOBIN '+bindir+'\n')

            elif (line[0:13] == 'setenv OUTDIR'):            
                destination.write('setenv OUTDIR '+outdir+'\n')            
                if not os.path.exists(outdir):
                    mkpath(outdir)

            elif (line[0:16] == 'setenv ADDOUTDIR'):   
                if (self.addoutdir is not None):       
                    destination.write('setenv ADDOUTDIR '+addoutdir+'\n')                          
                    if not os.path.exists(addoutdir):
                        mkpath(addoutdir)
                else:
                    destination.write('setenv ADDOUTDIR None \n') 

            elif (line[0:8] == '$RUN_CMD'):
                if (nodemax is not None and nodemax > 1):     
                    destination.write('$RUN_CMD  ./geo_vlidort.x geo_vlidort.rc ${SLURM_ARRAY_TASK_ID}'+'\n')    
                else:
                    destination.write('$RUN_CMD  ./geo_vlidort.x geo_vlidort.rc'+'\n')   

            else:
                destination.write(line)        
        source.close()
        destination.close()
        
        # Change wavelength number in Aod_EOS.rc
        source = open(self.cwd+'/Aod_EOS.rc','r')
        destination = open('Aod_EOS.rc','w')
        a = float(ch)*1e-3
        for line in source:
            if (line[0:11] == 'r_channels:'):
                destination.write('r_channels: '+'{:0.3f}e-6'.format(a)+'\n')
            else:
                destination.write(line) 
        source.close()
        destination.close()
        
        # Go back to original run directory 
        os.chdir(self.cwd) 
        
        if (self.addoutdir is not None):
            return dirname, outdir, addoutdir
        else:
            return dirname, outdir 

    def make_maiac_rcfile(dirname,date,ch,nodemax=None,i_band=None,layout=None):
        os.chdir(dirname)

        rcfile = open('geo_vlidort.rc','w')
        rcfile.write('INDIR: '+self.indir+'\n')
        rcfile.write('OUTDIR: .\n')
        rcfile.write('DATE: '+str(date.year)+str(date.month).zfill(2)+str(date.day).zfill(2)+'\n')
        rcfile.write('TIME: '+str(date.hour).zfill(2)+'\n')
        rcfile.write('INSTNAME: ' + self.instname.lower() + '\n')
        rcfile.write('SURFNAME: MAIACRTLS\n')
        rcfile.write('SURFMODEL: RTLS\n')

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


        rcfile.write('SURFBAND: ' + self.interp +'\n')
        if (self.interp.upper() == 'INTERPOLATE'):
            rcfile.write('SURFBAND_C: 645 858 469 555 1240 1640 2130 412\n')
        else:
            rcfile.write('SURFBAND_I: '+ i_band + '\n')

        rcfile.write('SURFBANDM: 8 \n')

        if (self.code == 'scalar'):
            rcfile.write('SCALAR: true\n')
        else:
            rcfile.write('SCALAR: false\n')


        rcfile.write('CHANNELS: '+ch+'\n')
        rcfile.write('CLDMAX: '+self.CLDMAX+'\n')
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
            
        rcfile.close()

        os.chdir(self.cwd)


    def make_ler_rcfile(dirname,date,ch,nodemax=None,i_band=None,layout=None):
        os.chdir(dirname)

        rcfile = open('geo_vlidort.rc','w')
        rcfile.write('INDIR: ' + self.indir + '\n')
        rcfile.write('OUTDIR: .\n')
        rcfile.write('DATE: ' + str(date.year) + str(date.month).zfill(2) + str(date.day).zfill(2) + '\n')
        rcfile.write('TIME: ' + str(date.hour).zfill(2) + '\n')
        rcfile.write('INSTNAME: ' + self.instname + '\n')
        rcfile.write('SURFNAME: LER\n')
        rcfile.write('SURFMODEL: LER\n')

        rcfile.write('SURFDATE: ' + str(date.month).zfill(2) +'\n')

        rcfile.write('SURFBAND: ' + self.interp + '\n')
        if (self.interp.upper() == 'INTERPOLATE'):
            rcfile.write('SURFBAND_C: 354 388\n')
        else:
            rcfile.write('SURFBAND_I: ' + i_band + '\n')

        rcfile.write('SURFBANDM: 2 \n')

        if (self.code == 'scalar'):
            rcfile.write('SCALAR: true\n')
        else:
            rcfile.write('SCALAR: false\n')


        rcfile.write('CHANNELS: ' + ch + '\n')
        rcfile.write('CLDMAX: '+self.CLDMAX+' \n')
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

        rcfile.close()

        os.chdir(self.cwd)    



    def destroy_workspace(self,i,jobid):
        os.chdir(self.dirstring[i])

        os.remove('Aod_EOS.rc')
        os.remove('Chem_MieRegistry.rc')
        os.remove('geo_vlidort.x')
        os.remove('clean_mem.sh')
        os.remove('ExtData')
        os.remove('geo_vlidort.rc')
        os.remove('geo_vlidort_run_array.j')

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
                    os.remove(errfile)
                    outfile = 'slurm_' +jobid + '_' + str(a) + '.out'
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
                combine_files(outfilelist)
                outfilelist = glob.glob('*.add.*.nc4')
                move_file(outfilelist,addoutdir)
            else:
                outfilelist = glob.glob('*.add.*.nc4')
                move_file(outfilelist,addoutdir)

        os.chdir(self.cwd)
        if profile is False:
            os.rmdir(self.dirstring[i])

    def combine_files(self,filelist):
        mergedfile = filelist[0]
        parts      = mergedfile.split('.')
        mergedfile = ''
        for p in parts[0:-2]:
            mergedfile = mergedfile + p + '.'
        mergedfile = mergedfile + 'nc4'
        #mergedfile = mergedfile[:-5] + 'nc4'

        os.rename(filelist[0],mergedfile)
        filelist = filelist[1:]

        if type(filelist) is str:
            filelist = filelist.split()

        ncmergedfile = Dataset(mergedfile, mode='r+')
        for filename in filelist:
            ncfile = Dataset(filename)
            for var in ncmergedfile.variables:
                vardata = ncfile.variables[var][:]
                mergedata = ncmergedfile.variables[var][:]
                if (np.ma.is_masked(vardata)):
                    mergedata[~vardata.mask] = vardata[~vardata.mask]
                    ncmergedfile.variables[var][:] = mergedata
            ncfile.close()
            os.remove(filename)

        ncmergedfile.close()






        
#########################################################

if __name__ == "__main__":
    # Defaults
    # ------------
    instname          = 'tempo'
    version           = '1.0'    
    channels          = '388'
    surface           = 'lambertian' #'lambertian' or 'MAIACRTLS'
    interp            = 'exact'  #'interpolate' or 'exact'
    i_band            = '2'    
    
    nodemax           = 10
    surf_version      = '1.0'
    layout            = None
    CLDMAX            = '0.01'
    nccs              = '/discover/nobackup/projects/gmao/osse2/pub/c1440_NR/OBS/'

    profile           = False
    runmode           = 'vector'
    
    #Flags
    # verbose           = False
    # additional_output = False

    # Parse command line options
    # ------------------------------
    parser = OptionParser(usage="Usage: %prog [options] startdate enddate")


    parser.add_option("-I", "--instname", dest="instname", default=instname,
                      help="Instrument name (default=%s)"\
                      %instname )

    parser.add_option("-V", "--version_string", dest="version", default=version,
                      help="Version name (default=%s)"\
                      %version )        

    parser.add_option("-c", "--channels", dest="channels", default=channels,
                      help="Channels (default=%s)"\
                      %channels )  

    parser.add_option("-s", "--surface", dest="surface", default=surface,
                      help="Surface Reflectance Dataset.  Choose from 'lambertian' or 'MAIACRTLS' "\
                      "(default=%s)"\
                      %surface )      

    parser.add_option("-S", "--surfversion", dest="surf_version", 
                      help="Surface Reflectance Dataset version.  Required if surface==MAIACRTLS' " )          

    parser.add_option("-i", "--interp", dest="interp", default=interp,
                      help="Surface reflectance channel interpolation. Choose from 'interpolate' or 'exact' "\
                      "(default=%s)"\
                      %interp )      

    parser.add_option("-b", "--band", dest="i_band",
                      help="Surface reflectance band index. Required if interp=='exact'" )  

    parser.add_option("-a", "--additional",
                      action="store_true", dest="additional_output",default=False,
                      help="Turn on writing additional_output. (default=False)")                           

    parser.add_option("-v", "--verbose",action="store_true",
                      dest="verbose", default=False,
                      help="Verbose (default=False)" )    


    parser.add_option("-n", "--nodemax", dest="nodemax", default=nodemax,
                      help="Max number of nodes to use. "\
                      "(default=%s)"\
                      %nodemax )    

    parser.add_option("-l", "--layout", dest="layout", default=layout,
                      help="Layout of domain. Used for breaking up high-res domains (e.g. GOES-R)"\
                      "(default=%s)"\
                      %layout )  

    parser.add_option("-C", "--cldmax", dest="CLDMAX", default=CLDMAX,
                      help="Maximum cloud fraction allowed. "\
                      "(default=%s)"\
                      %CLDMAX )   

    parser.add_option("-d", "--dir", dest="nccs", default=nccs,
                      help="Root level directory for inputs/outputs "\
                      "(default=%s)"\
                      %nccs )                                 

    parser.add_option("-p", "--profile", 
                      action="store_true",dest="profile", default=False,
                      help="Leave slurm output files in working directory."\
                      "(default=%s)"\
                      %profile )   

    parser.add_option("-r", "--runmode", dest="runmode", default=runmode,
                      help="VLIDORT run mode. Either 'scalar' or 'vector' "\
                      "(default=%s)"\
                      %runmode )          


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
    workspace = WORKSPACE(startdate,enddate,options)

    # Submit and Handle Jobs
    # -----------------------------
    if (workspace.dirstring) > 0:
        workspace.handle_jobs()
    else:
        print 'No model hours to run'

    


