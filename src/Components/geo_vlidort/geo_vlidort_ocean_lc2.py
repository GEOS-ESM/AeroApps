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
from   optparse        import OptionParser   # Command-line args

jobsmax   = 150
dt = timedelta(hours=1)
archive = '/archive/u/rgovinda/osse2/'

mr =  [1.396,1.362,1.349,1.345,1.339,1.335,1.334,1.333,1.332,1.331,1.329,1.326,
      1.323,1.318,1.312,1.306,1.292,1.261]

mr_ch = [200,250,300,337,400,488,515,550,633,694,860,1060,1300,1536,1800,2000,2250,2500]

mr = np.array(mr)
mr_ch = np.array(mr_ch)

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
                    for a in np.arange(self.nodemax_list[i]):
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
                    s = self.dirstring[i]
                    os.chdir(s)
                    jobid = np.append(jobid,subprocess.check_output(['qsub',runfile]))

                os.chdir(self.cwd)
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

        self.cloud = False
        for oo in options.__dict__:
            if (type(options.__dict__[oo]) is str) and (options.__dict__[oo].lower() == 'none'):
                self.__dict__[oo] = None
            else:
                self.__dict__[oo] = options.__dict__[oo]

        if self.nodemax is not None: self.nodemax = int(self.nodemax)
        self.nccs    = self.nccs + '/' + self.instname.upper() + '/DATA/' 
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

        if type(self.runmode) is str:
            if ',' in self.runmode:
                self.runmodeL = self.runmode.replace(' ','').split(',')
            else:
                self.runmodeL = self.runmode.split()

        self.startdateL = startdate
        self.enddateL   = enddate

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
        
        for startdate,enddate in zip(self.startdateL,self.enddateL):
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

                        for runmode in self.runmodeL:
                            # Runmode code
                            self.runmode = runmode
                            self.code = runmode + '.'
                            self.code += self.surface

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
                                if (self.surface.upper() == 'CX'):
                                    self.make_cx_rcfile(workdir,startdate,ch,nodemax=nodemax,
                                                   layout=laycode)
                            

            
                startdate = startdate + dt

        self.nodemax_list = np.array(self.nodemax_list)

    def get_from_archive(self,path):
        ii = path.find('/c1440_NR/')
        archive = self.archive + path[ii:]

        if not os.path.exists(os.path.dirname(path)):
            os.makedirs(os.path.dirname(path))

        # Copy from archive.  Note: overwrites file if it exists.
        try:
            shutil.copyfile(archive,path) 
        except IOError:
            print 'Could not find archive ',archive, ' to put in ', path
            sys.exit()

    def put_in_archive(self,path):
        ii = path.find('/c1440_NR/')
        archive = self.archive + path[ii:]
        # Copy to archive.  Note: overwrites file if it exists.
        if not os.path.exists(os.path.dirname(archive)):
            os.makedirs(os.path.dirname(archive))
        try:
            shutil.copyfile(path,archive) 
        except IOError:
            print 'Could not put '+path+' in archive ',archive
            sys.exit()


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
            aer   = g5dir + '/' + self.instname.lower() + '-g5nr.lb2.aer_Nv.' + nymd + '_' + hour + 'z.nc4'            

        else:
            met   = g5dir + '/' + self.instname.lower() + '-g5nr.lb2.met_Nv.' + nymd + '_' + hour + 'z_' + layout +'.nc4'
            geom  = g5dir + '/' + self.instname.lower() + '.lb2.angles.' + nymd + '_' + hour + 'z_' + layout +'.nc4'
            land  = self.indir + '/LevelB/invariant/' + self.instname.lower() + '-g5nr.lb2.asm_Nx_' + layout + '.nc4'  
            aer   = g5dir + '/' + self.instname.lower() + '-g5nr.lb2.aer_Nv.' + nymd + '_' + hour + 'z_' + layout +'.nc4'


        if self.verbose:
            print '++Opening metfile ',met

        if not os.path.exists(met):
            self.get_from_archive(met)
        ncMet = Dataset(met)
        Cld   = np.squeeze(ncMet.variables[u'CLDTOT'][:])
        f     = np.where(Cld <= float(self.CLDMAX))
        ncMet.close()
        if len(f[0]) == 0:
            return False, 0

        if not os.path.exists(geom):
            self.get_from_archive(geom)
        ncGeom = Dataset(geom)
        SZA    = np.squeeze(ncGeom.variables[u'solar_zenith'][:])
        VZA    = np.squeeze(ncGeom.variables[u'sensor_zenith'][:])
        ncGeom.close()

        if not os.path.exists(land):
            self.get_from_archive(land)        
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
        f      = np.where(FRLAND < 0.99)
        if len(f[0]) == 0:
            return 0

        return len(f[0])


    def make_workspace(self,date,ch,nodemax=None,layout=None):

        # Get necessary files from archive if needed
        g5dir = self.indir + '/LevelB/'+ 'Y'+ str(date.year) + '/M' + str(date.month).zfill(2) + '/D' + str(date.day).zfill(2) 
        nymd  = str(date.year) + str(date.month).zfill(2) + str(date.day).zfill(2)
        hour  = str(date.hour).zfill(2)

        if (layout is None) or (self.cloud is True) :
            chm   = g5dir + '/' + self.instname.lower() + '-g5nr.lb2.chm_Nv.' + nymd + '_' + hour + 'z.nc4'
            aer   = g5dir + '/' + self.instname.lower() + '-g5nr.lb2.aer_Nv.' + nymd + '_' + hour + 'z.nc4'            
        else:
            chm   = g5dir + '/' + self.instname.lower() + '-g5nr.lb2.chm_Nv.' + nymd + '_' + hour + 'z_' + layout +'.nc4'
            aer   = g5dir + '/' + self.instname.lower() + '-g5nr.lb2.aer_Nv.' + nymd + '_' + hour + 'z_' + layout +'.nc4'

        if (not os.path.exists(chm)) and (self.cloud is False):
            self.get_from_archive(chm)
        if not os.path.exists(aer):
            self.get_from_archive(aer)        

        # creating working directory
        dirname = '{}/{}.{}T{}.{}.{}'.format(self.prefix,self.instname.lower(),date.date(),str(date.hour).zfill(2),ch,self.code)
        jobname = '{}.{}T{}.{}'.format(self.instname.lower(),date.date(),str(date.hour).zfill(2),ch)

        if layout is not None:
            dirname = dirname + '.' + layout
            jobname = jobname + '.' + layout

        bindir = dirname
        if os.path.exists(dirname):
            shutil.rmtree(dirname)
        os.makedirs(dirname)    


        # outdirname
        outdir = '{}/Y{}/M{}/D{}'.format(self.outdir,date.year,str(date.month).zfill(2),str(date.day).zfill(2))
        if self.addoutdir is not None:
            addoutdir = '{}/Y{}/M{}/D{}'.format(self.addoutdir,date.year,str(date.month).zfill(2),str(date.day).zfill(2))

        # create some symlinks to shared files
        os.chdir(dirname)
        if not os.path.isfile(os.path.basename(self.execFile)):
            os.symlink(self.execFile,os.path.basename(self.execFile))
        if not os.path.exists('ExtData'):
            os.symlink(self.cwd+'/ExtData','ExtData')
        if self.cloud is True:
            if not os.path.exists('ExtDataCloud'):
                os.symlink(self.cwd+'/ExtDataCloud','ExtDataCloud')            
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
                    destination.write('$RUN_CMD  ./'+ os.path.basename(self.execFile) +' geo_vlidort.rc ${SLURM_ARRAY_TASK_ID}'+'\n')    
                else:
                    destination.write('$RUN_CMD  ./'+ os.path.basename(self.execFile) + ' geo_vlidort.rc'+'\n')   

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

    def make_cx_rcfile(self,dirname,date,ch,nodemax=None,layout=None):
        os.chdir(dirname)
        fch = float(ch)

        rcfile = open('geo_vlidort.rc','w')
        rcfile.write('INDIR: '+self.indir+'\n')
        rcfile.write('OUTDIR: .\n')
        rcfile.write('DATE: '+str(date.year)+str(date.month).zfill(2)+str(date.day).zfill(2)+'\n')
        rcfile.write('TIME: '+str(date.hour).zfill(2)+'\n')
        rcfile.write('INSTNAME: ' + self.instname.lower() + '\n')
        rcfile.write('SURFNAME: GISS_CoxMunk\n')
        rcfile.write('SURFMODEL: CX\n')

        mruse = np.interp(float(ch),mr_ch,mr)
        rcfile.write('SURFMR: {}\n'.format(mruse))

        if (self.runmode == 'scalar'):
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

        if layout is not None:
            rcfile.write('LAYOUT: '+layout+'\n')
            
        rcfile.close()

        os.chdir(self.cwd)
           

    def final_cleanup(self):
        # put LevelB aer files in archive or remove
        # --------------------
        # Loop over dates
        date = self.startdate
        while (date <= self.enddate):
            g5dir = self.indir + '/LevelB/'+ 'Y'+ str(date.year) + '/M' + str(date.month).zfill(2) + '/D' + str(date.day).zfill(2) 
            nymd  = str(date.year) + str(date.month).zfill(2) + str(date.day).zfill(2)
            hour  = str(date.hour).zfill(2)

            if self.layout is not None:
                # loop through tiles
                for tile in np.arange(self.ntiles):
                    laycode = self.layout + str(tile)
                    met   = g5dir + '/' + self.instname.lower() + '-g5nr.lb2.met_Nv.' + nymd + '_' + hour + 'z_' + laycode +'.nc4'
                    chm   = g5dir + '/' + self.instname.lower() + '-g5nr.lb2.chm_Nv.' + nymd + '_' + hour + 'z_' + laycode +'.nc4'
                    aer   = g5dir + '/' + self.instname.lower() + '-g5nr.lb2.aer_Nv.' + nymd + '_' + hour + 'z_' + laycode +'.nc4'
                    geom  = g5dir + '/' + self.instname.lower() + '.lb2.angles.' + nymd + '_' + hour + 'z_' + laycode +'.nc4'
                    land  = self.indir + '/LevelB/invariant/' + self.instname.lower() + '-g5nr.lb2.asm_Nx_' + laycode + '.nc4'  

                    if self.archive_lb:
                        self.put_in_archive(met)
                        self.put_in_archive(chm)
                        self.put_in_archive(aer)
                        self.put_in_archive(geom)
                        self.put_in_archive(land)
                        
                    if not self.keep_lb:
                        os.remove(met)
                        os.remove(chm)
                        os.remove(aer)
                        os.remove(geom)
                        os.remove(land)

            else:
                met   = g5dir + '/' + self.instname.lower() + '-g5nr.lb2.met_Nv.' + nymd + '_' + hour + 'z.nc4'
                chm   = g5dir + '/' + self.instname.lower() + '-g5nr.lb2.chm_Nv.' + nymd + '_' + hour + 'z.nc4'
                aer   = g5dir + '/' + self.instname.lower() + '-g5nr.lb2.aer_Nv.' + nymd + '_' + hour + 'z.nc4'                            
                geom  = g5dir + '/' + self.instname.lower() + '.lb2.angles.' + nymd + '_' + hour + 'z.nc4'
                land  = self.indir + '/LevelB/invariant/' + self.instname.lower() + '-g5nr.lb2.asm_Nx.nc4'  

                if self.archive_lb:
                    self.put_in_archive(met)
                    self.put_in_archive(chm)
                    self.put_in_archive(aer)
                    self.put_in_archive(geom)
                    self.put_in_archive(land)
                    
                if not self.keep_lb:
                    os.remove(met)
                    os.remove(chm)
                    os.remove(aer)
                    os.remove(geom)
                    os.remove(land)

            date += dt


    def destroy_workspace(self,i,jobid):
        os.chdir(self.dirstring[i])

        os.remove('Aod_EOS.rc')
        os.remove('Chem_MieRegistry.rc')
        os.remove(os.path.basename(self.execFile))
        os.remove('clean_mem.sh')
        os.remove('ExtData')
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
    surface           = 'CX' #'CX'
        
    nodemax           = 10
    layout            = None
    CLDMAX            = '0.01'
    nccs              = '/discover/nobackup/projects/gmao/osse2/pub/c1440_NR/OBS/'
    
    runmode           = 'vector'
    runfile           = 'geo_vlidort_lc2.j'
    execFile          = '/discover/nobackup/pcastell/workspace/GAAS/src/Components/geo_vlidort/geo_vlidort_ocean.x'
    
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

    parser.add_option("-V", "--version_string", dest="version", default=version,
                      help="Version name (default=%s)"\
                      %version )        

    parser.add_option("-c", "--channels", dest="channels", default=channels,
                      help="Channels (default=%s)"\
                      %channels )  

    parser.add_option("-s", "--surface", dest="surface", default=surface,
                      help="Surface Reflectance model.  Choose from 'CX' "\
                      "(default=%s)"\
                      %surface )      
      
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
    if ',' in startdate:
        startdate = startdate.split(',')
        enddate   = enddate.split(',')
        startL = []
        endL   = []
        for s in startdate:
            startL.append(parse(s))

        for e in enddate:
            endL.append(parse(e))
    else:
        startdate = [parse(startdate)]
        enddate   = [parse(enddate)]
    workspace = WORKSPACE(startdate,enddate,options)

    # Submit and Handle Jobs
    # -----------------------------
    if (workspace.dirstring) > 0:
        workspace.handle_jobs()
        workspace.final_cleanup()
    else:
        print 'No model hours to run'

    


