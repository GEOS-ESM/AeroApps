#!/usr/bin/env python
# -W ignore::DeprecationWarning

"""
Create and monitor jobs that call aeronet_lc.py
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
import argparse

jobsmax   = 150
archive = '/archive/u/rgovinda/osse2/'


class JOBS(object):
    def handle_jobs(self):
        # Figure out how many jobs you need to submit
        runlen  = len(self.dirstring)   
        numjobs = runlen
        
        devnull = open(os.devnull, 'w')
        if numjobs <= jobsmax:   
            countRun   = runlen    
            node_tally   = numjobs  
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
                node_tally  = node_tally - len(finishedJobs)

                workingJobs = np.delete(workingJobs,finishedJobs)

            # Add more jobs if needed
            # reinitialize stat variable
            if (runlen > countRun) and (node_tally < jobsmax):
                #print 'adding new jobs'
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
        errfile = 'slurm_' +jobid + '.err'
        statinfo = os.stat(errfile)
        if (statinfo.st_size != 0):
            error = True

        os.chdir(self.cwd)
        return error


class WORKSPACE(JOBS):
    """ Create Working Directories and slurm scripts """
    def __init__(self,options):
        self.startdate = parse(options.iso_t1)
        self.enddate   = parse(options.iso_t2)
        self.dt        = timedelta(hours=options.DT_hours)

        for oo in options.__dict__:
            if (type(options.__dict__[oo]) is str) and (options.__dict__[oo].lower() == 'none'):
                self.__dict__[oo] = None
            else:
                self.__dict__[oo] = options.__dict__[oo]

        self.nccs    = self.nccs + '/' + self.instname.upper() + '/' 
        self.prefix  = self.nccs + 'workdir/'
        self.execFile = 'aeronet_lc.py'

        ##################################################
        ####
        #    Loop through dates 
        #    Create working directories, and SLURM scripts
        ####
        ##################################################

        # save run directory
        if self.verbose:
            print '++Saving run directory',os.getcwd()
        self.cwd     = os.getcwd()

        #initialize arrays to hold directory names
        self.dirstring    = []

        # Loop over dates
        startdate = self.startdate
        while startdate <= self.enddate:
            # Create working directories for intermediate inputs
            # save directory names
            workdir = self.make_workspace(startdate)
            self.dirstring.append(workdir)

            startdate = startdate + self.dt


    def make_workspace(self,date):
        # creating working directory
        dirname = '{}/{}.{}'.format(self.prefix,self.instname.lower(),date.strftime('%Y%m%dT%H'))
        jobname = '{}.{}'.format(self.instname.lower(),date.strftime('%Y%m%dT%H'))

        bindir = dirname
        if os.path.exists(dirname):
            shutil.rmtree(dirname)
        os.makedirs(dirname)    

        # create some symlinks to shared files
        os.chdir(dirname)
        if not os.path.isfile(os.path.basename(self.execFile)):
            os.symlink(self.cwd+'/'+self.execFile,os.path.basename(self.execFile))
        if not os.path.exists('ExtData'):
            os.symlink(self.cwd+'/ExtData','ExtData')
        if not os.path.isfile('Chem_MieRegistry.rc'):
            os.symlink(self.cwd+'/Chem_MieRegistry.rc','Chem_MieRegistry.rc')
        if not os.path.isfile('Aod_EOS.rc'):
            os.symlink(self.cwd+'/Aod_EOS.rc','Aod_EOS.rc')
        if not os.path.exists('nccs'):
            os.symlink(self.cwd+'/nccs','nccs')
        
        # Create options string
        options = "-D {}".format(self.DT_hours)
        if self.verbose:
            options += ' --verbose'
        if self.dryrun:
            options += ' --dryrun'
        if self.channels is not None:
            options += ' --channels {}'.format(' '.join(self.channels))
        if self.extch is not None:
            options += ' --extch {}'.format(' '.join(self.extch))
        if self.norad:
            options += ' --norad'
        if self.noext:
            options += ' --noext'

        # Create slurm runfile - using self.runfile as template           
        source = open(self.cwd+'/'+self.runfile,'r')
        destination = open(self.runfile,'w')
        for line in source:
            if (line[0:18] == '#SBATCH --job-name'):            
                destination.write('#SBATCH --job-name='+jobname+'\n')

            elif (line[0:12] == 'setenv START'):
                destination.write('setenv START '+date.isoformat()+'\n')

            elif (line[0:10] == 'setenv END'):
                destination.write('setenv END '+date.isoformat()+'\n')

            elif (line[0:14] == 'setenv OPTIONS'):
                destination.write('setenv OPTIONS "'+options+'"\n')

            elif (line[0:14] == 'setenv AEROBIN'):
                destination.write('setenv AEROBIN '+bindir+'\n')

            else:
                destination.write(line)        
        source.close()
        destination.close()
        
        
        # Go back to original run directory 
        os.chdir(self.cwd) 
        
        return dirname



    def destroy_workspace(self,i,jobid):
        # put LevelB files in archive or remove
        # --------------------

        #parse date from distring
        date = parse(os.path.basename(self.dirstring[i]).split('.')[1])

        if self.layout is not None:
            # parse layout code
            layout = os.path.basename(self.dirstring[i]).split('.')[-1]
        else:
            layout = None

        g5dir = self.indir + '/LevelB/'+ 'Y'+ str(date.year) + '/M' + str(date.month).zfill(2) + '/D' + str(date.day).zfill(2) 
        nymd  = str(date.year) + str(date.month).zfill(2) + str(date.day).zfill(2)
        hour  = str(date.hour).zfill(2)

        if layout is None:
            met   = g5dir + '/' + self.instname.lower() + '-g5nr.lb2.met_Nv.' + nymd + '_' + hour + 'z.nc4'
            chm   = g5dir + '/' + self.instname.lower() + '-g5nr.lb2.chm_Nv.' + nymd + '_' + hour + 'z.nc4'
            aer   = g5dir + '/' + self.instname.lower() + '-g5nr.lb2.aer_Nv.' + nymd + '_' + hour + 'z.nc4'                            
            geom  = g5dir + '/' + self.instname.lower() + '.lb2.angles.' + nymd + '_' + hour + 'z.nc4'
            land  = self.indir + '/LevelB/invariant/' + self.instname.lower() + '-g5nr.lb2.asm_Nx.nc4'  
        else:
            met   = g5dir + '/' + self.instname.lower() + '-g5nr.lb2.met_Nv.' + nymd + '_' + hour + 'z_' + layout +'.nc4'
            chm   = g5dir + '/' + self.instname.lower() + '-g5nr.lb2.chm_Nv.' + nymd + '_' + hour + 'z_' + layout +'.nc4'
            aer   = g5dir + '/' + self.instname.lower() + '-g5nr.lb2.aer_Nv.' + nymd + '_' + hour + 'z_' + layout +'.nc4'
            geom  = g5dir + '/' + self.instname.lower() + '.lb2.angles.' + nymd + '_' + hour + 'z_' + layout +'.nc4'
            land  = self.indir + '/LevelB/invariant/' + self.instname.lower() + '-g5nr.lb2.asm_Nx_' + layout + '.nc4'  

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
        
#########################################################

if __name__ == "__main__":
    # Defaults
    # ------------
    DT_hours          = 1
    instname          = 'aeronet'    
    nccs              = '/discover/nobackup/projects/gmao/osse2/pub/c1440_NR/OBS/'
    
    runfile           = '{}_lc.j'.format(instname)
    
    #Flags
    # verbose           = False
    # keep_lb           = False
    # archive_lb        = False

    # Parse command line options
    # ------------------------------
#   Parse command line options
#   --------------------------

    parser = argparse.ArgumentParser()
    parser.add_argument("iso_t1",
                        help="starting iso time")
    parser.add_argument("iso_t2",
                        help="ending iso time")

    parser.add_argument("-i","--instname", default=instname,
                        help="Instrument name (default=%s)"%instname)   

    parser.add_argument("-j","--runfile", default=runfile,
                        help="slurm script template (default=%s)"%runfile)   

    parser.add_argument("-D","--DT_hours", default=DT_hours, type=int,
                        help="Timestep in hours for each file (default=%i)"%DT_hours)

    parser.add_argument("-v", "--verbose",action="store_true",
                        help="Verbose mode (default=False).")

    parser.add_argument("-r", "--dryrun",action="store_true",
                        help="do a dry run (default=False).")   

    parser.add_argument("-c","--channels", default=None,nargs='+',
                        help="channels to get BOA radiance")                            

    parser.add_argument("-e","--extch", default=None,nargs='+',
                        help="channels to run extinction sampler")  

    parser.add_argument("--norad",action="store_true",
                        help="No radiance calculations (default=False).")

    parser.add_argument("--noext",action="store_true",
                        help="No extinctions calculations (default=False).")

    parser.add_argument("-n", "--nccs", default=nccs,
                      help="Root level directory for inputs/outputs "\
                      "(default=%s)"\
                      %nccs )                                 


    ################
    ###
    #    End of uper inputs
    ###
    ################    
    args = parser.parse_args()


    # Setup workspace and runscript
    # ------------------------------
    workspace = WORKSPACE(args)

    # # Submit and Handle Jobs
    # # -----------------------------
    # if (workspace.dirstring) > 0:
    #     workspace.handle_jobs()
    # else:
    #     print 'No model hours to run'

 