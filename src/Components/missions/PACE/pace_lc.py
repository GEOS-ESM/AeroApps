#!/usr/bin/env python

"""
    Wrapper to loop through channels and submit pace_lc.j jobs to sbatch for one pace granule
"""

import os
import subprocess
import shutil
from   datetime        import datetime, timedelta
from   dateutil.parser import parse         as isoparser
import argparse
import numpy           as np
import time
from   pace            import granules
from   netCDF4         import Dataset

jobsmax   = 150
class JOBS(object):
    def handle_jobs(self):
        # Figure out how many jobs you need to submit
        runlen  = len(self.dirstring)   
        
        devnull = open(os.devnull, 'w')
        if runlen <= jobsmax:   
            countRun     = runlen  
            node_tally   = runlen 
        else:
            countRun   = jobsmax
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
                        print 'Jobid ',s,' in ',self.dirstring[i],' exited with errors'

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
                newRun     = jobsmax - node_tally
                node_tally = jobsmax

                if (newRun + countRun) > runlen:
                    newRun = runlen - countRun
                    node_tally = node_tally + newRun
                

                newjobs  = countRun + np.arange(newRun)
                workingJobs = np.append(workingJobs, newjobs)
                for i in newjobs:
                    s = self.dirstring[i]
                    os.chdir(s)
                    jobid = np.append(jobid,subprocess.check_output(['qsub',self.runfile]))

                os.chdir(self.cwd)
                countRun = countRun + newRun
                stat = subprocess.call(['qstat -u pcastell'], shell=True, stdout=devnull)


            print 'Waiting 20 minutes'
            time.sleep(60*20)
            

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
    """ Create slurm scripts for running pace_vlidort_lc.py """
    def __init__(self,args):
        self.Date    = isoparser(args.iso_t1)

        if not os.path.exists(args.tmp):
            os.makedirs(args.tmp)

        self.cwd = os.getcwd()
        self.runfile = args.slurm
        self.profile = args.profile
        self.rootdir = args.rootdir

        self.dirstring = []

        if args.channels is None:
            self.get_channels(self.Date)
        else:
            if type(args.channels) is float:
                self.channels = [args.channels]
            else:
                self.channels = args.channels

        if args.nodemax is not None: self.nodemax = int(args.nodemax)
        if args.nodemax is None: self.nodemax = 1

        for ch in self.channels:
            outpath = '{}/{}.{}'.format(args.tmp,self.Date.isoformat(),ch)

            # Copy over some files to working temp directory
            self.create_workdir(outpath,ch)

            # # Create pcf file
            # self.write_pcf()

            self.dirstring.append(outpath)




    def create_workdir(self,outpath,ch):
        # create directory     
        if os.path.exists(outpath):
            shutil.rmtree(outpath)
        os.makedirs(outpath)

        # Copy over slurm script    
        outfile = '{}/{}'.format(outpath,self.runfile)
        shutil.copyfile(self.runfile,outfile)

        # Copy over Aod_EOS.rc
        outfile = '{}/{}'.format(outpath,'Aod_EOS.rc')
        shutil.copyfile('Aod_EOS.rc',outfile)
        self.edit_AODrc(outfile,ch)

        # link over some scripts and files
        source = ['leo_vlidort_cloud.x','ExtData','Chem_MieRegistry.rc',
                  'ExtDataCloud']
        for src in source:
            os.symlink('{}/{}'.format(self.cwd,src),'{}/{}'.format(outpath,src))

        

    def edit_AODrc(self,filename,ch):
        f = open(filename)
        # read the file first
        text = []
        for l in f:
            text.append(l)

        f.close()

        # Edit text
        for i,l in enumerate(text):
            if 'r_channels:' in l[0:10]:
                newline = "r_channels: {}e-6".format(ch*0.001)
                text[i] = newline

        # Write edited text back to file
        f = open(filename,'w')
        for l in text:
            f.write(l)

        f.close()



    def write_pcf(self):
        # Change start and enddate
        iso1 = Date.isoformat()
        iso2 = edate.isoformat()

        newline  = 'setenv START {}\n'.format(iso1)
        text[18] = newline

        newline  = 'setenv END  {}\n'.format(iso2)
        text[19] = newline



    def get_channels(self,date):
        pdate = isoparser(date.strftime('2020-%m-%dT%H:%M:00'))
        inFile = granules(self.rootdir+'/L1B',pdate,pdate)
        if len(inFile) == 0:
            raise Exception('No PACE Infiles found, nothing to do, exiting.')
        
        nc = Dataset(inFile[0])
        group = 'sensor_band_parameters'
        SDS = ['blue_wavelength','red_wavelength','SWIR_wavelength']
        channels = []
        for sds in SDS:
            channels.append(nc.groups[group].variables[sds][:])


        self.gchannels = channels
        self.channels  = list(np.concatenate(channels))


    def destroy_workspace(self,i,jobid):
        os.chdir(self.dirstring[i])

        if self.profile is False:
            errfile = 'slurm_' +jobid + '.err'
            os.remove(errfile)        
            outfile = 'slurm_' +jobid + '.out'
            os.remove(outfile)        
            pyfile = 'slurm_' +jobid + '_py.out'
            os.remove(pyfile)        

            os.remove(self.runfile)

        # remove symlinks
        source = ['nccs','Aod_EOS.rc','Aod_EOS.440.rc','Aod_EOS.870.rc',
                  'run_aeronet_vlidort.py','aeronet_vlidort_lc.py','ExtData',
                  'Chem_MieRegistry.rc']
        for src in source:
            os.remove(src)

        os.chdir(self.cwd)
        if self.profile is False:
            os.rmdir(self.dirstring[i])



if __name__ == '__main__':
    
    #Defaults
    nodemax  = 20
    slurm    = 'pace_lc_array.j'
    rootdir  = '/discover/nobackup/projects/gmao/osse2/pub/c1440_NR/OBS/PACE'
    tmp      = '/discover/nobackup/projects/gmao/osse2/pub/c1440_NR/OBS/PACE/workdir'

    parser = argparse.ArgumentParser()
    parser.add_argument("iso_t1",help='starting iso time')

    parser.add_argument("--rootdir",default=rootdir,
                        help="root directory for PACE data (default=%s)."%rootdir)       

    parser.add_argument('-s',"--slurm",default=slurm,
                        help="slurm script template (default=%s)"%slurm)           

    parser.add_argument('-t','--tmp',default=tmp,
                        help="temp directory (default=%s)"%tmp)

    parser.add_argument("-p", "--profile",action="store_true",
                        help="Don't cleanup slurm files (default=False).")    

    parser.add_argument("-v", "--verbose",action="store_true",
                        help="Verbose mode (default=False).")

    parser.add_argument("-r", "--dryrun",action="store_true",
                        help="do a dry run (default=False).")   

    parser.add_argument("-c","--channels", default=None,nargs='+',
                        help="channels to get TOA radiance (default=None - read in from PACE L1B File)")                            

    parser.add_argument("-e","--extch", default=None,nargs='+',
                        help="channels to run extinction sampler (default=None - read in from PACE L1B File)")  

    parser.add_argument("--norad",action="store_true",
                        help="No radiance calculations (default=False).")

    parser.add_argument("--noext",action="store_true",
                        help="No extinctions calculations (default=False).")

    parser.add_argument("-n", "--nodemax", default=nodemax,
                        help="Max number of nodes to use. "\
                      "(default=%s)"\
                      %nodemax )       


    args = parser.parse_args()

    workspace = WORKSPACE(args)

    # Submit and monitor jobs
    # workspace.handle_jobs()