#!/usr/bin/env python3

"""
    Wrapper to submit jobs to sbatch
"""

import os
import subprocess
import shutil
from   datetime        import datetime, timedelta
from   dateutil.parser import parse         as isoparser
from   dateutil.relativedelta import relativedelta
import argparse
import numpy as np
import time
from   MAPL.config      import Config

class JOBS(object):
    def handle_jobs(self):
        jobsmax = 25
        # Figure out how many jobs you need to submit
        numjobs  = len(self.dirstring)   
        
        devnull = open(os.devnull, 'w')
        if numjobs <= jobsmax:   
            countRun     = numjobs 
            node_tally   = numjobs
        else:
            countRun   = jobsmax
            node_tally = jobsmax

        workingJobs = np.arange(countRun)
        

        # Submit first JOBSMAX jobs
        jobid = np.empty(0)
        for i in workingJobs:
            s = self.dirstring[i]
            os.chdir(s)
            result = subprocess.check_output(['sbatch',self.slurm])
            jobid = np.append(jobid,result.split()[-1])
        os.chdir(self.cwd)

        # launch subprocess that will monitor queue

        # Monitor jobs 1-by-1 
        # Add a new job when one finishes 
        # Until they are all done
        stat = 0
        countDone = 0
        while (stat == 0):
            finishedJobs = np.empty(0,dtype='int')
            for ii,i in enumerate(workingJobs):
                s = jobid[i]
                s = s.strip('\n')
                finished = False

                # Check to see if this job is finished
                result = subprocess.call(['qstat',s], stdout=devnull)
                if (result != 0):
                    finished = True   

                # if the job is finished add to finished jobs list
                # and clean up workspace
                if finished:
                    #print 'Job finished, cleaning up', s, i 
                    finishedJobs = np.append(finishedJobs,ii)
                    countDone += 1
                    errcheck = self.check_for_errors(i,s)               
                    if (errcheck is False):
                        self.errTally[i] = False

                        # Clean up workspaces
                        self.destroy_workspace(i,s)        
                    else:
                        print('Jobid ',s,' in ',self.dirstring[i],' exited with errors')

            # finished checking up on all the jobs
            # Remove finished jobs from the currently working list
            if len(finishedJobs) != 0:
                print('deleting finishedJobs',finishedJobs,jobid[workingJobs[finishedJobs]])
                node_tally  = node_tally - len(finishedJobs)

                workingJobs = np.delete(workingJobs,finishedJobs)

            # Add more jobs if needed
            # reinitialize stat variable
            if (numjobs > countRun) and (node_tally < jobsmax):
                #print 'adding new jobs'
                newRun     = jobsmax - node_tally
                node_tally = jobsmax

                if (newRun + countRun) > numjobs:
                    newRun = numjobs - countRun
                    node_tally = node_tally + newRun
                

                newjobs  = countRun + np.arange(newRun)
                workingJobs = np.append(workingJobs, newjobs)
                for i in newjobs:
                    s = self.dirstring[i]
                    os.chdir(s)
                    result = subprocess.check_output(['sbatch',self.slurm])
                    jobid = np.append(jobid,result.split()[-1])

                os.chdir(self.cwd)
                countRun = countRun + newRun

            # check if all the jobs are finished
            if countDone == numjobs:
                stat = 1
            else:
                print('Waiting 30 minutes')
                time.sleep(60*30)


        # Exited while loop
        print('All jobs done')

        # Postprocessing done
        print('Cleaned Up Worksapces')
        devnull.close()


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
    """ Create slurm scripts for running accp_polar_vlidort.py """
    def __init__(self,args):

        self.Date      = isoparser(args.iso_t1)
        self.enddate   = isoparser(args.iso_t2)
        self.Dt        = relativedelta(months=args.DT_months)

        
        self.track_pcf   = args.track_pcf
        self.orbit_pcf   = args.orbit_pcf
        self.inst_pcf    = args.inst_pcf
        self.cirrus_model = args.cirrus_model
        self.DT_months    = args.DT_months
        self.albedoType  = args.albedotype
        self.rcFile      = args.rcfile
        self.dryrun      = args.dryrun
        self.verbose     = args.verbose

        self.slurm       = args.slurm
        self.tmp         = args.tmp
        self.profile     = args.profile

        if not os.path.exists(args.tmp):
            os.makedirs(args.tmp)

        self.cwd         = os.getcwd()

        # figure out channels from instfile
        cf = Config(args.inst_pcf,delim=' = ')
        channels = cf('channels')
        if ',' in channels:
            channels = channels.split(',')    
        else:
            channels = [channels]    
        self.channels = np.array(channels).astype(int)
        self.nch = len(channels)
        
        # create working directories
        self.create_workdir()
        self.errTally    = np.ones(len(self.dirstring)).astype(bool)


    def create_workdir(self):
        sdate = self.Date
        self.dirstring = []
        while sdate < self.enddate:
            for ch in self.channels:
                # create directory

                workpath = '{}/{}.{}'.format(self.tmp,sdate.isoformat(),ch)
                if os.path.exists(workpath):
                    shutil.rmtree(workpath)

                os.makedirs(workpath)

                # copy over slurm scipt
                outfile = '{}/{}'.format(workpath,self.slurm)
                shutil.copyfile(self.slurm,outfile)

                # copy over pcf files                
                outfile = '{}/{}'.format(workpath,self.track_pcf)
                shutil.copyfile(self.track_pcf,outfile)

                outfile = '{}/{}'.format(workpath,self.orbit_pcf)
                shutil.copyfile(self.orbit_pcf,outfile)

                outfile = '{}/{}'.format(workpath,self.inst_pcf)
                shutil.copyfile(self.inst_pcf,outfile)

                #link over needed python scripts
                source = ['accp_polar_vlidort_cirrus.py'] 
                for src in source:
                    os.symlink('{}/{}'.format(self.cwd,src),'{}/{}'.format(workpath,src))

                # Copy over rc and edit
                outfile = '{}/{}'.format(workpath,self.rcFile)                

                source = open(self.rcFile,'r')
                destination = open(outfile,'w')
                a = float(ch)*1e-3
                for line in source:
                    if (line[0:11] == 'r_channels:'):
                        destination.write('r_channels: '+'{:0.3f}e-6'.format(a)+'\n')
                    else:
                        destination.write(line)
                source.close()
                destination.close()

                # edit slurm
                self.edit_slurm(sdate,workpath,ch)

                self.dirstring.append(workpath)
            sdate += self.Dt

    def edit_slurm(self,sdate,workpath,channel):
        edate = sdate + self.Dt
        outpath = '{}/{}'.format(workpath,self.slurm)

        # read file first
        f = open(outpath)

        text = []
        for l in f:
            text.append(l)

        # replace one line
        iso1 = sdate.isoformat()
        iso2 = edate.isoformat()
        Options = ' --DT_months {}'.format(self.DT_months) +\
                  ' --rcfile {}'.format(self.rcFile) 

        if self.albedoType is not None:
            Options += ' --albedotype {}'.format(self.albedoType)
        if self.dryrun:
            Options += ' -r'
        if self.verbose:
            Options +=  ' -v' 
        newline = 'python -u ./accp_polar_vlidort_cirrus.py {} {} {} {} {} {} {} {} >'.format(Options,iso1,iso2,self.track_pcf,self.orbit_pcf,self.inst_pcf,channel,self.cirrus_model) + ' slurm_${SLURM_JOBID}_py.out\n'
        text[-5] = newline
        f.close()

        #  write out
        f = open(outpath,'w')
        for l in text:
            f.write(l)
        f.close()             

    def destroy_workspace(self,i,jobid):
        os.chdir(self.dirstring[i])

        if self.profile is False:

            errfile = 'slurm_' +jobid + '.err'
            os.remove(errfile)        
            outfile = 'slurm_' +jobid + '.out'
            os.remove(outfile)        

            outfile = 'slurm_' +jobid + '_py.out'
            os.remove(outfile)     

            os.remove(self.slurm)
            os.remove(self.track_pcf)
            os.remove(self.orbit_pcf)
            os.remove(self.inst_pcf)
            os.remove(self.rcFile)

        # remove symlinks
        source = ['accp_polar_vlidort_cirrus.py'] 
        for src in source:
            os.remove(src)

        os.chdir(self.cwd)
        if self.profile is False:
            os.rmdir(self.dirstring[i])


if __name__ == '__main__':
    
    #Defaults
    DT_months = 1
    slurm    = 'run_accp_polar_vlidort_cirrus.j'
    tmp      = '/discover/nobackup/projects/gmao/osse2/pub/c1440_NR/OBS/A-CCP/workdir/vlidort'
    rcFile   = 'Aod_EOS.rc'
    albedoType = None

    parser = argparse.ArgumentParser()
    parser.add_argument("iso_t1",help='starting iso time')
    parser.add_argument("iso_t2",help='ending iso time')
    parser.add_argument("track_pcf",
                        help="prep config file with orbit variables")

    parser.add_argument("orbit_pcf",
                        help="prep config file with orbit variables")

    parser.add_argument("inst_pcf",
                        help="prep config file with instrument variables")

    parser.add_argument("cirrus_model",
                        help="aggregates, rosette0, rosette3, rosette50")

    parser.add_argument('-D',"--DT_months", default=DT_months, type=int,
                        help="Timestep in months for each file (default=%i)"%DT_months)

    parser.add_argument("-a","--albedotype", default=albedoType,
                        help="albedo type keyword. default is to figure out according to channel")

    parser.add_argument("--rcfile",default=rcFile,
                        help="rcFile (default=%s)"%rcFile)    

    parser.add_argument('-s',"--slurm",default=slurm,
                        help="slurm script template (default=%s)"%slurm)           

    parser.add_argument('-t','--tmp',default=tmp,
                        help="temp directory (default=%s)"%tmp)

    parser.add_argument("-r", "--dryrun",action="store_true",
                        help="do a dry run (default=False).") 

    parser.add_argument("-p", "--profile",action="store_true",
                        help="Don't cleanup slurm files (default=False).")   

    parser.add_argument("-v", "--verbose",action="store_true",
                        help="Verbose mode (default=False).")

    args = parser.parse_args()

    workspace = WORKSPACE(args)

    if not args.dryrun:
        workspace.handle_jobs()

