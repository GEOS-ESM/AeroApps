#!/usr/bin/env python3

"""
    Wrapper to submit jobs to sbatch
"""

import os
import subprocess
import shutil
from   datetime        import datetime, timedelta
from   dateutil.parser import parse         as isoparser
import argparse
import numpy as np
import time

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
    """ Create slurm scripts for running run_lidar_sampler.py """
    def __init__(self,args):

        self.Date      = isoparser(args.iso_t1)
        self.enddate   = isoparser(args.iso_t2)
        self.Dt        = args.DT_hours
        self.dt        = timedelta(days=args.dt_days)

        if not os.path.exists(args.tmp):
            os.makedirs(args.tmp)

        self.cwd         = os.getcwd()
        self.slurm       = args.slurm
        self.prep_config = args.prep_config
        self.tmp         = args.tmp
        self.profile     = args.profile
        self.nproc       = args.nproc
        self.exp         = args.exp
        self.tle_year    = args.tle_year

        # create working directories
        self.create_workdir()
        self.errTally    = np.ones(len(self.dirstring)).astype(bool)

        # modify slurm scripts
        self.edit_slurm()


    def create_workdir(self):
        sdate = self.Date
        self.dirstring = []
        # create directory
        while sdate < self.enddate:

            workpath = '{}/{}'.format(self.tmp,sdate.isoformat())
            if os.path.exists(workpath):
                shutil.rmtree(workpath)

            os.makedirs(workpath)

            # copy over slurm scipt
            outfile = '{}/{}'.format(workpath,self.slurm)
            shutil.copyfile(self.slurm,outfile)

            # copy over pcf file
            outfile = '{}/{}'.format(workpath,self.prep_config)
            shutil.copyfile(self.prep_config,outfile)

            #link over needed python scripts
            source = ['lidar_sampler.py','run_lidar_sampler.py'] #,'sampling','tle']
            for src in source:
                os.symlink('{}/{}'.format(self.cwd,src),'{}/{}'.format(workpath,src))

            self.dirstring.append(workpath)

            sdate += self.dt

    def edit_slurm(self):
        for workpath in self.dirstring:
            isodate = os.path.basename(workpath)
            sdate = isoparser(isodate)
            edate = sdate + self.dt
            if edate > self.enddate:
                edate = self.enddate
        

            outpath = '{}/{}'.format(workpath,self.slurm)

            # read file first
            f = open(outpath)

            text = []
            for l in f:
                text.append(l)

            # replace one line
            iso1 = sdate.isoformat()
            iso2 = edate.isoformat()
            if self.tle_year is not None:
                newline = 'nohup python -u run_lidar_sampler.py -v --exp {} --tle_year {} --nproc {} --DT_hours {} {} {} {} >'.format(self.exp,self.tle_year,self.nproc,self.Dt,iso1,iso2,self.prep_config) + ' slurm_${SLURM_JOBID}_py.out\n'
            else:
                newline = 'nohup python -u run_lidar_sampler.py -v --exp {} --nproc {} --DT_hours {} {} {} {} >'.format(self.exp,self.nproc,self.Dt,iso1,iso2,self.prep_config) + ' slurm_${SLURM_JOBID}_py.out\n'
            text[-4] = newline
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
            os.remove(self.prep_config)

        # remove symlinks
        source = ['lidar_sampler.py','run_lidar_sampler.py'] #'sampling','tle']
        for src in source:
            os.remove(src)

        os.chdir(self.cwd)
        if self.profile is False:
            os.rmdir(self.dirstring[i])


if __name__ == '__main__':
    
    #Defaults
    DT_hours = 1
    dt_days  = 75
    nproc    = 8
    slurm    = 'run_lidar_sampler.j'
    tmp      = '/discover/nobackup/projects/gmao/osse2/pub/c1440_NR/OBS/A-CCP/workdir/lidar_sampler'
    exp      = 'g5nr'
    tle_year = None

    parser = argparse.ArgumentParser()
    parser.add_argument("iso_t1",help='starting iso time')
    parser.add_argument("iso_t2",help='ending iso time')
    parser.add_argument("prep_config",
                        help="prep config filename")

    parser.add_argument('-e',"--exp", default=exp,
                        help="GEOS experiment name for outfiles (default=%s)"%exp)

    parser.add_argument("--tle_year", default=tle_year,
                        help="Year for TLE calc (default=None)")

    parser.add_argument('-D',"--DT_hours", default=DT_hours, type=int,
                        help="Timestep in hours for each file (default=%i)"%DT_hours)

    parser.add_argument('-d',"--dt_days", default=dt_days, type=int,
                        help="Timestep in days for each job (default=%i)"%dt_days)

    parser.add_argument('-s',"--slurm",default=slurm,
                        help="slurm script template (default=%s)"%slurm)           

    parser.add_argument('-t','--tmp',default=tmp,
                        help="temp directory (default=%s)"%tmp)

    parser.add_argument("-r", "--dryrun",action="store_true",
                        help="do a dry run (default=False).") 

    parser.add_argument("-p", "--profile",action="store_true",
                        help="Don't cleanup slurm files (default=False).")   

    parser.add_argument("-n", "--nproc",default=nproc,
                        help="Number of processors (default=%i)"%nproc)                            


    args = parser.parse_args()

    workspace = WORKSPACE(args)

    if not args.dryrun:
        workspace.handle_jobs()

