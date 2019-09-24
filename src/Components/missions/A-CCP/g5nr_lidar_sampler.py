#!/usr/bin/env python

"""
    Wrapper to submit jobs to sbatch
"""

import os
import subprocess
import shutil
from   datetime        import datetime, timedelta
from   dateutil.parser import parse         as isoparser
import argparse

class JOBS(object):
    def handle_jobs(self):
        jobsmax = 150
        # Figure out how many jobs you need to submit
        runlen  = len(self.dirstring)   

        if self.nodemax is not None:
            numjobs = self.nodemax*runlen
        else:
            numjobs = runlen
        
        devnull = open(os.devnull, 'w')
        if numjobs <= jobsmax:   
            countRun     = runlen  
            node_tally   = runlen 
        else:
            if self.nodemax is not None:
                countRun   = int(jobsmax/self.nodemax)
                node_tally = self.nodemax*int(jobsmax/self.nodemax)
                jobsmax    = node_tally
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

        # launch subprocess that will monitor queue and do file merging

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
                if self.nodemax is not None:
                    # Loop through the nodes working on this job
                    finishedCNT = 0
                    for a in np.arange(self.nodemax):
                        a = a + 1

                        try:
                            result = subprocess.check_output(['squeue','-j',s+'_'+ str(a)],stderr=subprocess.STDOUT)
                            if (s+'_'+ str(a) not in result):
                                #print ' Job Not Found! '+s+'_'+ str(a)
                                finishedCNT = finishedCNT + 1
                        except subprocess.CalledProcessError as e:
                            #print ' ERROR Job Not Found! '+s+'_'+ str(a)
                            finishedCNT = finishedCNT + 1

                    if (finishedCNT == self.nodemax):
                        finished = True
                else:
                    result = subprocess.call(['qstat',s], stdout=devnull)
                    if (result != 0):
                        finished = True   

                # if the job is finished add to finished jobs list
                if finished:
                    #print 'Job finished, cleaning up', s, i 
                    finishedJobs = np.append(finishedJobs,ii)
                    errcheck = self.check_for_errors(i,s)               
                    if (errcheck is False):
                        self.errTally[i] = False
                    else:
                        print 'Jobid ',s,' in ',self.dirstring[i],' exited with errors'

            # finished checking up on all the jobs
            # Remove finished jobs from the currently working list
            if len(finishedJobs) != 0:
                print 'deleting finishedJobs',finishedJobs,jobid[workingJobs[finishedJobs]]
                if self.nodemax is not None:
                    node_tally = node_tally - self.nodemax*len(finishedJobs)
                else:                
                    node_tally  = node_tally - len(finishedJobs)

                workingJobs = np.delete(workingJobs,finishedJobs)

            # Add more jobs if needed
            # reinitialize stat variable
            if (runlen > countRun) and (node_tally < jobsmax):
                #print 'adding new jobs'
                if self.nodemax is not None:
                    newRun     = (jobsmax - node_tally)/self.nodemax
                    node_tally = jobsmax
                else:
                    newRun     = jobsmax - node_tally
                    node_tally = jobsmax

                if (newRun + countRun) > runlen:
                    newRun = runlen - countRun
                    if self.nodemax is not None:
                        node_tally = newRun*self.nodemax
                    else:               
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


            print 'Waiting 5 minutes'
            time.sleep(60*5)
            

        # Exited while loop
        print 'All jobs done'

        # Clean up workspaces for completed jobs
        for i,s in enumerate(jobid):
            s = s.strip('\n')
            if not self.errTally[i]:
                self.destroy_workspace(i,s)
                if self.nodemax is not None:
                    self.combine_files(i)

                    # Compress outfiles
                    self.compress(i,devnull)

        # Postprocessing done
        print 'Cleaned Up Worksapces'
        devnull.close()


    def check_for_errors(self,i,jobid):
        os.chdir(self.dirstring[i])  

        error = False 
        if self.nodemax is not None:
            for a in np.arange(self.nodemax):
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
    """ Create slurm scripts for running run_lidar_sampler.py """
    def __init__(self,args):

        self.Date      = isoparser(args.iso_t1)
        self.eddate    = isoparser(args.iso_t2)
        self.Dt        = timedelta(hours=args.DT_hours)

        if not os.path.exists(args.tmp):
            os.makedirs(args.tmp)

        self.cwd         = os.getcwd()
        self.slurm       = arg.slurm
        self.prep_config = args.prep_config
        self.tmp         = args.tmp


    def create_workdir(self):
        sdate = self.Date
        while sdate < self.enddate:
            # create directory

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
            source = ['lidar_sampler.py','run_lidar_sampler.py','sampling']
            for src in source:
                



            sdate += self.Dt

if __name__ == '__main__':
    
    #Defaults
    DT_hours = 24
    slurm    = 'run_lidar_sampler.j'
    tmp      = './tmp'

    parser = argparse.ArgumentParser()
    parser.add_argument("iso_t1",help='starting iso time')
    parser.add_argument("iso_t2",help='ending iso time')
    parser.add_argument("prep_config",
                        help="prep config filename")

    parser.add_argument('-D',"--DT_hours", default=DT_hours, type=int,
                        help="Timestep in hours for each file (default=%i)"%DT_hours)

    parser.add_argument('-s',"--slurm",default=slurm,
                        help="slurm script template (default=%s)"%slurm)           

    parser.add_argument('-t','--tmp',default=tmp,
                        help="temp directory (default=%s)"%tmp)

    args = parser.parse_args()

    Date    = isoparser(args.iso_t1)
    enddate = isoparser(args.iso_t2)
    Dt      = timedelta(hours=args.DT_hours)

    if not os.path.exists(args.tmp):
        os.makedirs(args.tmp)

    cwd = os.getcwd()
    while Date < enddate:
        edate = Date + Dt

        # copy template to temp
        outfile = '{}_{}.j'.format(args.slurm[:-2],Date.date())
        outpath = '{}/{}'.format(args.tmp,outfile)

        shutil.copyfile(args.slurm,outpath)

        f = open(outpath)
        text = []
        for l in f:
            text.append(l)

        iso1 = Date.isoformat()
        iso2 = edate.isoformat()
        newline = 'python -u run_lidar_sampler.py -v --nproc 6 --DT_hours {} {} {} {} >'.format(args.DT_hours,iso1,iso2,args.prep_config) + ' slurm_${SLURM_JOBID}_py.out\n'
        text[-2] = newline
        f.close()

        f = open(outpath,'w')
        for l in text:
            f.write(l)
        f.close()


        os.chdir(args.tmp)
        subprocess.call('sbatch {}'.format(outfile),shell=True)
        os.chdir(cwd)

        Date += Dt

