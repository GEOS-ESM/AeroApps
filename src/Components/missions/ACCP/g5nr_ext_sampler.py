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
import numpy as np
import time
from MAPL     import Config

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
            jobid = np.append(jobid,subprocess.check_output(['qsub',self.slurm]))
        os.chdir(self.cwd)

        # launch subprocess that will monitor queue

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
                    jobid = np.append(jobid,subprocess.check_output(['qsub',self.slurm]))

                os.chdir(self.cwd)
                countRun = countRun + newRun
                stat = subprocess.call(['qstat -u pcastell'], shell=True, stdout=devnull)


            print 'Waiting 5 minutes'
            time.sleep(60*1)
            

        # Exited while loop
        print 'All jobs done'

        # Clean up workspaces for completed jobs
        for i,s in enumerate(jobid):
            s = s.strip('\n')
            if not self.errTally[i]:
                self.destroy_workspace(i,s)

        # Postprocessing done
        print 'Cleaned Up Worksapces'
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
        self.DT_hours  = args.DT_hours
        self.dt        = timedelta(hours=args.dt_hours)

        if not os.path.exists(args.tmp):
            os.makedirs(args.tmp)

        self.cwd         = os.getcwd()
        self.slurm       = args.slurm
        self.prep_config = args.prep_config
        self.tmp         = args.tmp
        self.profile     = args.profile
        self.instname    = args.instname
        self.channels    = args.channels
        self.rcFile      = args.rcFile
        self.extData     = args.extData
        self.inDir       = args.inDir
        self.outDir      = args.outDir
        self.case        = args.case

        # create working directories
        self.create_workdir()
        self.errTally    = np.ones(len(self.dirstring)).astype(bool)

        # modify slurm scripts
        self.edit_slurm()


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

                # Copy over Aod_EOS.rc
                os.makedirs(workpath+'/rc')
                outfile = '{}/{}'.format(workpath,self.rcFile)
                shutil.copyfile(self.rcFile,outfile)
                self.edit_AODrc(outfile,ch)            

                #cp over needed extdata
                source = ['Chem_MieRegistry.rc','run_ext_sampler.py']
                for src in source:
                    shutil.copyfile('{}/{}'.format(self.cwd,src),'{}/{}'.format(workpath,src))
               
                source = ['ExtData'] #,self.extData]
                for src in source:
                    os.symlink('{}/{}'.format(self.cwd,src),'{}/{}'.format(workpath,src))                                        

                self.dirstring.append(workpath)
            sdate += self.dt

    def edit_AODrc(self,filename,ch):
        f = open(filename)
        # read the file first
        text = []
        for l in f:
            text.append(l)   

        f.close()    

        # Edit text
        for i,l in enumerate(text):
            if 'r_channels:' in l[0:11]:
                newline = "r_channels: {}e-6".format(ch*0.001)
                text[i] = newline

        # Write edited text back to file
        f = open(filename,'w')
        for l in text:
            f.write(l)

        f.close()

    def edit_slurm(self):

        for workpath in self.dirstring:
            isodate, ch = os.path.basename(workpath).split('.')
            sdate = isoparser(isodate)
            edate = sdate + self.dt
            isoedate = edate.isoformat()

            # get inFile
            inFile = self.inDir + '/Y{}/M{}/D{}/' + self.instname + '-g5nr.lb2.aer_Nv.{}.nc4'

            # outFile
            outFile = self.outDir + '/Y{}/M{}/D{}/'+ self.instname +'-g5nr.lc.ext.{}.{}nm.nc4'
            
            # read file first
            outpath = '{}/{}'.format(workpath,self.slurm)
            f = open(outpath)

            text = []
            for l in f:
                text.append(l)


            # replace one line
            if len(self.case) == 0:
                command = 'python -u ./run_ext_sampler.py --DT_hours {} --rc {} {} {} '.format(self.DT_hours,self.rcFile,isodate,isoedate)
            else:
                options = ''
                for case in self.case:
                    options = options + ' --' + case
                command = 'python -u ./run_ext_sampler.py {} --DT_hours {} --rc {} {} {} '.format(options,self.DT_hours,self.rcFile,isodate,isoedate)

            command = command + " '" + inFile + "' "
            command = command + " '" + outFile + "' "
            command = command + ch
            endcommand = ' > slurm_${SLURM_JOBID}_py.out\n'
            newline = command  + endcommand
            text[-2] = newline
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

        # remove symlinks
        source = ['ExtData','Chem_MieRegistry.rc','run_ext_sampler.py']
        for src in source:
            os.remove(src)
	    # remove rc directory
	    #shutil.rmtree('rc')
        #shutil.rmtree(self.extData)

        os.chdir(self.cwd)
        if self.profile is False:
            #os.rmdir(self.dirstring[i])
            shutil.rmtree(self.dirstring[i])


if __name__ == '__main__':
    
    #Defaults
    DT_hours = 1
    dt_hours = 24
    slurm    = 'run_ext_sampler.j'
    tmp      = '/discover/nobackup/projects/gmao/osse2/pub/c1440_NR/OBS/A-CCP/workdir'

    parser = argparse.ArgumentParser()
    parser.add_argument("iso_t1",help='starting iso time')
    parser.add_argument("iso_t2",help='ending iso time')
    parser.add_argument("prep_config",
                        help="prep config filename")

    parser.add_argument('-D',"--DT_hours", default=DT_hours, type=int,
                        help="Timestep in hours for each file (default=%i)"%DT_hours)

    parser.add_argument('-d',"--dt_hours", default=dt_hours, type=int,
                        help="Timestep in hours for job (default=%i)"%dt_hours)    

    parser.add_argument('-s',"--slurm",default=slurm,
                        help="slurm script template (default=%s)"%slurm)           

    parser.add_argument('-t','--tmp',default=tmp,
                        help="temp directory (default=%s)"%tmp)

    parser.add_argument("-r", "--dryrun",action="store_true",
                        help="do a dry run (default=False).") 

    parser.add_argument("-p", "--profile",action="store_true",
                        help="Don't cleanup slurm files (default=False).")    


    args = parser.parse_args()

    # Parse prep config
    # -----------------
    cf = Config(args.prep_config,delim=' = ')
    args.instname = cf('INSTNAME')
    if ',' in cf('CHANNELS'):
        args.channels = np.array(cf('CHANNELS').split(',')).astype(int)
    else:
        args.channels = np.array([int(cf('CHANNELS'))])
    args.rcFile   = cf('RCFILE')
    args.extData  = cf('EXTDATA')
    args.inDir    = cf('INDIR')
    args.outDir   = cf('OUTDIR')

    if ',' in cf('CASE'):
        args.case     = cf('CASE').split(',')
    else:
        args.case     = list(cf('CASE'))

    if 'all' in args.case:
        args.case.remove('all')


    workspace = WORKSPACE(args)

    if not args.dryrun:
        workspace.handle_jobs()
