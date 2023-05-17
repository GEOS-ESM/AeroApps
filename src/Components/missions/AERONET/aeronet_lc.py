#!/usr/bin/env python3

"""
    Wrapper to submit aeronet_lc.j jobs to sbatch
"""

import os
import subprocess
import shutil
from   datetime        import datetime, timedelta
from   dateutil.parser import parse         as isoparser
import argparse
import numpy           as np
import time

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
                        print('Jobid ',s,' in ',self.dirstring[i],' exited with errors')

            # finished checking up on all the jobs
            # Remove finished jobs from the currently working list
            if len(finishedJobs) != 0:
                print('deleting finishedJobs',finishedJobs,jobid[workingJobs[finishedJobs]])
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


            print('Waiting 20 minutes')
            time.sleep(60*20)
            

        # Exited while loop
        print('All jobs done')


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
    """ Create slurm scripts for running aeronet_vlidort_lc.py """
    def __init__(self,args):
        Date    = isoparser(args.iso_t1)
        enddate = isoparser(args.iso_t2)
        Dt      = timedelta(hours=args.DT_hours)
        Djobs   = timedelta(hours=args.DT_jobs)

        if not os.path.exists(args.tmp):
            os.makedirs(args.tmp)

        self.cwd = os.getcwd()
        self.runfile = args.slurm
        self.profile = args.profile

        self.dirstring = []
        while Date <= enddate:
            edate = Date + Djobs - Dt
            if edate > enddate: edate = enddate

            # copy template to temp
            outpath = '{}/{}'.format(args.tmp,Date.isoformat())
            if os.path.exists(outpath):
                shutil.rmtree(outpath)
            os.makedirs(outpath)    
            outfile = '{}/{}'.format(outpath,self.runfile)
            shutil.copyfile(self.runfile,outfile)

            self.dirstring.append(outpath)

            # Read file first - this is the copied template
            f = open(outfile)
            text = []
            for l in f:
                text.append(l)
            f.close()

            # Edit text
            iso1 = Date.isoformat()
            iso2 = edate.isoformat()

            newline  = 'setenv START {}\n'.format(iso1)
            text[18] = newline

            newline  = 'setenv END  {}\n'.format(iso2)
            text[19] = newline

            # Create options string
            options = "-D {}".format(args.DT_hours)
            if args.verbose:
                options += ' --verbose'
            if args.dryrun:
                options += ' --dryrun'
            if args.channels is not None:
                options += ' --channels {}'.format(' '.join(args.channels))
            if args.extch is not None:
                options += ' --extch {}'.format(' '.join(args.extch))
            if args.norad:
                options += ' --norad'
            if args.noext:
                options += ' --noext'

            newline  = 'setenv OPTIONS "{}"\n'.format(options)
            text[20] = newline

            newline  = 'setenv AEROBIN {}\n'.format(outpath)
            text[21] = newline

            #Write edited text back to file
            f = open(outfile,'w')
            for l in text:
                f.write(l)
            f.close()

            # link over some scripts and files
            source = ['nccs','Aod_EOS.rc','Aod_EOS.440.rc','Aod_EOS.870.rc',
                      'run_aeronet_vlidort.py','aeronet_vlidort_lc.py','ExtData',
                      'Chem_MieRegistry.rc']
            for src in source:
                os.symlink('{}/{}'.format(self.cwd,src),'{}/{}'.format(outpath,src))


            Date += Djobs


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
    DT_hours = 1
    DT_jobs  = 1    #hours
    slurm    = 'aeronet_lc.j'
    tmp      = '/discover/nobackup/projects/gmao/osse2/pub/c1440_NR/OBS/AERONET/workdir'

    parser = argparse.ArgumentParser()
    parser.add_argument("iso_t1",help='starting iso time')
    parser.add_argument("iso_t2",help='ending iso time')


    parser.add_argument('-D',"--DT_hours", default=DT_hours, type=int,
                        help="Timestep in hours for each file (default=%i)"%DT_hours)

    parser.add_argument("--DT_jobs", default=DT_jobs, type=int,
                        help="Timestep in days for each job (default=%i)"%DT_jobs)    

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
                        help="channels to get BOA radiance (default=None)")                            

    parser.add_argument("-e","--extch", default=None,nargs='+',
                        help="channels to run extinction sampler (default=None)")  

    parser.add_argument("--norad",action="store_true",
                        help="No radiance calculations (default=False).")

    parser.add_argument("--noext",action="store_true",
                        help="No extinctions calculations (default=False).")


    args = parser.parse_args()

    workspace = WORKSPACE(args)

    # Submit and monitor jobs
    workspace.handle_jobs()
