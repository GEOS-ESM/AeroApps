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
from   netCDF4         import Dataset
from   MAPL            import Config
from py_leo_vlidort.vlidort import get_chd

mr =  [1.396,1.362,1.349,1.345,1.339,1.335,1.334,1.333,1.332,1.331,1.329,1.326,
      1.323,1.318,1.312,1.306,1.292,1.261]

mr_ch = [200,250,300,337,400,488,515,550,633,694,860,1060,1300,1536,1800,2000,2250,2500]

mr = np.array(mr)
mr_ch = np.array(mr_ch)

mr_ACCP = 1.333

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

                        # Clean up workspaces
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


            print 'Waiting 30 minutes'
            time.sleep(60*30)
            

        # Exited while loop
        print 'All jobs done'

        # Compress Files
        for i,s in enumerate(jobid):
            s = s.strip('\n')
            if not self.errTally[i]:
                # Compress outfiles
                self.compress(i,devnull)

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
    """ Create slurm scripts for running pace_vlidort_lc.py """
    def __init__(self,args):
        self.Date    = isoparser(args.iso_t1)
        self.enddate   = isoparser(args.iso_t2)
        self.Dt        = timedelta(hours=args.DT_hours)        
        self.dt        = timedelta(hours=args.dt_hours)

        if not os.path.exists(args.tmp):
            os.makedirs(args.tmp)

        self.cwd           = os.getcwd()
        self.track_pcf   = args.track_pcf
        self.orbit_pcf   = args.orbit_pcf
        self.inst_pcf    = args.inst_pcf
        self.DT_hours    = args.DT_hours
        self.albedoType  = args.albedoType
        self.rcFile      = args.rcfile
        self.dryrun      = args.dryrun
        self.verbose     = args.verbose

        self.slurm       = args.slurm
        self.tmp         = args.tmp
        self.profile     = args.profile
        self.planeparallel = args.planeparallel

        # Parse prep config files
        # --------------------------

        # figure out channels from instfile
        cf             = Config(args.inst_pcf,delim=' = ')
        self.instname       = cf('instname')
        channels = cf('channels')
        if ',' in channels:
            channels = channels.split(',')
        else:
            channels = [channels]
        self.channels = np.array(channels).astype(int)
        self.nch = len(channels)

        # orbit name
        cf             = Config(args.orbit_pcf,delim=' = ')
        self.orbitname      = cf('orbitname')
        self.ORBITNAME      = self.orbitname.upper()

        # filename templates
        cf             = Config(args.track_pcf,delim=' = ')
        self.inTemplate     = cf('inDir')     + '/' + cf('inFile')
        self.outTemplate    = cf('outDir')    + '/' + cf('outFile')

        try:
            brdfTemplate = cf('brdfDir') + '/' + cf('brdfFile')
        except:
            brdfTemplate = None

        try:
            ndviTemplate   = cf('ndviDir')   + '/' + cf('ndviFile')
            lcTemplate     = cf('lcDir')     + '/' + cf('lcFile')
        except:
            ndviTemplate = None
            lcTemplate   = None

        try:
            lerTemplate    = cf('lerDir')    + '/' + cf('lerFile')
        except:
            lerTemplate    = None

        self.brdfTemplate = brdfTemplate
        self.ndviTemplate = ndviTemplate
        self.lcTemplate = lcTemplate
        self.lerTemplate = lerTemplate

        # create working directories
        self.create_workdir()
        self.errTally    = np.ones(len(self.dirstring)).astype(bool)


    def create_workdir(self):
        sdate = self.Date
        self.dirstring = []
        self.outfilelist = []
        while sdate < self.enddate:
            # create directory     
            workpath = '{}/{}'.format(self.tmp,sdate.isoformat())
            if os.path.exists(workpath):
                shutil.rmtree(workpath)

            os.makedirs(workpath)

            # Copy over slurm script    
            outfile = '{}/{}'.format(workpath,self.slurm)
            shutil.copyfile(self.slurm,outfile)

            # link over some scripts and files
            source = ['accp_polar_vlidort.x','clean_mem.sh']
            for src in source:
                os.symlink('{}/{}'.format(self.cwd,src),'{}/{}'.format(workpath,src))
            
            for ch in self.channels:
                subenddate = sdate + self.dt
                if subenddate > self.enddate:
                    subenddate = self.enddate
                subdate = sdate
                while subdate < subenddate:
                    # Create rc File
                    self.write_rc(workpath,subdate,ch)
                    
                    subdate += self.Dt

                # Copy over and edit Aod_EOS.rc
                outfile = '{}/Aod_EOS.{}.rc'.format(workpath,ch)
                shutil.copyfile(self.rcFile,outfile)
                self.edit_AODrc(outfile,ch)

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


    def write_rc(self,outpath,date,channel):
        nymd  = str(date.date()).replace('-','')
        year  = str(date.year)
        month = str(date.month).zfill(2)
        day   = str(date.day).zfill(2)
        hour  = str(date.hour).zfill(2)
        orbitname = self.orbitname
        ORBITNAME = self.ORBITNAME
        instname  = self.instname

        inFile     = self.inTemplate.replace('%year',year).replace('%month',month).replace('%day',day).replace('%nymd',nymd).replace('%hour',hour).replace('%orbitname',orbitname).replace('%ORBITNAME',ORBITNAME)
        outFile    = self.outTemplate.replace('%year',year).replace('%month',month).replace('%day',day).replace('%nymd',nymd).replace('%hour',hour).replace('%chd',get_chd(channel)).replace('%orbitname',orbitname).replace('%ORBITNAME',ORBITNAME).replace('%instname',instname)
        
        # create oufile directory
        if not os.path.exists(os.path.dirname(outFile)):
            os.makedirs(os.path.dirname(outFile))


        if self.brdfTemplate is None:
            brdfFile = None
        else:
            brdfFile = self.brdfTemplate.replace('%year',year).replace('%month',month).replace('%day',day).replace('%nymd',nymd).replace('%hour',hour).replace('%orbitname',orbitname).replace('%ORBITNAME',ORBITNAME)

        if self.ndviTemplate is None:
            ndviFile = None
            lcFile   = None
        else:
            ndviFile   = self.ndviTemplate.replace('%year',year).replace('%month',month).replace('%day',day).replace('%nymd',nymd).replace('%hour',hour).replace('%orbitname',orbitname).replace('%ORBITNAME',ORBITNAME)
            lcFile     = self.lcTemplate.replace('%year',year).replace('%month',month).replace('%day',day).replace('%nymd',nymd).replace('%hour',hour).replace('%orbitname',orbitname).replace('%ORBITNAME',ORBITNAME)

        if self.lerTemplate is None:
            lerFile = None
        else:
            lerFile   = self.lerTemplate.replace('%year',year).replace('%month',month).replace('%day',day).replace('%nymd',nymd).replace('%hour',hour).replace('%orbitname',orbitname).replace('%ORBITNAME',ORBITNAME)

        # figure out albedoType keyword
        if self.albedoType is 'BRDF':
            watermodel = 'CX'
            watermr    = mr_ACCP
            if channel <= 388:
                landname = 'LER'
                landmodel = 'LAMBERTIAN'
            # HYBRID
            elif (channel > 388.) and (channel < 470):
                landname = 'LER-MCD43C'
                landmodel = 'RTLS-HYBRID'
            else:
                landname = 'MCD43C'
                landmodel = 'RTLS'
        elif self.albedoType is 'BPDF':
            watermodel = 'GissCX'
            watermr = mr_ACCP
            if channel <= 388:
                landname = 'LER'
                landmodel = 'LAMBERTIAN'
            # HYBRID
            elif (channel > 388.) and (channel < 470):
                landname = 'LER-MCD43C-MAIGNAN'
                landmodel = 'RTLS-HYBRID-BPDF'
            else:
                landname = 'MCD43C-MAIGNAN'
                landmodel = 'RTLS-BPDF'


        filename = '{}/accp_polar_vlidort.{}.{}.rc'.format(outpath,date.isoformat(),channel)
        f = open(filename,'w')
        text = []

        # META STUFF
        nymd = date.strftime('%Y%m%d')
        hms  = date.strftime('%H%M00')

        newline  = 'DATE: {}\n'.format(nymd)
        text.append(newline)

        newline  = 'TIME: {}\n'.format(hms)
        text.append(newline)

        newline  = 'INSTNAME: {}-{}\n'.format(self.orbitname,self.instname)
        text.append(newline)

        newline  = 'CHANNELS: {}\n'.format(channel)
        text.append(newline)

        newline  = 'RCFILE: Aod_EOS.{}.rc\n'.format(channel)
        text.append(newline)

        newline  = 'SCALAR: false\n'
        text.append(newline)

        newline = 'PLANE_PARALLEL: {}\n'.format(self.planeparallel)
        text.append(newline)
        text.append('\n')

        # LAND STUFF

        newline = 'LANDNAME: {}\n'.format(landname)
        text.append(newline)
        newline = 'LANDMODEL: {}\n'.format(landmodel)
        text.append(newline)
        text.append('\n')

        # LER
        newline = 'LANDBAND_C_LER: 354 388\n'
        text.append(newline)
        newline = 'LER_file: {}\n'.format(lerFile)
        text.append(newline)
        text.append('\n')

        # BRDF
        newline = 'LANDBAND_C_BRDF: 470 550 650 850 1200 1600 2100\n'
        text.append(newline)
        newline = 'BRDF_file: {}\n'.format(brdfFile)
        text.append(newline)
        text.append('\n')

        # Land Cover
        newline = 'BPDF_file: {}\n'.format(lcFile)
        text.append(newline)

        # NDVI File
        newline = 'NDVI_file: {}\n'.format(ndviFile)
        text.append(newline)
        text.append('\n')

        # WATER STUFF
        newline = 'WATERMODEL: {}\n'.format(watermodel)
        text.append(newline)

        # #refractive index
        # if self.mr_spectral:
        #     mruse = np.interp(float(ch),mr_ch,mr)
        # else:
        #     mruse = mr_ACCP
        newline = 'WATERMR: {}\n'.format(watermr)
        text.append(newline)
        text.append('\n')


        # INPUT FILENAMES
        AER_file = inFile.replace('%col','aer_Nv')
        newline = 'AER_file: {}\n'.format(AER_file)
        text.append(newline)

        MET_file = inFile.replace('%col','met_Nv')
        newline = 'MET_file: {}\n'.format(MET_file)
        text.append(newline)

        ANG_file = inFile.replace('%col',self.instname)
        newline = 'ANG_file: {}\n'.format(ANG_file)
        text.append(newline)

        INV_file = inFile.replace('%col','asm_Nx')
        newline = 'INV_file: {}\n'.format(INV_file)
        text.append(newline)

        newline = 'OUT_file: {}\n'.format(outFile)
        text.append(newline)
        self.outfilelist.append(outFile)

        # write text to file
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

            os.remove(self.slurm)
            os.remove(self.rcFile)
            os.remove('accp_polar_vlidort.rc')

        # remove symlinks
        source = ['accp_polar_vlidort.x','clean_mem.sh']
        for src in source:
            os.remove(src)

        os.chdir(self.cwd)
        if self.profile is False:
            os.rmdir(self.dirstring[i])

    def compress(self,i,devnull):
        outfile = self.outfilelist[i]

        cmd = "$BASEDIR/Linux/bin/ncks --hst --abc -O -4 -L 2 --cnk_plc=g2d --cnk_dmn time,3600 --cnk_dmn view_angles,10 --cnk_dmn lev,72 {} {}"
        newcmd = cmd.format(outfile,outfile)
        
        stat = subprocess.call(newcmd, shell=True, stdout=devnull)
        





if __name__ == '__main__':
    
    #Defaults
    DT_hours  = 1
    dt_hours  = 48
    slurm     = 'run_accp_polar_vlidort_singlenode.j'
    tmp       = '/discover/nobackup/projects/gmao/osse2/pub/c1440_NR/OBS/A-CCP/workdir/vlidort'
    rcFile   = 'Aod_EOS.rc'
    albedoType = 'BRDF'    

    parser = argparse.ArgumentParser()
    parser.add_argument("iso_t1",help='starting iso time')
    parser.add_argument("iso_t2",help='ending iso time')

    parser.add_argument("track_pcf",
                        help="prep config file with orbit variables")

    parser.add_argument("orbit_pcf",
                        help="prep config file with orbit variables")

    parser.add_argument("inst_pcf",
                        help="prep config file with instrument variables")

    parser.add_argument('-D',"--DT_hours", default=DT_hours, type=int,
                        help="Timestep in hours for each file (default=%i)"%DT_hours)

    parser.add_argument('-d',"--dt_hours", default=dt_hours, type=int,
                        help="Timestep in hours for each job (default=%i)"%dt_hours)

    parser.add_argument("-a","--albedoType", default=albedoType,
                        help="albedo type keyword. either BRDF of BPDF (default=%s)"%albedoType)

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
    args.planeparallel = True

    workspace = WORKSPACE(args)

    if not args.dryrun:
        workspace.handle_jobs()

