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
from   glob            import glob
from MAPL.constants import MAPL_UNDEF


mr =  [1.396,1.362,1.349,1.345,1.339,1.335,1.334,1.333,1.332,1.331,1.329,1.326,
      1.323,1.318,1.312,1.306,1.292,1.261]

mr_ch = [200,250,300,337,400,488,515,550,633,694,860,1060,1300,1536,1800,2000,2250,2500]

mr = np.array(mr)
mr_ch = np.array(mr_ch)

mr_OCI = 1.334

nlev = 72

#---
class JOBS(object):
    def handle_jobs(self):
        jobsmax = 50
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
            result = subprocess.check_output(['sbatch',self.runfile])
            jobid = np.append(jobid,result.split()[-1])
        os.chdir(self.cwd)

        # launch subprocess that will monitor queue and do file merging

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
                    countDone += 1
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
                    result = subprocess.check_output(['sbatch',self.runfile])
                    jobid = np.append(jobid,result.split()[-1])

                os.chdir(self.cwd)
                countRun = countRun + newRun

            # check if all the jobs are finished
            if countDone == runlen:
                stat = 1
            else:
                print 'Waiting 30 minutes'
                time.sleep(60*3)
            

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
    """ Create slurm scripts for running pace_vlidort_lc.py """
    def __init__(self,args):
        self.Date    = isoparser(args.iso_t1)

        if not os.path.exists(args.tmp):
            os.makedirs(args.tmp)

        self.cwd           = os.getcwd()
        self.runfile       = args.slurm
        self.runfile_m     = args.slurm_m
        self.profile       = args.profile
        self.rootdir       = args.rootdir
        self.verbose       = args.verbose
        self.run_name      = args.run_name
        self.bpdf_name     = args.bpdf_name
        self.water_name    = args.water_name
        self.rtls_name     = args.rtls_name
        self.write_add     = args.write_add
        self.write_aer     = args.write_aer
        self.write_cld     = args.write_cld
        self.planeparallel = args.planeparallel
        self.do_diagnostic_only = args.do_diagnostic_only
        self.do_sleave_adjust = args.do_sleave_adjust
        self.do_sleave_iso = args.do_sleave_iso
        self.do_single_xtrack = args.do_single_xtrack
        self.cloudfree     = args.cloudfree
        self.aerosolfree   = args.aerosolfree
        self.gasfree       = args.gasfree
        self.mr_spectral   = args.mr_spectral
        self.version       = args.version

        self.ALPHA_file    = args.ALPHA_file

        self.dirstring = []
        self.outfilelist = []
        self.addfilelist = []
        self.aerfilelist = []
        self.cldfilelist = []

        if (args.channels is None) and (args.ichannels is None):
            self.get_channels()
        elif args.channels is not None:
            if ',' in args.channels:
                makelist=lambda s: map(int, s.split(","))
                self.channels = makelist(args.channels)
            elif ':' in args.channels:
                makelist=lambda s: map(int, s.split(":"))
                start,stop,delta = makelist(args.channels)
                self.channels = range(start,stop+delta,delta)
            else:
                self.channels = [float(args.channels)]

            self.get_ichannels()
        elif args.ichannels is not None:
            if ',' in args.ichannels:
                makelist=lambda s: map(int, s.split(","))
                self.ichannels = makelist(args.ichannels)
            elif ':' in args.ichannels:
                makelist=lambda s: map(int, s.split(":"))
                start,stop,delta = makelist(args.ichannels)
                self.ichannels = range(start,stop+delta,delta)
            else:
                self.ichannels = [int(args.ichannels)]

            self.get_channels(ichannels=self.ichannels)

        if int(args.nodemax) > 1 : 
            self.nodemax = int(args.nodemax)
        else:
            self.nodemax = None

        for ich,ch,nbin in zip(self.ichannels,self.channels,self.nbins):
            fch = "{:.2f}".format(ch).replace('.','d')
            outpath = '{}/{}-{}.{}'.format(args.tmp,str(ich).zfill(4),self.Date.isoformat(),fch)

            # Copy over some files to working temp directory
            self.create_workdir(outpath,ch,nbin)

            # # Create pcf file
            self.write_rc(outpath,ich,ch)

            self.dirstring.append(outpath)


        self.errTally    = np.ones(len(self.dirstring)).astype(bool)

    def create_workdir(self,outpath,ch,nbin):
        # create directory     
        if os.path.exists(outpath):
            shutil.rmtree(outpath)
        os.makedirs(outpath)
        os.makedirs(outpath+'/oci_tables')
        alphadir = os.path.dirname(self.ALPHA_file).split('/')[-1]
        os.makedirs(outpath+'/'+alphadir)

        # Copy over slurm script    
        outfile = '{}/{}'.format(outpath,self.runfile)
        if nbin > 1:
            shutil.copyfile(self.runfile_m,outfile)
        else:
            shutil.copyfile(self.runfile,outfile)

        if self.nodemax is not None:
            self.edit_runFile(outfile)

        # Copy over Aod_EOS.rc
        outfile = '{}/{}'.format(outpath,'Aod_EOS.rc')
        shutil.copyfile('Aod_EOS.rc',outfile)
        self.edit_AODrc(outfile,ch)

        # link over some scripts and files
        source = ['ExtData','Chem_MieRegistry.rc',
                  'ExtDataCloud','ExtDataOsku','clean_mem.sh',
                   self.ALPHA_file]
        for src in source:
            os.symlink('{}/{}'.format(self.cwd,src),'{}/{}'.format(outpath,src))

        os.symlink('{}/{}'.format(self.cwd,'pace_vlidort_gasabs_multinode_v0.x'),'{}/{}'.format(outpath,'pace_vlidort_gasabs_multinode.x'))

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

    def edit_runFile(self,filename):
        f = open(filename)
        # read the file first
        text = []
        for l in f:
            text.append(l)

        f.close()

        # Edit text
        for i,l in enumerate(text):
            if (l[0:15] == '#SBATCH --array'):
                newline = '#SBATCH --array=1-'+str(self.nodemax)+'\n'
                text[i] = newline

        # Write edited text back to file
        f = open(filename,'w')
        for l in text:
            f.write(l)

        f.close()


    def write_rc(self,outpath,ich_d,ch):
        ich       = str(ich_d).zfill(4)
        YMDdir    = self.Date.strftime('Y%Y/M%m/D%d')
        LbDir     = '{}/LevelB/{}'.format(self.rootdir,YMDdir)
        LcDir     = '{}/LevelC/{}'.format(self.rootdir,YMDdir)
        LcDirCh   = '{}/LevelC/{}/v{}/channel'.format(self.rootdir,YMDdir,self.version)
        surfDir   = '{}/LevelB/surface'.format(self.rootdir)

        pYMDdir   = self.Date.strftime('Y2020/M%m/D%d')
        L1bDir    = '{}/L1B/{}'.format(self.rootdir,pYMDdir)


        filename = '{}/pace_vlidort.rc'.format(outpath)
        f = open(filename,'w')
        text = []

        # META STUFF
        nymd = self.Date.strftime('%Y%m%d')
        hms  = self.Date.strftime('%H%M00')

        newline  = 'DATE: {}\n'.format(nymd)
        text.append(newline)

        newline  = 'TIME: {}\n'.format(hms)
        text.append(newline)

        newline  = 'INSTNAME: pace\n'
        text.append(newline)

        newline  = 'CHANNELS: {}\n'.format(ch)
        text.append(newline)

        newline  = 'IOCI: {}\n'.format(ich_d+1)
        text.append(newline)

        newline  = 'SCALAR: false\n'
        text.append(newline)

        newline = 'PLANE_PARALLEL: {}\n'.format(self.planeparallel)
        text.append(newline)

        newline = 'DO_SLEAVE_ADJUST: {}\n'.format(self.do_sleave_adjust)
        text.append(newline)

        newline = 'DO_SLEAVE_ISO: {}\n'.format(self.do_sleave_iso)
        text.append(newline)

        newline = 'DO_SINGLE_XTRACK: {}\n'.format(self.do_single_xtrack)
        text.append(newline)

        newline = 'DO_DIAGNOSTIC_ONLY: {}\n'.format(self.do_diagnostic_only)
        text.append(newline)

        newline = 'CLOUD_FREE: {}\n'.format(self.cloudfree)
        text.append(newline)

        newline = 'AEROSOL_FREE: {}\n'.format(self.aerosolfree)
        text.append(newline)

        newline = 'GAS_FREE: {}\n'.format(self.gasfree)
        text.append(newline)

        if self.nodemax is not None:
            newline = 'NODEMAX: {}\n'.format(self.nodemax)
        else:
            newline = 'NODEMAX: 1\n'

        text.append(newline)

        if self.write_add:
            newline = 'ADDITIONAL_OUTPUT: true\n'
        else:
            newline = 'ADDITIONAL_OUTPUT: false\n'
        text.append(newline)
        newline = 'VERSION: {}\n'.format(self.version)
        text.append(newline)
        newline = 'LAYOUT: None\n'
        text.append(newline)
        text.append('\n')

        # OCI SPECIFIC STUFF
        newline = 'ALPHA_file: {}\n'.format(self.ALPHA_file)
        text.append(newline)
        text.append('\n')

        # LAND STUFF

        # LER
        if ch <= 388.0:
            newline = 'LANDNAME: LER\n'
            text.append(newline)
            newline = 'LANDMODEL: LAMBERTIAN\n'
            text.append(newline)
            newline = 'LANDBAND_C_LER: 354 388\n'
            text.append(newline)
            LER_file = '{}/OMI_LER/{}/pace-g5nr.lb.omi_ler.{}_{}.nc4'.format(surfDir,YMDdir,nymd,hms)
            newline = 'LER_file: {}\n'.format(LER_file)
            text.append(newline)
        # HYBRID
        elif (ch > 388.) and (ch < 470):
            if self.bpdf_name is None:
                newline = 'LANDNAME: LER-{}\n'.format(self.rtls_name)
            else:
                newline = 'LANDNAME: LER-{}-{}\n'.format(self.rtls_name,self.bpdf_name)                

            text.append(newline)
            if self.bpdf_name is None:
                newline = 'LANDMODEL: RTLS-HYBRID\n'
            else:
                newline = 'LANDMODEL: RTLS-HYBRID-BPDF\n'
            text.append(newline)
            newline = 'LANDBAND_C_LER: 354 388\n'
            text.append(newline)
            LER_file = '{}/OMI_LER/{}/pace-g5nr.lb.omi_ler.{}_{}.nc4'.format(surfDir,YMDdir,nymd,hms)
            newline = 'LER_file: {}\n'.format(LER_file)
            text.append(newline)
            newline = 'LANDBAND_C_BRDF: 470 550 650 850 1200 1600 2100\n'
            text.append(newline)
            BRDF_file = '{}/BRDF/MCD43C1/006/{}/pace-g5nr.lb.brdf.{}_{}.nc4'.format(surfDir,YMDdir,nymd,hms)
            newline = 'BRDF_file: {}\n'.format(BRDF_file)
            text.append(newline)

            if self.bpdf_name is not None:
                BPDF_file = '{}/BPDF/LAND_COVER/MCD12C1/051/{}/pace-g5nr.lb.land_cover.{}_{}.nc4'.format(surfDir,YMDdir,nymd,hms)
                newline = 'BPDF_file: {}\n'.format(BPDF_file)
                text.append(newline)
                NDVI_file = '{}/BPDF/NDVI/MYD13C2/006/pace-g5nr.lb.myd13c2.{}_{}.nc4'.format(surfDir,YMDdir,nymd,hms)
                newline = 'NDVI_file: {}\n'.format(NDVI_file)
                text.append(newline)

        else:
            if self.bpdf_name is None:
                newline = 'LANDNAME: {}\n'.format(self.rtls_name)
            else:
                newline = 'LANDNAME: {}-{}\n'.format(self.rtls_name,self.bpdf_name)                
            text.append(newline)
            if self.bpdf_name is None:
                newline = 'LANDMODEL: RTLS\n'
            else:
                newline = 'LANDMODEL: RTLS-BPDF\n'
            text.append(newline)
            newline = 'LANDBAND_C_BRDF: 470 550 650 850 1200 1600 2100\n'
            text.append(newline)
            BRDF_file = '{}/BRDF/MCD43C1/006/{}/pace-g5nr.lb.brdf.{}_{}.nc4'.format(surfDir,YMDdir,nymd,hms)
            newline = 'BRDF_file: {}\n'.format(BRDF_file)
            text.append(newline)

            if self.bpdf_name is not None:            
                BPDF_file = '{}/BPDF/LAND_COVER/MCD12C1/051/{}/pace-g5nr.lb.land_cover.{}_{}.nc4'.format(surfDir,YMDdir,nymd,hms)
                newline = 'BPDF_file: {}\n'.format(BPDF_file)
                text.append(newline)
                NDVI_file = '{}/BPDF/NDVI/MYD13C2/006/pace-g5nr.lb.myd13c2.{}_{}.nc4'.format(surfDir,YMDdir,nymd,hms)
                newline = 'NDVI_file: {}\n'.format(NDVI_file)
                text.append(newline)

        text.append('\n')
        # WATER STUFF
        if self.water_name is not None:
            newline = 'WATERNAME: {}\n'.format(self.water_name)
            text.append(newline)
            if self.water_name == 'NOBM-CX':
                WAT_file = '{}/SLEAVE/NOBM/{}/pace-g5nr.lb.sleave.{}_{}.nc4'.format(surfDir,YMDdir,nymd,hms)
                newline = 'WAT_file: {}\n'.format(WAT_file)
                text.append(newline)                
        else:
            if ch <= 775.0:
                newline = 'WATERNAME: NOBM-CX\n'
                text.append(newline)
                WAT_file = '{}/SLEAVE/NOBM/{}/pace-g5nr.lb.sleave.{}_{}.nc4'.format(surfDir,YMDdir,nymd,hms)
                newline = 'WAT_file: {}\n'.format(WAT_file)
                text.append(newline)

                if not self.do_sleave_iso:
                    FQ_file = '{}/SLEAVE/NOBM/{}/pace-g5nr.lb.fq.{}_{}.nc4'.format(surfDir,YMDdir,nymd,hms)
                    newline = 'FQ_file: {}\n'.format(FQ_file)
                    text.append(newline)
            else:
                newline = 'WATERNAME: CX\n'
                text.append(newline)

        newline = 'WATERMODEL: CX\n'
        text.append(newline)

        #refractive index
        if self.mr_spectral:
            mruse = np.interp(float(ch),mr_ch,mr)
        else:
            mruse = mr_OCI
        newline = 'WATERMR: {}\n'.format(mruse)
        text.append(newline)
        text.append('\n')

        # CLOUD STUFF
        newline = 'ICLDTABLE:ExtDataCloud/IceLegendreCoeffs.nc4\n'
        text.append(newline)
        newline = 'LCLDTABLE:ExtDataCloud/WaterLegendreCoeffs.nc4\n'
        text.append(newline)
        text.append('\n')

        # OTHER FILENAMES
        AER_file = '{}/pace-g5nr-std.lb.aer_Nv.{}_{}.nc4'.format(LbDir,nymd,hms)
        newline = 'AER_file: {}\n'.format(AER_file)
        text.append(newline)

        MET_file = '{}/pace-g5nr-std.lb.met_Nv.{}_{}.nc4'.format(LbDir,nymd,hms)
        newline = 'MET_file: {}\n'.format(MET_file)
        text.append(newline)

        CHM_file = '{}/pace-g5nr-std.lb.chm_Nv.{}_{}.nc4'.format(LbDir,nymd,hms)
        newline = 'CHM_file: {}\n'.format(CHM_file)
        text.append(newline)

        pDate = isoparser(self.Date.strftime('2020-%m-%dT%H:%M:00'))
        nyj = pDate.strftime('%Y%j')
        ANG_file = '{}/OCI{}{}.L1B_PACE.nc'.format(L1bDir,nyj,hms)
        newline = 'ANG_file: {}\n'.format(ANG_file)
        text.append(newline)

        INV_file = '{}/pace-g5nr.lb.asm_Nx.{}_{}.nc4'.format(LbDir,nymd,hms)
        newline = 'INV_file: {}\n'.format(INV_file)
        text.append(newline)

        fch = "{:.2f}".format(ch).replace('.','d')
        if self.run_name is None:
            OUT_file = '{}/pace-g5nr.lc.vlidort.{}_{}.{}-{}.nc4'.format(LcDirCh,nymd,hms,ich,fch)
        else:
            OUT_file = '{}/pace-g5nr.lc.vlidort.{}.{}_{}.{}-{}.nc4'.format(LcDirCh,self.run_name,nymd,hms,ich,fch)
        newline = 'OUT_file: {}\n'.format(OUT_file)
        text.append(newline)
        self.outfilelist.append(OUT_file)

        if self.write_add:
            if self.run_name is None:
                ADD_file = '{}/pace-g5nr.lc.add.{}_{}.{}-{}.nc4'.format(LcDirCh,nymd,hms,ich,fch)
            else:
                ADD_file = '{}/pace-g5nr.lc.add.{}.{}_{}.{}-{}.nc4'.format(LcDirCh,self.run_name,nymd,hms,ich,fch)
            newline = 'ADD_file: {}\n'.format(ADD_file)
            text.append(newline)
            self.addfilelist.append(ADD_file)

        CLD_file = '{}/pace-g5nr.TOTWPDF-GCOP-SKEWT.{}_{}.nc4'.format(LcDir,nymd,hms)
        newline = 'CLD_file: {}\n'.format(CLD_file)
        text.append(newline)

        if self.write_cld:
            newline = 'CLOUD_OUTPUT: true\n'
            text.append(newline)
            if self.run_name is None:
                CLDO_file = '{}/pace-g5nr.lc.cloud.{}_{}.{}-{}.nc4'.format(LcDirCh,nymd,hms,ich,fch)
            else:
                CLDO_file = '{}/pace-g5nr.lc.cloud.{}.{}_{}.{}-{}.nc4'.format(LcDirCh,self.run_name,nymd,hms,ich,fch)
            newline = 'CLDO_file: {}\n'.format(CLDO_file)
            text.append(newline)
            self.cldfilelist.append(CLDO_file)   

        if self.write_aer:
            newline = 'AEROSOL_OUTPUT: true\n'
            text.append(newline)
            if self.run_name is None:
                AERO_file = '{}/pace-g5nr.lc.aerosol.{}_{}.{}-{}.nc4'.format(LcDirCh,nymd,hms,ich,fch)
            else:
                AERO_file = '{}/pace-g5nr.lc.aerosol.{}.{}_{}.{}-{}.nc4'.format(LcDirCh,self.run_name,nymd,hms,ich,fch)
            newline = 'AERO_file: {}\n'.format(AERO_file)
            text.append(newline)
            self.aerfilelist.append(AERO_file)                        
        # write text to file
        for l in text:
            f.write(l)

        f.close()
# ---
    def get_channels(self,ichannels=None):
        nc = Dataset(self.ALPHA_file)
        self.channels  = nc.variables['channels'][:]
        self.nbins     = nc.variables['nbins'][:]
        nc.close()
        if ichannels is None:
            self.ichannels = np.arange(len(self.channels))
        else:
            self.channels = self.channels[ichannels]
            self.nbins    = self.nbins[ichannels]
# ---
    def get_ichannels(self):
        nc = Dataset(self.ALPHA_file)
        channels = nc.variables['channels'][:]
        nbins    = nc.variables['nbins'][:]
        nc.close()
        nch = len(channels)
        self.ichannels = []
        self.nbins     = []
        for ch in self.channels:
          ii = np.arange(nch)[channels==ch]
          if len(ii) == 1:
              self.ichannels.append(ii[0])
              self.nbins.append(nbins[ii[0]])
          else:
              stop = False
              for i in ii:
                  if not stop:
                      if i not in self.ichannels:
                          self.ichannels.append(i)
                          self.nbins.append(nbins[i])
                          stop = True
# ---
    def destroy_workspace(self,i,jobid):
        os.chdir(self.dirstring[i])

        if self.profile is False:

            if self.nodemax is not None:
                for a in np.arange(self.nodemax):
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

            os.remove(self.runfile)

        # remove symlinks
        source = ['Aod_EOS.rc','ExtData','ExtDataCloud','ExtDataOsku','pace_vlidort.rc',
                  'pace_vlidort_gasabs_multinode.x','Chem_MieRegistry.rc','clean_mem.sh']
        for src in source:
            os.remove(src)

        shutil.rmtree('oci_tables')
        alphadir = os.path.dirname(self.ALPHA_file).split('/')[-1]
        shutil.rmtree(alphadir)

        os.chdir(self.cwd)
        if self.profile is False:
            os.rmdir(self.dirstring[i])

    def compress(self,i,devnull):

        cmd = "$BASEDIR/Linux/bin/ncks --hst --no_abc -O -4 -L 2 --cnk_plc=g2d --cnk_dmn ccd_pixels,91 --cnk_dmn number_of_scans,144 --cnk_dmn lev,1 {} {}"

        if not self.do_diagnostic_only:
            outfile = self.outfilelist[i]
            newcmd = cmd.format(outfile,outfile)        
            stat = subprocess.call(newcmd, shell=True, stdout=devnull)

        if self.write_add:
            addfile = self.addfilelist[i]
            newcmd = cmd.format(addfile,addfile)
            stat = subprocess.call(newcmd, shell=True, stdout=devnull)

        if self.write_aer:
            aerfile = self.aerfilelist[i]
            newcmd = cmd.format(aerfile,aerfile)
            stat = subprocess.call(newcmd, shell=True, stdout=devnull)

        if self.write_cld:
            cldfile = self.cldfilelist[i]
            newcmd = cmd.format(cldfile,cldfile)
            stat = subprocess.call(newcmd, shell=True, stdout=devnull)
        
    def combine_files(self,i):
        outfile = self.outfilelist[i]
        filelist = []
        for a in np.arange(self.nodemax):
            a = a + 1
            filelist.append(outfile + '_' + str(a) )

        self.do_merge(outfile,filelist)

        if self.write_add:
            addfile = self.addfilelist[i]
            filelist = []
            for a in np.arange(self.nodemax):
                a = a + 1
                filelist.append(addfile + '_' + str(a) )
            
            self.do_merge(addfile,filelist)

        if self.write_aer:
            aerfile = self.aerfilelist[i]
            filelist = []
            for a in np.arange(self.nodemax):
                a = a + 1
                filelist.append(aerfile + '_' + str(a) )
            
            self.do_merge(aerfile,filelist)

        if self.write_cld:
            cldfile = self.cldfilelist[i]
            filelist = []
            for a in np.arange(self.nodemax):
                a = a + 1
                filelist.append(cldfile + '_' + str(a) )
            
            self.do_merge(cldfile,filelist)            


    def do_merge(self,mergedfile,filelist):
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



if __name__ == '__main__':
    
    #Defaults
    nodemax    = 1 
    slurm      = 'pace_lc_gas.j'
    slurm_m    = 'pace_lc_gas_multinode.j'
    rootdir    = '/discover/nobackup/projects/gmao/osse2/pub/c1440_NR/OBS/PACE'
    tmp        = '/discover/nobackup/projects/gmao/osse2/pub/c1440_NR/OBS/PACE/v0_workdir'
    rtls_name  = 'MCD43C'
    bpdf_name  = None # 'MAIGNAN'
    water_name = None
    ALPHA_file = 'alphaTable_v0/alpha_CK_Thuillier_o3.nc4'
    version    = 2.0

    parser = argparse.ArgumentParser()
    parser.add_argument("iso_t1",help='starting iso time')

    parser.add_argument("--rootdir",default=rootdir,
                        help="root directory for PACE data (default=%s)."%rootdir)       

    parser.add_argument('-s',"--slurm",default=slurm,
                        help="slurm script template (default=%s)"%slurm)           

    parser.add_argument("--slurm_m",default=slurm_m,
                        help="multinode slurm script template (default=%s)"%slurm_m)

    parser.add_argument('-t','--tmp',default=tmp,
                        help="temp directory (default=%s)"%tmp)

    parser.add_argument("-p", "--profile",action="store_true",
                        help="Don't cleanup slurm files (default=False).")    

    parser.add_argument("-v", "--verbose",action="store_true",
                        help="Verbose mode (default=False).")

    parser.add_argument("-r", "--dryrun",action="store_true",
                        help="do a dry run (default=False).")   

    parser.add_argument("-c","--channels", default=None,
                        help="channels to get TOA radiance. Can be a list (1,2) or a range (1:2:1) (default=None - read in from PACE L1B File)")                            

    parser.add_argument("--ichannels", default=None,
                        help="index of channels to get TOA radiance. Can be a list (1,2) or a range (1:2:1) (default=None - read in from PACE L1B File)")

    parser.add_argument("-e","--extch", default=None,type=lambda s: map(int, s.split(",")),
                        help="channels to run extinction sampler (default=None - read in from PACE L1B File)")  

    parser.add_argument('--version',default=version,
                        help="data version numbering (default={})".format(version))

    parser.add_argument('--run_name',default=None,
                        help="Optional version naming for outfiles (default=None)")    

    parser.add_argument('--bpdf_name',default=bpdf_name,
                        help="BPDF model name (e.g. MAIGNAN) (default=%s)"%bpdf_name)   

    parser.add_argument('--water_name',default=water_name,
                        help="Ocean reflectance model name (e.g. CX-NOBM) (default=%s, set by wavelength)"%water_name)                                                  

    parser.add_argument('--rtls_name',default=rtls_name,
                        help="RTLS data name (e.g. MCD43C) (default=%s)"%rtls_name)     

    parser.add_argument('--ALPHA_file',default=ALPHA_file,
                        help="gas cross section file (default=%s)"%ALPHA_file)

    parser.add_argument("--withcloud",action="store_true",
                        help="Include cloud effects (default=False).")    

    parser.add_argument("--withaerosol",action="store_true",
                        help="Include aerosol effects (default=False).")   

    parser.add_argument("--withgas",action="store_true",
                        help="Include trace gas absorption effects (default=False).")

    parser.add_argument("--mr_spectral",action="store_true",
                        help="Use spectral refractive index for water (default=False).")       

    parser.add_argument("--pseudospherical",action="store_true",
                        help="Use pseudospherical approximation, instead of plane parallel (default=False).")                                
    parser.add_argument("--do_single_xtrack",action="store_true",
                        help="for testing - mask all but 1 xtrack (default=False).")

    parser.add_argument("--do_sleave_adjust",action="store_true",
                        help="Do transmittance adjustment of surface leaving radiance (default=False).")

    parser.add_argument("--do_sleave_iso",action="store_true",
                        help="Do isotropic water leaving radiance (default=False).")

    parser.add_argument("--do_diagnostic_only",action="store_true",
                        help="Do diagnostics only, no radiance calculations (default=False).")

    parser.add_argument("--no_add",action="store_true",
                        help="Do NOT write additional output in VLIDORT call (default=False).")  

    parser.add_argument("--no_write_aer",action="store_true",
                        help="Do NOT aerosol output in VLIDORT call (default=False).")  

    parser.add_argument("--no_write_cld",action="store_true",
                        help="Do NOT cloud output in VLIDORT call (default=False).")                                                    
    parser.add_argument("-n", "--nodemax", default=nodemax,
                        help="Max number of nodes to use. "\
                      "(default=%s)"\
                      %nodemax )       


    args = parser.parse_args()
    args.write_add = True
    if args.no_add:
        args.write_add = False
        
    args.write_cld = True
    if args.no_write_cld:
        args.write_cld = False

    args.write_aer = True
    if args.no_write_aer:
        args.write_aer = False     

    args.cloudfree = True
    if args.withcloud:
        args.cloudfree = False

    args.aerosolfree = True
    if args.withaerosol:
        args.aerosolfree = False

    args.gasfree = True
    if args.withgas:
        args.gasfree = False

    args.planeparallel = True
    if args.pseudospherical:
        args.planeparallel = False

    workspace = WORKSPACE(args)

    if not args.dryrun:
        # Submit and monitor jobs
        workspace.handle_jobs()

