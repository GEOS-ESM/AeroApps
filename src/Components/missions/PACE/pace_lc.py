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

mr =  [1.396,1.362,1.349,1.345,1.339,1.335,1.334,1.333,1.332,1.331,1.329,1.326,
      1.323,1.318,1.312,1.306,1.292,1.261]

mr_ch = [200,250,300,337,400,488,515,550,633,694,860,1060,1300,1536,1800,2000,2250,2500]

mr = np.array(mr)
mr_ch = np.array(mr_ch)



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

                # if the job is finished clean up the workspace
                if finished:
                    #print 'Job finished, cleaning up', s, i 
                    finishedJobs = np.append(finishedJobs,ii)
                    errcheck = self.check_for_errors(i,s)               
                    if (errcheck is False):
                        self.destroy_workspace(i,s)
                        if self.nodemax is not None:
                            self.combine_files(i)
                    else:
                        print 'Jobid ',s,' in ',self.dirstring[i],' exited with errors'

            # finished checking up on all the jobs
            # Remove finished jobs from the currently working list
            if len(finishedJobs) != 0:
                print 'deleting finishedJobs',finishedJobs,jobid[workingJobs[finishedJobs]]
                if self.nodexmax is not None:
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
        self.outfilelist = []
        self.addfilelist = []

        if args.channels is None:
            self.get_channels(self.Date)
        else:
            if type(args.channels) is float:
                self.channels = [args.channels]
            else:
                self.channels = args.channels

        if args.nodemax > 1 : 
            self.nodemax = int(args.nodemax)
        else:
            self.nodemax = None

        for ch in self.channels:
            outpath = '{}/{}.{}'.format(args.tmp,self.Date.isoformat(),ch)

            # Copy over some files to working temp directory
            self.create_workdir(outpath,ch)

            # # Create pcf file
            self.write_rc(outpath,ch)

            self.dirstring.append(outpath)




    def create_workdir(self,outpath,ch):
        # create directory     
        if os.path.exists(outpath):
            shutil.rmtree(outpath)
        os.makedirs(outpath)

        # Copy over slurm script    
        outfile = '{}/{}'.format(outpath,self.runfile)
        shutil.copyfile(self.runfile,outfile)

        if self.nodemax is not None:
            self.edit_runFile(outfile)

        # Copy over Aod_EOS.rc
        outfile = '{}/{}'.format(outpath,'Aod_EOS.rc')
        shutil.copyfile('Aod_EOS.rc',outfile)
        self.edit_AODrc(outfile,ch)

        # link over some scripts and files
        source = ['leo_vlidort_cloud.x','ExtData','Chem_MieRegistry.rc',
                  'ExtDataCloud','clean_mem.sh']
        for src in source:
            os.symlink('{}/{}'.format(self.cwd,src),'{}/{}'.format(outpath,src))

    def edit_runFile(self,filename):
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

    def edit_AODrc(self,filename,ch):
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



    def write_rc(self,outpath,ch):
        YMDdir    = self.Date.strftime('Y%Y/M%m/D%d')
        LbDir     = '{}/LevelB/{}'.format(self.rootdir,YMDdir)
        LcDir     = '{}/LevelC/{}'.format(self.rootdir,YMDdir)
        surfDir   = '{}/LevelB/surface'.format(self.rootdir)

        pYMDdir   = self.Date.strftime('Y2020/M%m/D%d')
        L1bDir    = '{}/L1B/{}'.format(self.rootdir,pYMDdir)


        filename = '{}/leo_vlidort_cloud.rc'.format(outpath)
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

        newline  = 'SCALAR: false\n'
        text.append(newline)

        if self.nodemax is not None:
            newline = 'NODEMAX: {}\n'.format(self.nodemax)
        else:
            newline = 'NODEMAX: 1'

        text.append(newline)

        newline = 'ADDITIONAL_OUTPUT: true\n'
        text.append(newline)
        newline = 'VERSION: 1.0\n'
        text.append(newline)
        newline = 'LAYOUT: None\n'
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
            LER_file = '{}/pace-g5nr.lb.omi_ler.{}_{}.nc4'.format(LbDir,nymd,hms)
            newline = 'LER_file: {}\n'.format(LER_file)
            text.append(newline)
        # HYBRID
        elif (ch > 388.) and (ch < 470):
            newline = 'LANDNAME: LER-MCD43C\n'
            text.append(newline)
            newline = 'LANDMODEL: RTLS-HYBRID\n'
            text.append(newline)
            newline = 'LANDBAND_C_LER: 354 388\n'
            text.append(newline)
            LER_file = '{}/pace-g5nr.lb.omi_ler-discover.{}_{}.nc4'.format(LbDir,nymd,hms)
            newline = 'LER_file: {}\n'.format(LER_file)
            newline = 'LANDBAND_C_BRDF: 470 550 650 850 1200 1600 2100\n'
            text.append(newline)
            BRDF_file = '{}/BRDF/MCD43C1/006/{}/pace-g5nr.lb.brdf.{}_{}.nc4'.format(surfDir,YMDdir,nymd,hms)
            newline = 'BRDF_file: {}\n'.format(BRDF_file)
            text.append(newline)
        else:
            newline = 'LANDNAME: MCD43C-BPDF\n'
            text.append(newline)
            newline = 'LANDMODEL: RTLS\n'
            text.append(newline)
            newline = 'LANDBAND_C_BRDF: 470 550 650 850 1200 1600 2100\n'
            text.append(newline)
            BRDF_file = '{}/BRDF/MCD43C1/006/{}/pace-g5nr.lb.brdf.{}_{}.nc4'.format(surfDir,YMDdir,nymd,hms)
            newline = 'BRDF_file: {}\n'.format(BRDF_file)
            BPDF_file = '{}/BPDF/LAND_COVER/MCD12C1/051/{}/pace-g5nr.lb.land_cover.{}_{}.nc4'.format(surfDir,YMDdir,nymd,hms)
            newline = 'BPDF_file: {}\n'.format(BPDF_file)
            NDVI_file = '{}/{}/pace-g5nr.lb-nearest.myd13c2.{}_{}.nc4'.format(LbDir,YMDdir,nymd,hms)
            newline = 'NDVI_file: {}\n'.format(NDVI_file)
            text.append(newline)

        text.append('\n')
        # WATER STUFF
        newline = 'WATERNAME: NOBM-CX\n'
        text.append(newline)
        newline = 'WATERMODEL: CX\n'
        text.append(newline)
        WAT_file = '{}/SLEAVE/NOBM/{}/pace-g5nr.lb.sleave.{}_{}.nc4'.format(surfDir,YMDdir,nymd,hms)
        newline = 'WAT_file: {}\n'.format(WAT_file)
        text.append(newline)

        #refractive index
        mruse = np.interp(float(ch),mr_ch,mr)
        newline = 'WATERMR NOBM-CX\n'
        text.append(newline)
        newline = 'WATERMR: {}\n'.format(mruse)
        text.append(newline)

        text.append('\n')
        # CLOUD STUFF
        newline = 'ICLDTABLE:ExtDataCloud/IceLegendreCoeffs.nc4\n'
        text.append(newline)
        newline = 'LCLDTABLE:ExtDataCloud/WaterLegendreCoeffs.nc4\n'
        text.append(newline)

        cldFile = 'ExtDataCloud/IceLegendreCoeffs.nc4'
        nc = Dataset(cldFile)
        mch = nc.variables['MODIS_wavelength'][:]
        nc.close()
        idx = 1+ np.argmin(np.abs(1e3*mch - ch))
        newline = 'IDXCLD: {}\n'.format(idx)
        text.append(newline)

        text.append('\n')
        # OTHER FILENAMES
        AER_file = '{}/pace-g5nr.lb.aer_Nv.{}_{}.nc4'.format(LbDir,nymd,hms)
        newline = 'AER_file: {}\n'.format(AER_file)
        text.append(newline)

        MET_file = '{}/pace-g5nr.lb.met_Nv.{}_{}.nc4'.format(LbDir,nymd,hms)
        newline = 'MET_file: {}\n'.format(MET_file)
        text.append(newline)

        pDate = isoparser(self.Date.strftime('2020-%m-%dT%H:%M:00'))
        nyj = pDate.strftime('%Y%j')
        ANG_file = '{}/OCI{}{}.L1B_PACE.nc'.format(L1bDir,nyj,hms)
        newline = 'ANG_file: {}\n'.format(ANG_file)
        text.append(newline)

        INV_file = '{}/pace-g5nr.lb.asm_Nx.{}_{}.nc4'.format(LbDir,nymd,hms)
        newline = 'INV_file: {}\n'.format(INV_file)
        text.append(newline)

        OUT_file = '{}/pace-g5nr.lc.vlidort.{}_{}.{}.nc4'.format(LcDir,nymd,hms,int(ch))
        newline = 'OUT_file: {}\n'.format(OUT_file)
        text.append(newline)
        self.outfilelist.append(OUT_file)

        ADD_file = '{}/pace-g5nr.lc.add.{}_{}.nc4'.format(LcDir,nymd,hms)
        newline = 'ADD_file: {}\n'.format(ADD_file)
        text.append(newline)
        self.addfilelist.append(ADD_file)

        CLD_file = '{}/pace-g5nr.TOTWPDF-GCOP-SKEWT.{}_{}.nc4'.format(LcDir,nymd,hms)
        newline = 'CLD_file: {}\n'.format(CLD_file)
        text.append(newline)

        # write text to file
        for l in text:
            f.write(l)

        f.close()

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
        source = ['Aod_EOS.rc','ExtData','ExtDataCloud','leo_vlidort_cloud.rc',
                  'leo_vlidort_cloud.x','Chem_MieRegistry.rc','clean_mem.sh']
        for src in source:
            os.remove(src)

        os.chdir(self.cwd)
        if self.profile is False:
            os.rmdir(self.dirstring[i])

    def combine_files(self,i):
        date_ch = self.dirstring[i].split('/')[-1]
        date,ch = date_ch.split('.')

        outfile = self.outfilelist[i]
        addfile = self.addfilelist[i]

        filelist = []
        for a in np.arange(self.nodemax):
            a = a + 1
            filelist.append(outfile + '_' + str(a) )

        self.do_merge(outfile,filelist)

        filelist = []
        for a in np.arange(self.nodemax):
            a = a + 1
            filelist.append(addfile + '_' + str(a) )
        
        self.do_merge(addfile,filelist)

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

    parser.add_argument("-c","--channels", default=None,type=lambda s: map(int, s.split(",")),
                        help="channels to get TOA radiance (default=None - read in from PACE L1B File)")                            

    parser.add_argument("-e","--extch", default=None,type=lambda s: map(int, s.split(",")),
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
    workspace.handle_jobs()
