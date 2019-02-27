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
from   MAPL.ShaveMantissa_ import shave32
from MAPL.constants import MAPL_UNDEF


mr =  [1.396,1.362,1.349,1.345,1.339,1.335,1.334,1.333,1.332,1.331,1.329,1.326,
      1.323,1.318,1.312,1.306,1.292,1.261]

mr_ch = [200,250,300,337,400,488,515,550,633,694,860,1060,1300,1536,1800,2000,2250,2500]

mr = np.array(mr)
mr_ch = np.array(mr_ch)

nlev = 72

SDS_RT_blue = {'surf_ref_I_blue': ['surface reflectance for blue CCD','None',None]}
SDS_RT_red  = {'surf_ref_I_red':  ['surface reflectance for red CCD','None',None]}
SDS_RT_SWIR = {'surf_ref_I_SWIR': ['surface reflectance for SWIR bands','None',None]}

SDS_ADD_blue = {'rot_blue': ['rayleigh optical thickness for blue CCD','None',nlev]}
SDS_ADD_red  = {'rot_red':  ['rayleigh optical thickness for red CCD','None',nlev]}
SDS_ADD_SWIR = {'rot_SWIR': ['rayleigh optical thickness for SWIR bands','None',nlev]}
          
          # 'pe': ['layer edge pressure','Pa',nlev+1],
          # 'ze': ['layer edge altitude','m',nlev+1],
          # 'te': ['layer edge temperature','K',nlev+1]} 

SDS_CLD_blue = {'lcot_blue': ['liquid cloud optical depth for blue CCD','None',nlev],
                'lc_ssa_blue':  ['liquid cloud single scattering albedo for blue CCD','None',nlev],
                # 'lc_g_blue':  ['liquid cloud assymetry parameter for blue CCD','None',nlev],
                'icot_blue': ['ice cloud optical depth for blue CCD','None',nlev],
                'ic_ssa_blue':  ['ice  cloud single scattering albedo for blue CCD','None',nlev],
                # 'ic_g_blue':  ['ice  cloud assymetry parameter for blue CCD','None',nlev],
                }

SDS_CLD_red  = {'lcot_red': ['liquid cloud optical depth for red CCD','None',nlev],
                'lc_ssa_red':  ['liquid cloud single scattering albedo for red CCD','None',nlev],
                # 'lc_g_red':  ['liquid cloud assymetry parameter for red CCD','None',nlev],    
                'icot_red': ['ice  cloud optical depth for red CCD','None',nlev],
                'ic_ssa_red':  ['ice  cloud single scattering albedo for red CCD','None',nlev],
                # 'ic_g_red':  ['ice  cloud assymetry parameter for red CCD','None',nlev],    
                }
SDS_CLD_SWIR = {'lcot_SWIR': ['liquid cloud optical depth for SWIR bands','None',nlev],
                'lc_ssa_SWIR':  ['liquid cloud single scattering albedo for SWIR bands','None',nlev],
                # 'lc_g_SWIR':  ['liquid cloud assymetry parameter for SWIR bands','None',nlev],
                'icot_SWIR': ['ice  cloud optical depth for SWIR bands','None',nlev],
                'ic_ssa_SWIR':  ['ice  cloud single scattering albedo for SWIR bands','None',nlev],
                # 'ic_g_SWIR':  ['ice  cloud assymetry parameter for SWIR bands','None',nlev],
                }

SDS_RT = {'red': SDS_RT_red,
          'blue': SDS_RT_blue,
          'SWIR': SDS_RT_SWIR}

SDS_ADD = {'red': SDS_ADD_red,
          'blue': SDS_ADD_blue,
          'SWIR': SDS_ADD_SWIR}          

SDS_CLD = {'red': SDS_CLD_red,
          'blue': SDS_CLD_blue,
          'SWIR': SDS_CLD_SWIR}                

#---
def shave(q,undef=MAPL_UNDEF,has_undef=1,nbits=12):
    """
    Shave variable. On input, nbits is the number of mantissa bits to keep
    out of maximum of 24.
    """

    # Determine shaving parameters
    # ----------------------------
    xbits = 24 - nbits
    shp = q.shape
    rank = len(shp)
    if rank == 2:  # yx
        chunksize = shp[0]*shp[1] 
    elif rank == 3: # zyx
        chunksize = shp[1]*shp[2]
    else:
        raise ValueError, "invalid rank=%d"%rank

    # Shave it
    # --------
    qs, rc = shave32(q.ravel(),xbits,has_undef,undef,chunksize)
    if rc:
        raise ValueError, "error on return from shave32, rc=%d"%rc

    return qs.reshape(shp)

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
                        self.errTally[i] = False
                        self.destroy_workspace(i,s)
                        if self.nodemax is not None:
                            self.combine_files(i)
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


    def check_for_errors(self,i,jobid):
        os.chdir(self.dirstring[i])  

        error = False 
        if self.nodemax is not None:
            for a in np.arange(nodemax):
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

        self.cwd = os.getcwd()
        self.runfile   = args.slurm
        self.profile   = args.profile
        self.rootdir   = args.rootdir
        self.verbose   = args.verbose
        self.bpdf_name = args.bpdf_name
        self.rtls_name = args.rtls_name
        self.write_add = args.write_add
        self.write_aer = args.write_aer
        self.write_cld = args.write_cld

        self.dirstring = []
        self.outfilelist = []
        self.addfilelist = []
        self.aerfilelist = []
        self.cldfilelist = []

        if args.channels is None:
            self.get_channels(self.Date)
        else:
            if ',' in args.channels:
                makelist=lambda s: map(int, s.split(","))
                self.channels = makelist(args.channels)
            elif ':' in args.channels:
                makelist=lambda s: map(int, s.split(":"))
                start,stop,delta = makelist(args.channels)
                self.channels = range(start,stop+delta,delta)
            else:
                self.channels = [int(args.channels)]

        if int(args.nodemax) > 1 : 
            self.nodemax = int(args.nodemax)
        else:
            self.nodemax = None

        for ch in self.channels:
            fch = "{:.2f}".format(ch).replace('.','d')
            outpath = '{}/{}.{}'.format(args.tmp,self.Date.isoformat(),fch)

            # Copy over some files to working temp directory
            self.create_workdir(outpath,ch)

            # # Create pcf file
            self.write_rc(outpath,ch)

            self.dirstring.append(outpath)


        self.errTally    = np.ones(len(self.dirstring)).astype(bool)

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



    def write_rc(self,outpath,ch):
        YMDdir    = self.Date.strftime('Y%Y/M%m/D%d')
        LbDir     = '{}/LevelB/{}'.format(self.rootdir,YMDdir)
        LcDir     = '{}/LevelC/{}'.format(self.rootdir,YMDdir)
        LcDirCh   = '{}/LevelC/{}/channel'.format(self.rootdir,YMDdir)
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
            newline = 'NODEMAX: 1\n'

        text.append(newline)

        if self.write_add:
            newline = 'ADDITIONAL_OUTPUT: true\n'
        else:
            newline = 'ADDITIONAL_OUTPUT: false\n'
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
        if ch <= 775.0:
            newline = 'WATERNAME: NOBM-CX\n'
            text.append(newline)
            WAT_file = '{}/SLEAVE/NOBM/{}/pace-g5nr.lb.sleave.{}_{}.nc4'.format(surfDir,YMDdir,nymd,hms)
            newline = 'WAT_file: {}\n'.format(WAT_file)
            text.append(newline)
        else:
            newline = 'WATERNAME: CX\n'
            text.append(newline)

        newline = 'WATERMODEL: CX\n'
        text.append(newline)

        #refractive index
        mruse = np.interp(float(ch),mr_ch,mr)
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

        fch = "{:.2f}".format(ch).replace('.','d')
        OUT_file = '{}/pace-g5nr.lc.vlidort.{}_{}.{}.nc4'.format(LcDirCh,nymd,hms,fch)
        newline = 'OUT_file: {}\n'.format(OUT_file)
        text.append(newline)
        self.outfilelist.append(OUT_file)

        if self.write_add:
            ADD_file = '{}/pace-g5nr.lc.add.{}_{}.{}.nc4'.format(LcDirCh,nymd,hms,fch)
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
        self.channels  = list(np.unique(np.concatenate(channels)))


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



def populate_L1B(outfilelist,rootdir,channels,Date,force=False):
        YMDdir    = Date.strftime('Y%Y/M%m/D%d')
        pYMDdir   = Date.strftime('Y2020/M%m/D%d')
        L1bDir    = '{}/L1B/{}'.format(rootdir,pYMDdir)
        Lc2Dir    = '{}/LevelC2/{}'.format(rootdir,YMDdir)
        hms       = Date.strftime('%H%M00')

        pDate = isoparser(Date.strftime('2020-%m-%dT%H:%M:00'))
        nyj = pDate.strftime('%Y%j')
        L1B_file = '{}/OCI{}{}.L1B_PACE.nc'.format(L1bDir,nyj,hms)

        outfile = '{}/OCI{}{}.L1B_PACE.nc'.format(Lc2Dir,nyj,hms)
        exists  = os.path.isfile(outfile)
        if force or (not exists):
            shutil.copyfile(L1B_file,outfile)


        ncmerge = Dataset(outfile,mode='r+')
        SDS = {'blue_wavelength' : 'Lt_blue',
               'red_wavelength'  : 'Lt_red',
               'SWIR_wavelength' : 'Lt_SWIR'}
        for sds in SDS:
            pchannels = ncmerge.groups['sensor_band_parameters'].variables[sds][:]
            pvar      = ncmerge.groups['observation_data'].variables[SDS[sds]]

            for ch,filename in zip(channels,outfilelist):
                for i,pch in enumerate(pchannels):
                    if float(ch) == pch:
                        print 'inserting ',filename
                        nc = Dataset(filename)
                        fch = "{:.2f}".format(ch)
                        data = np.squeeze(nc.variables['I_'+fch][:])

                        ii = data == -500.0
                        if np.sum(ii) > 0:
                            if (not hasattr(data,'mask')):
                                data = np.ma.array(data)
                                data.mask = ii
                            elif (type(data.mask) is np.bool_):
                                data.mask = ii
                            else:
                                data.mask[ii] = True
                        pvar[i,:,:] = data

                        nc.close()

        ncmerge.close()

def condense_LC(outfilelist,addfilelist,cldfilelist,aerfilelist,rootdir,channels,Date,force=False):
        YMDdir    = Date.strftime('Y%Y/M%m/D%d')
        pYMDdir   = Date.strftime('Y2020/M%m/D%d')
        LcDir     = '{}/LevelC/{}'.format(rootdir,YMDdir)
        L1bDir    = '{}/L1B/{}'.format(rootdir,pYMDdir)
        hms       = Date.strftime('%H%M00')
        nymd      = Date.strftime('%Y%m%d')

        pDate = isoparser(Date.strftime('2020-%m-%dT%H:%M:00'))
        nyj = pDate.strftime('%Y%j')
        L1B_file = '{}/OCI{}{}.L1B_PACE.nc'.format(L1bDir,nyj,hms)

        # Create file if you need to
        for wav in ['blue','red','SWIR']:
            outfile = '{}/pace-g5nr.lc.vlidort_{}.{}_{}.nc4'.format(LcDir,wav,nymd,hms)
            exists  = os.path.isfile(outfile)
            if force or (not exists):
                # create new outfile
                if addfilelist is None:
                    SDS = SDS_RT[wav]
                else:
                    SDS = dict(SDS_RT[wav], **SDS_ADD[wav])
                create_condenseFile(L1B_file,outfile,Date,SDS)

            # Condense RT stuff
            SDS = SDS_RT[wav]        
            insert_condenseVar(outfile,SDS,channels,outfilelist)

            # Condense ADD stuff
            if addfilelist is not None:
                SDS = SDS_ADD[wav]
                insert_condenseVar(outfile,SDS,channels,addfilelist)

            if cldfilelist is not None:
                # Create file if you need to
                outfile = '{}/pace-g5nr.cloud_{wav}.{}_{}.nc4'.format(LcDir,wav,nymd,hms)
                exists  = os.path.isfile(outfile)
                if force or (not exists):
                    # create new outfile
                    SDS = SDS_CLD[wav]
                    create_condenseFile(L1B_file,outfile,Date,SDS)

                # Condense Cloud stuff
                SDS = SDS_CLD[wav]        
                insert_condenseVar(outfile,SDS,channels,cldfilelist)



def insert_condenseVar(outfile,SDS,channels,outfilelist):
    # insert data into correct place in outfile
    
    for sds in SDS:

        if 'blue' in sds:
            chname = 'blue_wavelength'
        elif 'red' in sds:
            chname = 'red_wavelength'
        elif 'SWIR' in sds:
            chname = 'SWIR_wavelength'
        else:
            chname = None

        if chname is not None:
            ncmerge = Dataset(outfile,mode='r')
            pchannels = ncmerge.variables[chname][:]
            ncmerge.close()
        else:
            pchannels = None

        if '_' in sds:
            oname = '_'.join(sds.split('_')[:-1])
        else:
            oname = sds


        if pchannels is not None:
            for ch,filename in zip(channels,outfilelist):
                for i,pch in enumerate(pchannels):
                    if float(ch) == pch:
                        ncmerge = Dataset(outfile,mode='r+')
                        pvar      = ncmerge.variables[sds]
                        print 'inserting ',oname, filename, ' to ',sds
                        nc = Dataset(filename)
                        fch = "{:.2f}".format(ch)
                        ovar = nc.variables[oname + '_'+fch]
                        missing_value = ovar._FillValue
                        data = np.squeeze(ovar[:])
                        rank = len(data.shape)
                        if rank == 1:
                            pvar[:,i] = data
                        if rank == 2:
                            pvar[:,:,i] = shave(data,undef=missing_value)
                        if rank == 3:
                            pvar[:,:,:,i] = shave(data,undef=missing_value)

                        nc.close()
                        ncmerge.close()


        if pchannels is None:
            ncmerge = Dataset(outfile,mode='r+')
            pvar      = ncmerge.variables[sds]
            filename = outfilelist[0]
            print 'inserting ',oname,filename, ' to ',sds
            nc = Dataset(filename)
            ovar = nc.variables[oname]
            data = np.squeeze(ovar[:])
            missing_value = ovar._FillValue
            rank = len(data.shape)
            if rank == 1:
                pvar[:] = data
            else:
                pvar[:] = shave(data,undef=missing_value)

            nc.close()
            ncmerge.close()
                        
    

def create_condenseFile(L1B_file,outfile,Date,SDS):

    # Open NC file
    # ------------
    print 'create outfile',outfile
    nc = Dataset(outfile,'w',format='NETCDF4_CLASSIC')
    print 'copy from',L1B_file
    nctrj = Dataset(L1B_file)

    # Get Dimensions
    # ----------------------
    npixel = len(nctrj.dimensions['ccd_pixels'])
    nscan  = len(nctrj.dimensions['number_of_scans'])
    nblue  = len(nctrj.dimensions['blue_bands'])
    nred   = len(nctrj.dimensions['red_bands'])
    nswir  = len(nctrj.dimensions['SWIR_bands'])


    # Set global attributes
    # ---------------------
    nc.title = "RT inputs for VLIDORT Simulation of GEOS-5 PACE Sampled Data"
    nc.institution = 'NASA/Goddard Space Flight Center'
    nc.source = 'Global Model and Assimilation Office'
    nc.history = 'Condensed outputs from channel files'
    nc.references = 'n/a'
    nc.comment = 'n/a'
    nc.contact = 'Patricia Castellanos <patricia.castellanos@nasa.gov>'
    nc.Conventions = 'CF'    

    # Create dimensions
    # -----------------
    l = nc.createDimension('lev',nlev)
    le = nc.createDimension('leve',nlev+1)
    x = nc.createDimension('ccd_pixels',npixel)
    xs = nc.createDimension('SWIR_pixels',npixel)
    y = nc.createDimension('number_of_scans',nscan)
    b = nc.createDimension('blue_bands',nblue)
    r = nc.createDimension('red_bands',nred)
    s = nc.createDimension('SWIR_bands',nswir)

    # Save lon/lat
    # --------------------------
    _copyVar(nctrj,nc,u'longitude','geolocation_data',dtype='f4',zlib=True,verbose=True)
    _copyVar(nctrj,nc,u'latitude','geolocation_data',dtype='f4',zlib=True,verbose=True)
    _copyVar(nctrj,nc,u'ev_mid_time','scan_line_attributes', dtype='f4',zlib=True,verbose=True)
    _copyVar(nctrj,nc,u'blue_wavelength','sensor_band_parameters', dtype='f4',zlib=True,verbose=True)
    _copyVar(nctrj,nc,u'red_wavelength','sensor_band_parameters', dtype='f4',zlib=True,verbose=True)
    _copyVar(nctrj,nc,u'SWIR_wavelength','sensor_band_parameters', dtype='f4',zlib=True,verbose=True)


    # Loop over SDS creating each dataset
    #---------------------------------------
    for sds in SDS:
        lname,unit,levs = SDS[sds]
        if levs == nlev:
            if 'red' in sds:
                dim = ('lev','number_of_scans','ccd_pixels','red_bands',)
            elif 'blue' in sds:
                dim = ('lev','number_of_scans','ccd_pixels','blue_bands',)
            elif 'SWIR' in sds:
                dim = ('lev','number_of_scans','ccd_pixels','SWIR_bands',)
            else:
                dim = ('lev','number_of_scans','ccd_pixels')

        elif levs == nlev+1:
            if 'red' in sds:
                dim = ('leve','number_of_scans','ccd_pixels','red_bands',)
            elif 'blue' in sds:
                dim = ('leve','number_of_scans','ccd_pixels','blue_bands',)
            elif 'SWIR' in sds:
                dim = ('leve','number_of_scans','ccd_pixels','SWIR_bands',)
            else:
                dim = ('leve','number_of_scans','ccd_pixels')

        else:
            if 'red' in sds:
                dim = ('number_of_scans','ccd_pixels','red_bands',)
            elif 'blue' in sds:
                dim = ('number_of_scans','ccd_pixels','blue_bands',)
            elif 'SWIR' in sds:
                dim = ('number_of_scans','ccd_pixels','SWIR_bands',)
            else:
                dim = ('number_of_scans','ccd_pixels')


        this = nc.createVariable(sds,'f4',dim,zlib=True,fill_value=-999.0)  

        this.long_name = lname
        this.missing_value = -999.0
        this.unit = unit  


    nc.close()
    nctrj.close()          


#----
def _copyVar(ncIn,ncOut,name,group,dtype='f4',zlib=False,verbose=False):
    """
    Create variable *name* in output file and copy its
    content over,
    """
    x = ncIn.groups[group].variables[name]
    if verbose:
        print 'copy variable ',name,x.dimensions
    y = ncOut.createVariable(name,dtype,x.dimensions,zlib=zlib)
    if hasattr(x,'long_name'): y.long_name = x.long_name
    if hasattr(x,'units'): y.units = x.units
    if hasattr(x,'missing_value'): y.missing_value = x.missing_value
    rank = len(x.shape)

    if rank == 1:
        y[:] = x[:]
    elif rank == 2:
        if hasattr(x,'missing_value'):
            y[:,:] = shave(x[:,:],undef=x.missing_value)
        else:
            y[:,:] = shave(x[:,:],has_undef=0)
    elif rank == 3:
        if hasattr(x,'missing_value'):
            y[:,:,:] = shave(x[:,:,:],undef=x.missing_value)
        else:
            y[:,:,:] = shave(x[:,:,:],has_undef=0)
    else:
        raise ValueError, "invalid rank of <%s>: %d"%(name,rank)


if __name__ == '__main__':
    
    #Defaults
    nodemax   = 20
    slurm     = 'pace_lc_array.j'
    rootdir   = '/discover/nobackup/projects/gmao/osse2/pub/c1440_NR/OBS/PACE'
    tmp       = '/discover/nobackup/projects/gmao/osse2/pub/c1440_NR/OBS/PACE/workdir'
    rtls_name = 'MCD43C'
    bpdf_name = None # 'MAIGNAN'

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

    parser.add_argument("-c","--channels", default=None,
                        help="channels to get TOA radiance. Can be a list (1,2) or a range (1:2:1) (default=None - read in from PACE L1B File)")                            

    parser.add_argument("-e","--extch", default=None,type=lambda s: map(int, s.split(",")),
                        help="channels to run extinction sampler (default=None - read in from PACE L1B File)")  

    parser.add_argument('--bpdf_name',default=bpdf_name,
                        help="BPDF model name (default=%s)"%bpdf_name)                          

    parser.add_argument('--rtls_name',default=rtls_name,
                        help="RTLS data name (default=%s)"%rtls_name)                          

    parser.add_argument("--norad",action="store_true",
                        help="No radiance calculations (default=False).")

    parser.add_argument("--no_add",action="store_true",
                        help="Do NOT write additional output in VLIDORT call (default=False).")  

    parser.add_argument("--write_aer",action="store_true",
                        help="Do aerosol output in VLIDORT call (default=False).")  

    parser.add_argument("--write_cld",action="store_true",
                        help="Do cloud output in VLIDORT call (default=False).")                                                    

    parser.add_argument("--doext",action="store_true",
                        help="Do extinctions calculations (default=False).")

    parser.add_argument("--force",action="store_true",
                        help="Overwrite existing L1B file in LevelC directory (default=False).")    

    parser.add_argument("-n", "--nodemax", default=nodemax,
                        help="Max number of nodes to use. "\
                      "(default=%s)"\
                      %nodemax )       


    args = parser.parse_args()
    if args.no_add:
        args.write_add = False
    else:
        args.write_add = True


    workspace = WORKSPACE(args)

    if not args.dryrun:
        # Submit and monitor jobs
        if not args.norad:
            workspace.handle_jobs()

            # Take VLIDORT outputs and populate PACE L1b File
            I = ~workspace.errTally
            if any(I):
                outfilelist = np.array(workspace.outfilelist)[I]                
                channels    = np.array(workspace.channels)[I]
                rootdir     = workspace.rootdir
                Date        = workspace.Date
                if workspace.write_add:
                    addfilelist = np.array(workspace.addfilelist)[I]
                else:
                    addfilelist = None
                if workspace.write_aer:
                    aerfilelist = np.array(workspace.aerfilelist)[I]
                else:
                    aerfilelist = None
                if workspace.write_cld:
                    cldfilelist = np.array(workspace.cldfilelist)[I]
                else:
                    cldfilelist = None                    
                populate_L1B(outfilelist,rootdir,channels,Date,force=args.force)
                condense_LC(outfilelist,addfilelist,cldfilelist,aerfilelist,rootdir,channels,Date,force=args.force)


        if args.doext:
            # calls extinction sampler
            pass