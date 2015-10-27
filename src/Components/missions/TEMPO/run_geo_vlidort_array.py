#!/bin/python
# -*- coding: utf-8 -*-
""" Runscript for geo_vlidort episodes"""
from datetime import datetime, timedelta 
from dateutil.parser import parse
import os
import subprocess
from distutils.dir_util import mkpath
import numpy as np
import math
import time
import glob
import shutil
from netCDF4 import Dataset

def make_workspace(date,ch,code,outdir,runfile,prefix='workdir',nodemax=None,addoutdir=None):
    cwd     = os.getcwd()
    dirname = prefix + '/'+str(date.date())+'T'+str(date.hour).zfill(2)+'.'+ch+'.'+code
    jobname = str(date.date())+'T'+str(date.hour).zfill(2)+'.'+ch
    outdir = outdir + '/Y'+str(date.year)+'/M'+str(date.month).zfill(2)+'/D'+str(date.day).zfill(2)
    if addoutdir is not None:
        addoutdir = addoutdir + '/Y'+str(date.year)+'/M'+str(date.month).zfill(2)+'/D'+str(date.day).zfill(2)
    bindir = os.getcwd() + '/' + dirname
    if not os.path.exists(dirname):
        os.makedirs(dirname)
    
    os.chdir(dirname)
    if not os.path.isfile('geo_vlidort.x'):
        os.symlink(cwd+'/geo_vlidort.x','geo_vlidort.x')
    if not os.path.exists('ExtData'):
        os.symlink(cwd+'/ExtData','ExtData')
    if not os.path.isfile('Chem_MieRegistry.rc'):
        os.symlink(cwd+'/Chem_MieRegistry.rc','Chem_MieRegistry.rc')
    if not os.path.isfile('clean_mem.sh'):
        os.symlink(cwd+'/clean_mem.sh','clean_mem.sh')
        
    source = open(cwd+'/'+runfile,'r')
    destination = open(runfile,'w')
    for line in source:
        if (line[0:18] == '#SBATCH --job-name'):            
            destination.write('#SBATCH --job-name='+jobname+'\n')

        elif (line[0:16] == '#SBATCH --output'):
            if (nodemax is not None and nodemax > 1):
                destination.write('#SBATCH --output='+'slurm'+"_%A_%a.out"+'\n')
            else:
                destination.write('#SBATCH --output='+'slurm'+"_%j.out"+'\n')

        elif (line[0:15] == '#SBATCH --error'):
            if (nodemax is not None and nodemax > 1):
                destination.write('#SBATCH --error='+'slurm'+"_%A_%a.err"+'\n')      
            else:
                destination.write('#SBATCH --error='+'slurm'+"_%j.err"+'\n')

        elif (line[0:15] == '#SBATCH --array'):
            if (nodemax is not None and nodemax > 1):
                destination.write('#SBATCH --array=1-'+str(nodemax)+'\n') 

        elif (line[0:15] == 'setenv TEMPOBIN'):
            destination.write('setenv TEMPOBIN '+bindir+'\n')

        elif (line[0:13] == 'setenv OUTDIR'):            
            destination.write('setenv OUTDIR '+outdir+'\n')            
            if not os.path.exists(outdir):
                mkpath(outdir)

        elif (line[0:16] == 'setenv ADDOUTDIR'):   
            if (addoutdir is not None):       
                destination.write('setenv ADDOUTDIR '+addoutdir+'\n')                          
                if not os.path.exists(addoutdir):
                    mkpath(addoutdir)
            else:
                destination.write('setenv ADDOUTDIR None \n') 

        elif (line[0:8] == '$RUN_CMD'):
            if (nodemax is not None and nodemax > 1):     
                destination.write('$RUN_CMD  ./geo_vlidort.x geo_vlidort.rc ${SLURM_ARRAY_TASK_ID}'+'\n')    
            else:
                destination.write('$RUN_CMD  ./geo_vlidort.x geo_vlidort.rc'+'\n')   

        else:
            destination.write(line)        
    source.close()
    destination.close()
    
    source = open(cwd+'/Aod_EOS.rc','r')
    destination = open('Aod_EOS.rc','w')
    a = float(ch)*1e-3
    for line in source:
        if (line[0:11] == 'r_channels:'):
            destination.write('r_channels: '+'{:0.3f}e-6'.format(a)+'\n')
        else:
            destination.write(line) 
    source.close()
    destination.close()
     
    os.chdir(cwd) 
    
    if (addoutdir is not None):
        return dirname, outdir, addoutdir
    else:
        return dirname, outdir 

def check_for_errors(dirname,jobid,nodemax=None):
    cwd     = os.getcwd()
    os.chdir(dirname)  

    error = False  
    if nodemax is not None and nodemax > 1:
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

    os.chdir(cwd)
    return error


def destroy_workspace(jobid,dirname,outdir,addoutdir=None, nodemax=None):
    cwd     = os.getcwd()
    os.chdir(dirname)

    os.remove('Aod_EOS.rc')
    os.remove('Chem_MieRegistry.rc')
    os.remove('geo_vlidort.x')
    os.remove('clean_mem.sh')
    os.remove('ExtData')
    os.remove('geo_vlidort.rc')
    os.remove('geo_vlidort_run_array.j')

    if nodemax is not None and nodemax > 1:
        for a in np.arange(nodemax):
            a = a + 1
            errfile = 'slurm_' +jobid + '_' + str(a) + '.err'
            os.remove(errfile)
            outfile = 'slurm_' +jobid + '_' + str(a) + '.out'
            os.remove(outfile)
        os.remove('slurm_%A_%a.out')

    else:
        errfile = 'slurm_' +jobid + '.err'
        os.remove(errfile)        
        outfile = 'slurm_' +jobid + '.out'
        os.remove(outfile)        

    def move_file(filelist,dest):
        for movefile in filelist:
            if os.path.exists(dest+'/'+movefile):
                os.remove(dest+'/'+movefile)
            shutil.move(movefile,dest)

    #runscript to combine files
    if nodemax is not None and nodemax > 1:
        outfilelist = glob.glob('*.lc2.*.nc4')
        combine_files(outfilelist)
        outfilelist = glob.glob('*.lc2.*.nc4')
        move_file(outfilelist,outdir)
    else:
        outfilelist = glob.glob('*.lc2.*.nc4')
        move_file(outfilelist,outdir)


    if (addoutdir is not None):
        if nodemax is not None and nodemax > 1:
            outfilelist = glob.glob('*.lc.*.nc4')
            combine_files(outfilelist)
            outfilelist = glob.glob('*.lc.*.nc4')
            move_file(outfilelist,addoutdir)
        else:
            outfilelist = glob.glob('*.lc.*.nc4')
            move_file(outfilelist,addoutdir)

    os.chdir(cwd)
    os.rmdir(dirname)

def combine_files(filelist):
    mergedfile = filelist[0]
    mergedfile = mergedfile[:-5] + 'nc4'

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


def make_maiac_rcfile(dirname,indir,date,ch,code,interp,additional_output,nodemax=None,i_band=None):
    cwd = os.getcwd()
    os.chdir(dirname)

    rcfile = open('geo_vlidort.rc','w')
    rcfile.write('INDIR: '+indir+'\n')
    rcfile.write('OUTDIR: .\n')
    rcfile.write('DATE: '+str(date.year)+str(date.month).zfill(2)+str(date.day).zfill(2)+'\n')
    rcfile.write('TIME: '+str(date.hour).zfill(2)+'\n')
    rcfile.write('INSTNAME: tempo\n')
    rcfile.write('SURFNAME: MAIACRTLS\n')

    #figure out correct MODIS doy
    if (str(startdate.date()) == '2005-12-31'):
        rcfile.write('SURFDATE: 2006008\n')
    elif ((str(startdate.year) == '2007') and (str(date.month).zfill(2) == '03')):
        rcfile.write('SURFDATE: 2007072\n')
    elif ((str(startdate.year) == '2007') and (str(date.month).zfill(2) == '04')):
        rcfile.write('SURFDATE: 2007120\n')
    elif ((str(startdate.year) == '2006') and (str(date.month).zfill(2) == '07')):
        rcfile.write('SURFDATE: 2006216\n')
    elif ((str(startdate.year) == '2006') and (str(date.month).zfill(2) == '08')):
        rcfile.write('SURFDATE: 2006216\n')        
    else:
        doy = date.toordinal() - datetime(date.year-1,12,31).toordinal()
        DOY = 8*(int(doy/8) + 1)
        if (DOY > 365):
            DOY = 8
            rcfile.write('SURFDATE: '+str(date.year+1)+str(DOY).zfill(3)+'\n')
        else:
            rcfile.write('SURFDATE: '+str(date.year)+str(DOY).zfill(3)+'\n')


    rcfile.write('SURFBAND: ' + interp +'\n')
    if (interp.upper() == 'INTERPOLATE'):
        rcfile.write('SURFBAND_C: 645 858 469 555 1240 1640 2130\n')
    else:
        rcfile.write('SURFBAND_I: '+ i_band)

    rcfile.write('SURFBANDM: 7 \n')

    if (code == 'scalar'):
        rcfile.write('SCALAR: true\n')
    else:
        rcfile.write('SCALAR: false\n')


    rcfile.write('CHANNELS: '+ch+'\n')
    rcfile.write('CLDMAX: 0.01\n')
    if nodemax is not None:
        rcfile.write('NODEMAX: '+ str(nodemax) + '\n')

    if additional_output:
        rcfile.write('ADDITIONAL_OUTPUT: true\n')
    else:
        rcfile.write('ADDITIONAL_OUTPUT: false\n')
    rcfile.close()

    os.chdir(cwd)

def make_ler_rcfile(dirname,indir,date,ch,code,interp,additional_output,nodemax=None,i_band=None):
    cwd = os.getcwd()
    os.chdir(dirname)

    rcfile = open('geo_vlidort.rc','w')
    rcfile.write('INDIR: ' + indir + '\n')
    rcfile.write('OUTDIR: .\n')
    rcfile.write('DATE: ' + str(date.year) + str(date.month).zfill(2) + str(date.day).zfill(2) + '\n')
    rcfile.write('TIME: ' + str(date.hour).zfill(2) + '\n')
    rcfile.write('INSTNAME: tempo\n')
    rcfile.write('SURFNAME: LER\n')

    rcfile.write('SURFDATE: ' + str(date.month).zfill(2) +'\n')

    rcfile.write('SURFBAND: ' + interp + '\n')
    if (interp.upper() == 'INTERPOLATE'):
        rcfile.write('SURFBAND_C: 354 388\n')
    else:
        rcfile.write('SURFBAND_I: ' + i_band + '\n')

    rcfile.write('SURFBANDM: 2 \n')

    if (code == 'scalar'):
        rcfile.write('SCALAR: true\n')
    else:
        rcfile.write('SCALAR: false\n')


    rcfile.write('CHANNELS: ' + ch + '\n')
    rcfile.write('CLDMAX: 0.01\n')
    if nodemax is not None:
        rcfile.write('NODEMAX: '+ str(nodemax) + '\n')

    if additional_output:
        rcfile.write('ADDITIONAL_OUTPUT: true\n')
    else:
        rcfile.write('ADDITIONAL_OUTPUT: false\n')
    rcfile.close()

    os.chdir(cwd)    

def prefilter(date,indir):
    g5dir = indir + '/LevelB/'+ 'Y'+ str(date.year) + '/M' + str(date.month).zfill(2) + '/D' + str(date.day).zfill(2) 
    nymd  = str(date.year) + str(date.month).zfill(2) + str(date.day).zfill(2)
    hour  = str(date.hour).zfill(2)
    met   = g5dir + '/tempo-g5nr.lb2.met_Nv.' + nymd + '_' + hour + 'z.nc4'
    geom  = g5dir + '/tempo.lb2.angles.' + nymd + '_' + hour + 'z.nc4'
    land  = indir + '/LevelB/invariant/tempo-g5nr.lb2.asm_Nx.nc4'  

    ncMet = Dataset(met)
    Cld   = np.squeeze(ncMet.variables[u'CLDTOT'][:])
    f     = np.where(Cld > 0.01)
    ncMet.close()
    if len(f[0]) == 0:
        return False, 0

    ncGeom = Dataset(geom)
    SZA    = np.squeeze(ncGeom.variables[u'solar_zenith'][:])
    VZA    = np.squeeze(ncGeom.variables[u'sensor_zenith'][:])
    ncGeom.close()
    ncLand = Dataset(land)
    FRLAND = np.squeeze(ncLand.variables[u'FRLAND'][:])
    ncLand.close()

    SZA = SZA[f]
    VZA = VZA[f]
    FRLAND = FRLAND[f]

    f   = np.where(VZA < 80)
    if len(f[0]) == 0:
        return False, 0

    SZA    = SZA[f]
    FRLAND = FRLAND[f]
    f      = np.where(SZA < 80)
    if len(f[0]) == 0:
        return False, 0

    FRLAND = FRLAND[f]
    f      = np.where(FRLAND >= 0.99)
    if len(f[0]) == 0:
        return False, 0

    return True, len(f[0])

        
#########################################################

if __name__ == "__main__":
    
    startdate         = '2005-12-31T17:00:00'
    enddate           = '2005-12-31T17:00:00'
    channels          = '550'
    surface           = 'MAIACRTLS'
    interp            = 'interpolate'
    i_band            = None    
    additional_output = False
    nodemax           = 6

    runfile           = 'geo_vlidort_run_array.j'
    nccs              = '/discover/nobackup/projects/gmao/osse2/pub/c1440_NR/OBS/TEMPO/DATA/'

    ################
    ###
    #    End of uper inputs
    ###
    ################
    indir             = nccs 
    outdir            = nccs + 'LevelC2'
    if (additional_output):
        addoutdir         = nccs + 'LevelC'
    else:
        addoutdir         = None
    

    dt = timedelta(hours=1)
    startdate = parse(startdate)
    enddate   = parse(enddate)

    # Create working directories, SLURM scripts, and RC-files
    dirstring = np.empty(0)
    outdirstring = np.empty(0)
    if (additional_output):
        addoutdirstring = np.empty(0)

    if type(channels) is str:
        channels = channels.split()

    if type(i_band) is str:
        i_band = i_band.split()

    while (startdate <= enddate):
        anypixels, numpixels = prefilter(startdate,indir) 
        if (anypixels):

            if (numpixels <= 1000 and nodemax is not None):
                nodemax = 1

            for i, ch in enumerate(channels):
                band_i = None
                if (i_band is not None):
                    band_i = i_band[i]

                # Vector Case
                code = 'vector.'
                if (interp.lower() == 'interpolate'):
                    code = code + 'i'                 

                code = code + surface
                
                dirlist = make_workspace(startdate,ch,code,outdir,runfile,nodemax=nodemax,addoutdir=addoutdir)

                if (additional_output):
                    workdir, outdir_, addoutdir_ = dirlist
                else:
                    workdir, outdir_ = dirlist

                if (surface.upper() == 'MAIACRTLS'):
                    make_maiac_rcfile(workdir,indir,startdate,ch,'vector',\
                                      interp,additional_output,nodemax=nodemax,i_band=band_i)
                else:
                    make_ler_rcfile(workdir,indir,startdate,ch,'vector',\
                                    interp,additional_output,nodemax=nodemax,i_band=band_i)
                
                dirstring = np.append(dirstring,workdir)
                outdirstring = np.append(outdirstring,outdir_)
                if (additional_output):
                    addoutdirstring = np.append(addoutdirstring, addoutdir_)

       
        startdate = startdate + dt

    # Submite Jobs      
    runlen  = len(dirstring)   
    if nodemax is not None:
        numjobs = runlen*nodemax
    else:
        numjobs = runlen

    #submit first 250 jobs
    #monitor qstat every hour until all jobs are done   
    cwd = os.getcwd()
    jobid = np.empty(0)
    devnull = open(os.devnull, 'w')
    if runlen > 0:
        if numjobs <= 250:
            # Submit all jobs
            for i,s in enumerate(dirstring):
                os.chdir(s)
                jobid = np.append(jobid,subprocess.check_output(['qsub',runfile]))
                os.chdir(cwd)

            # Monitor jobs
            stat = subprocess.call(['qstat -u pcastell'], shell=True, stdout=devnull)
            while (stat == 0):
                print 'Waiting 1 minutes'
                time.sleep(60)
                stat = subprocess.call(['qstat -u pcastell'], shell=True, stdout=devnull)

            print 'All jobs done'
            # Check for errors
            # Clean things up
            for i,s in enumerate(jobid):
                s = s.strip('\n')
                errcheck = check_for_errors(dirstring[i],s,nodemax=nodemax)
                if (errcheck is False):
                    if (additional_output):
                        destroy_workspace(s,dirstring[i],outdirstring[i],
                                      addoutdir=addoutdirstring[i],nodemax=nodemax)
                    else:
                        destroy_workspace(s,dirstring[i],outdirstring[i],
                                      addoutdir=None,nodemax=nodemax)
                else:
                    print 'Jobid ',s,' exited with errors'

            # for i,s in enumerate(jobid):
            #     s = s.strip('\n')
            #     if nodemax is not None and nodemax > 1:
            #         # Loop through the nodes working on this job
            #         finished = 0
            #         for a in np.arange(nodemax):
            #             a = a + 1
            #             stat = subprocess.check_output(['squeue','-j',s+'_'+ str(a)])
            #             if (s+'_'+ str(a) not in stat):
            #                 print ' Job Not Found! '+s+'_'+ str(a)
            #                 finished = finished + 1

            #         if (finished == nodemax):
            #             #all nodes done for this job
            #             #clean up work space!
            #             errcheck = check_for_errors(dirstring[i],s,nodemax=nodemax)
            #             if (errcheck is False):
            #                 destroy_workspace(s,dirstring[i],outdirstring[i],
            #                                   addoutdir=addoutdirstring[i],nodemax=nodemax)

            #     else:
            #         stat = subprocess.call(['qstat',s], stdout=devnull)
            #         if (stat != 0):
            #             #all done for this job
            #             #clean up work space!
            #             errcheck = check_for_errors(dirstring[i],s,nodemax=nodemax)
            #             if (errcheck is False):
            #                 destroy_workspace(s,dirstring[i],outdirstring[i],
            #                                   addoutdir=addoutdirstring[i],nodemax=nodemax)                                  
        else:
            print 'too many jobs!  reduce run length'

    else:
        print 'No model hours to run'

    


