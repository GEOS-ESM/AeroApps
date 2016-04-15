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


def make_workspace(date,ch,code,outdir,runfile,instname,prefix='workdir',nodemax=None,
                   addoutdir=None,layout=None):
    cwd     = os.getcwd()
    dirname = prefix + '/'+ instname.lower() + '.' + str(date.date())+'T'+str(date.hour).zfill(2)+'.'+ch+'.'+code        
    jobname = instname.lower() + '.' + str(date.date())+'T'+str(date.hour).zfill(2)+'.'+ch
    if layout is not None:
        dirname = dirname + '.' + layout
        jobname = jobname + '.' + layout

    outdir = outdir + '/Y'+str(date.year)+'/M'+str(date.month).zfill(2)+'/D'+str(date.day).zfill(2)
    if addoutdir is not None:
        addoutdir = addoutdir + '/Y'+str(date.year)+'/M'+str(date.month).zfill(2)+'/D'+str(date.day).zfill(2)

    bindir = dirname
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

        elif (line[0:13] == 'setenv GEOBIN'):
            destination.write('setenv GEOBIN '+bindir+'\n')

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


def destroy_workspace(jobid,dirname,outdir,addoutdir=None,nodemax=None,profile=False):
    cwd     = os.getcwd()
    os.chdir(dirname)

    os.remove('Aod_EOS.rc')
    os.remove('Chem_MieRegistry.rc')
    os.remove('geo_vlidort.x')
    os.remove('clean_mem.sh')
    os.remove('ExtData')
    os.remove('geo_vlidort.rc')
    os.remove('geo_vlidort_run_array.j')

    if profile is False:
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
    if profile is False:
        os.rmdir(dirname)

def combine_files(filelist):
    mergedfile = filelist[0]
    parts      = mergedfile.split('.')
    mergedfile = ''
    for p in parts[0:-2]:
        mergedfile = mergedfile + p + '.'
    mergedfile = mergedfile + 'nc4'
    #mergedfile = mergedfile[:-5] + 'nc4'

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


def make_maiac_rcfile(dirname,indir,date,ch,code,interp,additional_output, instname,
                      nodemax=None,i_band=None,version=None,surf_version=None,
                      layout=None):
    cwd = os.getcwd()
    os.chdir(dirname)

    rcfile = open('geo_vlidort.rc','w')
    rcfile.write('INDIR: '+indir+'\n')
    rcfile.write('OUTDIR: .\n')
    rcfile.write('DATE: '+str(date.year)+str(date.month).zfill(2)+str(date.day).zfill(2)+'\n')
    rcfile.write('TIME: '+str(date.hour).zfill(2)+'\n')
    rcfile.write('INSTNAME: ' + instname.lower() + '\n')
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
        rcfile.write('SURFBAND_C: 645 858 469 555 1240 1640 2130 412\n')
    else:
        rcfile.write('SURFBAND_I: '+ i_band + '\n')

    rcfile.write('SURFBANDM: 8 \n')

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

    if version is not None:
        rcfile.write('VERSION: '+version+'\n')

    if surf_version is not None:
        rcfile.write('SURF_VERSION: '+surf_version+'\n')

    if layout is not None:
        rcfile.write('LAYOUT: '+layout+'\n')
        
    rcfile.close()

    os.chdir(cwd)

def make_mcd43c_rcfile(dirname,indir,date,ch,code,interp,additional_output, instname,
                      nodemax=None,i_band=None,version=None,surf_version=None,
                      layout=None):
    cwd = os.getcwd()
    os.chdir(dirname)

    rcfile = open('geo_vlidort.rc','w')
    rcfile.write('INDIR: '+indir+'\n')
    rcfile.write('OUTDIR: .\n')
    rcfile.write('DATE: '+str(date.year)+str(date.month).zfill(2)+str(date.day).zfill(2)+'\n')
    rcfile.write('TIME: '+str(date.hour).zfill(2)+'\n')
    rcfile.write('INSTNAME: ' + instname.lower() + '\n')
    rcfile.write('SURFNAME: MCD43C\n')

    #figure out correct MODIS doy
    if (str(startdate.date()) == '2005-12-31'):
        rcfile.write('SURFDATE: 2005361\n')
    elif (str(startdate.date()) == '2006-01-01'):
        rcfile.write('SURFDATE: 2005361\n')
    elif ((str(startdate.year) == '2007') and (str(date.month).zfill(2) == '03')):
        rcfile.write('SURFDATE: 2007081\n')
    elif ((str(startdate.year) == '2007') and (str(date.month).zfill(2) == '04')):
        rcfile.write('SURFDATE: 2007097\n')
    elif ((str(startdate.year) == '2006') and (str(date.month).zfill(2) == '07')):
        rcfile.write('SURFDATE: 2006201\n')
    elif ((str(startdate.year) == '2006') and (str(date.month).zfill(2) == '08')):
        rcfile.write('SURFDATE: 2006201\n')        
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
        rcfile.write('SURFBAND_I: '+ i_band + '\n')

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

    if version is not None:
        rcfile.write('VERSION: '+version+'\n')

    if surf_version is not None:
        rcfile.write('SURF_VERSION: '+surf_version+'\n')

    if layout is not None:
        rcfile.write('LAYOUT: '+layout+'\n')
        
    rcfile.close()

    os.chdir(cwd)    

def make_ler_rcfile(dirname,indir,date,ch,code,interp,additional_output, instname,
                    nodemax=None,i_band=None,version=None,surf_version=None,layout=None):
    cwd = os.getcwd()
    os.chdir(dirname)

    rcfile = open('geo_vlidort.rc','w')
    rcfile.write('INDIR: ' + indir + '\n')
    rcfile.write('OUTDIR: .\n')
    rcfile.write('DATE: ' + str(date.year) + str(date.month).zfill(2) + str(date.day).zfill(2) + '\n')
    rcfile.write('TIME: ' + str(date.hour).zfill(2) + '\n')
    rcfile.write('INSTNAME: ' + instname + '\n')
    rcfile.write('SURFNAME: LER\n')
    rcfile.write('SURFMODEL: LER\n')

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

    if version is not None:
        rcfile.write('VERSION: '+version+'\n')

    if surf_version is not None:
        rcfile.write('SURF_VERSION: '+surf_version+'\n')

    if layout is not None:
        rcfile.write('LAYOUT: '+layout+'\n')

    rcfile.close()

    os.chdir(cwd)    

def prefilter(date,indir,instname,layout=None):
    g5dir = indir + '/LevelB/'+ 'Y'+ str(date.year) + '/M' + str(date.month).zfill(2) + '/D' + str(date.day).zfill(2) 
    nymd  = str(date.year) + str(date.month).zfill(2) + str(date.day).zfill(2)
    hour  = str(date.hour).zfill(2)

    if layout is None:
        met   = g5dir + '/' + instname.lower() + '-g5nr.lb2.met_Nv.' + nymd + '_' + hour + 'z.nc4'
        geom  = g5dir + '/' + instname.lower() + '.lb2.angles.' + nymd + '_' + hour + 'z.nc4'
        land  = indir + '/LevelB/invariant/' + instname.lower() + '-g5nr.lb2.asm_Nx.nc4'  
    else:
        met   = g5dir + '/' + instname.lower() + '-g5nr.lb2.met_Nv.' + nymd + '_' + hour + 'z_' + laycode +'.nc4'
        geom  = g5dir + '/' + instname.lower() + '.lb2.angles.' + nymd + '_' + hour + 'z_' + laycode +'.nc4'
        land  = indir + '/LevelB/invariant/' + instname.lower() + '-g5nr.lb2.asm_Nx_' + laycode + '.nc4'  


    ncMet = Dataset(met)
    Cld   = np.squeeze(ncMet.variables[u'CLDTOT'][:])
    mask  = Cld.mask
    Cld   = Cld[~mask]    
    f     = np.where(Cld <= 0.01)
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

    SZA = SZA[~mask]
    VZA = VZA[~mask]
    FRLAND = FRLAND[~mask]
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
    instname          = 'sentinel-4'
    version           = '1.0'    
    startdate         = '2006-07-27T00:00:00'
    enddate           = '2006-07-27T00:00:00'
    episode           = None
    channels          = '550'
    surface           = 'MCD43C'
    interp            = 'interpolate'
    i_band            = '7'    
    additional_output = False
    nodemax           = 1
    surf_version      = 'MCD43C'
    layout            = None

    runfile           = 'geo_vlidort_run_array.j'
    nccs              = '/discover/nobackup/projects/gmao/osse2/pub/c1440_NR/OBS/'+ \
                         instname.upper() + '/DATA/'

    prefix            = nccs + 'workdir/'
    ################
    ###
    #    End of uper inputs
    ###
    ################
    profile           = False
    runmode           = 'scalar'
    indir             = nccs 
    outdir            = nccs + 'LevelC2'
    if (additional_output):
        addoutdir         = nccs + 'LevelC'
    else:
        addoutdir         = None
    

    if episode is not None:
        if episode == 1:
            startdate  = '2005-12-31T00:00:00Z'
            enddate    = '2006-01-01T23:00:00Z'
        elif episode == 2:
            startdate  = '2006-07-27T00:00:00Z'
            enddate    = '2006-08-09T23:00:00Z'
        elif episode == 3:
            startdate  = '2007-03-28T00:00:00Z'
            enddate    = '2007-03-29T23:00:00Z'
        elif episode == 4:
            startdate  = '2007-04-10T00:00:00Z'
            enddate    = '2007-04-11T23:00:00Z'

    dt = timedelta(hours=1)
    startdate = parse(startdate)
    enddate   = parse(enddate)
    jobsmax   = 150
    ##################################################
    ####
    #    Create working directories, SLURM scripts, and RC-files
    ####
    ##################################################
    dirstring    = np.empty(0)
    outdirstring = np.empty(0)
    nodemax_list = np.empty(0,dtype='int')
    if (additional_output):
        addoutdirstring = np.empty(0)

    if type(channels) is str:
        channels = channels.split()

    if type(i_band) is str:
        i_band = i_band.split()

    # Loop over dates
    while (startdate <= enddate):

        # check for layout keyword. 
        # figure out number of tiles        
        if layout is None:
            ntiles = 1
        else:
            ntiles = int(layout[0])*int(layout[1])

        # loop through tiles
        for tile in np.arange(ntiles):
        #for tile in [2]:
            if layout is not None:
                laycode = layout + str(tile)
            else:
                laycode = None

            # check to see if there is any work to do
            anypixels, numpixels = prefilter(startdate,indir,instname,layout=laycode) 
            if (anypixels):

                if (numpixels <= 1000 and nodemax is not None):
                    nodemax_ = 1
                else:
                    nodemax_ = nodemax

                for i, ch in enumerate(channels):
                    band_i = None
                    if (i_band is not None):
                        band_i = i_band[i]

                    # Vector Case
                    code = runmode + '.'
                    if (interp.lower() == 'interpolate'):
                        code = code + 'i'                 

                    code = code + surface

                    dirlist = make_workspace(startdate,ch,code,outdir,runfile,instname,nodemax=nodemax_,
                                             addoutdir=addoutdir,layout=laycode,prefix=prefix)

                    if (additional_output):
                        workdir, outdir_, addoutdir_ = dirlist
                    else:
                        workdir, outdir_ = dirlist

                    if surface.upper() == 'MAIACRTLS':
                        make_maiac_rcfile(workdir,indir,startdate,ch,runmode, interp,additional_output,
                                          instname, nodemax=nodemax_,i_band=band_i,
                                          version=version,surf_version=surf_version,layout=laycode)
                    elif surface.upper() == 'MCD43C':
                        make_mcd43c_rcfile(workdir,indir,startdate,ch,runmode, interp,additional_output,
                                          instname, nodemax=nodemax_,i_band=band_i,
                                          version=version,surf_version=surf_version,layout=laycode)
                    else:
                        make_ler_rcfile(workdir,indir,startdate,ch,runmode, interp,additional_output,
                                        instname, nodemax=nodemax_,i_band=band_i,
                                        version=version,surf_version=surf_version,layout=laycode)
                    
                    dirstring    = np.append(dirstring,workdir)
                    outdirstring = np.append(outdirstring,outdir_)
                    nodemax_list = np.append(nodemax_list,int(nodemax_))
                    if (additional_output):
                        addoutdirstring = np.append(addoutdirstring, addoutdir_)

       
        startdate = startdate + dt


    ##################################################
    ####
    #    Submit & Handle Jobs  
    ####
    ##################################################

    runlen  = len(dirstring)   
    if nodemax is not None:
        numjobs = sum(nodemax_list)
    else:
        numjobs = runlen

    cwd = os.getcwd()
    jobid = np.empty(0)
    devnull = open(os.devnull, 'w')
    if runlen > 0:
        if numjobs <= jobsmax:   
            countRun   = runlen    
            node_tally   = numjobs  
        else:
            if nodemax is not None:
                keep_adding = True
                node_tally = 0
                countRun = 0
                while(keep_adding):
                    if (node_tally + nodemax_list[countRun] <= jobsmax):
                        node_tally = node_tally + nodemax_list[countRun]
                        countRun = countRun + 1
                    else:
                        keep_adding = False                
            else:
                countRun = jobsmax
                node_tally = jobsmax

        workingJobs = np.arange(countRun)

        # Submit first JOBSMAX jobs
        for i in workingJobs:
            s = dirstring[i]
            os.chdir(s)
            jobid = np.append(jobid,subprocess.check_output(['qsub',runfile]))
            
        os.chdir(cwd)
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
                if nodemax is not None and nodemax_list[i] > 1:
                    # Loop through the nodes working on this job
                    finishedCNT = 0
                    for a in np.arange(nodemax_list[i]):
                        a = a + 1

                        try:
                            result = subprocess.check_output(['squeue','-j',s+'_'+ str(a)],stderr=subprocess.STDOUT)
                            if (s+'_'+ str(a) not in result):
                                #print ' Job Not Found! '+s+'_'+ str(a)
                                finishedCNT = finishedCNT + 1
                        except subprocess.CalledProcessError as e:
                            #print ' ERROR Job Not Found! '+s+'_'+ str(a)
                            finishedCNT = finishedCNT + 1

                    if (finishedCNT == nodemax_list[i]): 
                        finished = True
                else:
                    result = subprocess.call(['qstat',s], stdout=devnull)
                    if (result != 0):
                        finished = True   

                # if the job is finished clean up the workspace
                if finished:
                    #print 'Job finished, cleaning up', s, i 
                    finishedJobs = np.append(finishedJobs,ii)
                    errcheck = check_for_errors(dirstring[i],s,nodemax=nodemax_list[i])               
                    if (errcheck is False):
                        if (additional_output):
                            destroy_workspace(s,dirstring[i],outdirstring[i],
                                          addoutdir=addoutdirstring[i],nodemax=nodemax_list[i],profile=profile)
                        else:
                            destroy_workspace(s,dirstring[i],outdirstring[i],
                                          addoutdir=None,nodemax=nodemax_list[i],profile=profile)                       
                    else:
                        print 'Jobid ',s,' exited with errors'

            # finished checking up on all the jobs
            # Remove finished jobs from the currently working list
            if len(finishedJobs) != 0:
                print 'deleting finishedJobs',finishedJobs,jobid[workingJobs[finishedJobs]]
                if nodemax is not None:
                    node_tally  = node_tally - sum(nodemax_list[workingJobs[finishedJobs]])
                else:
                    node_tally  = node_tally - len(finishedJobs)

                workingJobs = np.delete(workingJobs,finishedJobs)

            # Add more jobs if needed
            # reinitialize stat variable
            if (runlen > countRun) and (node_tally < jobsmax):
                #print 'adding new jobs'
                if nodemax is not None:
                    keep_adding = True
                    newRun      = 0
                    while (keep_adding):
                        if (node_tally + nodemax_list[countRun + newRun] <= jobsmax):
                            node_tally = node_tally + nodemax_list[countRun + newRun]
                            newRun = newRun + 1
                            if (countRun + newRun == runlen):
                                keep_adding = False
                        else:
                            keep_adding = False
                else:
                    newRun = jobsmax - node_tally
                    node_tally = jobsmax

                newjobs  = countRun + np.arange(newRun)
                workingJobs = np.append(workingJobs, newjobs)
                for i in newjobs:
                    s = dirstring[i]
                    os.chdir(s)
                    jobid = np.append(jobid,subprocess.check_output(['qsub',runfile]))

                os.chdir(cwd)
                countRun = countRun + newRun
                stat = subprocess.call(['qstat -u pcastell'], shell=True, stdout=devnull)


            print 'Waiting 1 minutes'
            time.sleep(60)
            

        # Exited while loop
        print 'All jobs done'

    else:
        print 'No model hours to run'

    


