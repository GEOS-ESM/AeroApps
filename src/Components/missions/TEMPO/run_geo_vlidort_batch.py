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

def make_workspace(date,ch,code,outdir):
    dirname = str(date.date())+'T'+str(date.time())+'.'+ch+'.'+code
    outdir = outdir + '/Y'+str(date.year)+'/M'+str(date.month).zfill(2)+'/D'+str(date.day).zfill(2)
    if not os.path.exists(dirname):
        os.makedirs(dirname)
    
    os.chdir(dirname)
    if not os.path.isfile('geo_vlidort.x'):
        os.symlink('../geo_vlidort.x','geo_vlidort.x')
    if not os.path.exists('ExtData'):
        os.symlink('../ExtData','ExtData')
    if not os.path.isfile('Chem_MieRegistry.rc'):
        os.symlink('../Chem_MieRegistry.rc','Chem_MieRegistry.rc')
        
    source = open('../geo_vlidort_run.j','r')
    destination = open('geo_vlidort_run.j','w')
    for line in source:
        if (line[0:18] == '#SBATCH --job-name'):
            destination.write('#SBATCH --job-name '+dirname+'\n')
        elif (line[0:15] == 'setenv TEMPOBIN'):
            destination.write('setenv TEMPOBIN '+dirname+'\n')
        elif (line[0:13] == 'setenv OUTDIR'):            
            destination.write('setenv OUTDIR '+outdir+'\n')            
            if not os.path.exists(outdir):
                mkpath(outdir)
        else:
            destination.write(line)        
    source.close()
    destination.close()
    
    source = open('../Aod_EOS.rc','r')
    destination = open('Aod_EOS.rc','w')
    a = float(ch)*1e-3
    for line in source:
        if (line[0:11] == 'r_channels:'):
            destination.write('r_channels: '+'{:0.3f}e-6'.format(a)+'\n')
        else:
            destination.write(line) 
    source.close()
    destination.close()
     
    os.chdir('../') 
    
    return dirname 

def make_rcfile(dirname,indir,date,ch,code):
    os.chdir(dirname)

    rcfile = open('geo_vlidort.rc','w')
    rcfile.write('INDIR: '+indir+'\n')
    rcfile.write('OUTDIR: .\n')
    rcfile.write('DATE: '+str(date.year)+str(date.month).zfill(2)+str(date.day).zfill(2)+'\n')
    rcfile.write('TIME: '+str(date.hour).zfill(2)+'\n')
    rcfile.write('INSTNAME: tempo\n')
    rcfile.write('BRDFNAME: MAIACRTLS\n')

    #figure out correct MODIS doy
    if (str(startdate.date()) == '2005-12-31'):
        rcfile.write('BRDFDATE: 2006008\n')
    else:
        doy = date.toordinal() - datetime(date.year-1,12,31).toordinal()
        DOY = 8*(int(doy/8) + 1)
        if (DOY > 365):
            DOY = 8
            rcfile.write('BRDFDATE: '+str(date.year+1)+str(DOY).zfill(3)+'\n')
        else:
            rcfile.write('BRDFDATE: '+str(date.year)+str(DOY).zfill(3)+'\n')

    rcfile.write('BRDFBAND: interpolate\n')
    rcfile.write('BRDFBAND_C: 645 858 469 555 1240 1640 2130\n')
    if (code == 'scalar'):
        rcfile.write('SCALAR: true\n')
    else:
        rcfile.write('SCALAR: false\n')
    rcfile.write('CHANNELS: '+ch+'\n')
    rcfile.write('CLDMAX: 0.01\n')
    rcfile.write('SZAMAX: 90.0\n')
    rcfile.close()

    os.chdir('../')
        
#########################################################

if __name__ == "__main__":
    
    startdate = '2005-12-31T00:00:00'
    enddate   = '2005-12-31T02:00:00'
    channels  = '550','670'

    calculon  = '/nobackup/TEMPO'
    nccs      = '/discover/nobackup'
    
    if (os.path.exits(calculon)):
        indir     = calculon + '/DATA'
        outdir    = '/nobackup/3/pcastell/TEMPO/DATA'
    else:
        indir     = nccs + '/projects/gmao/osse2/pub/c1440_NR/OBS/TEMPO/DATA'
        outdir    = nccs + '/pcastell/TEMPO/DATA'
    
    dt = timedelta(hours=1)
    startdate = parse(startdate)
    enddate   = parse(enddate)

    runstring = np.empty(0)
    while (startdate <= enddate):
        for ch in channels:
            # Vector Case
            workdir = make_workspace(startdate,ch,'vector',outdir)
            make_rcfile(workdir,indir,startdate,ch,'vector')
            runstring = np.append(runstring,'./'+workdir+'/geo_vlidort_run.j')

            # Scalar Case
            workdir = make_workspace(startdate,ch,'scalar',outdir)
            make_rcfile(workdir,indir,startdate,ch,'scalar')
            runstring = np.append(runstring,'./'+workdir+'/geo_vlidort_run.j')
       
        startdate = startdate + dt
     
    runlen  = len(runstring)   
    streams = min(10, runlen)
    drun    = np.array([runlen/streams for i in range(streams)])
    drun[streams-1] = drun[streams-1] + math.fmod(runlen,streams)
    for s,i in enumerate(np.linspace(0,streams*(runlen/streams),streams,endpoint=False)):
        for j in range(drun[s]):
            if (j == 0 ):
                jobid = int(subprocess.check_output(['qsub',runstring[i+j]]))
                #print 'qsub '+runstring[i+j]
            else:
                jobid = int(subprocess.check_output(['qsub','-W','depend=afterok:'+str(jobid),runstring[i+j]]))
                #print 'qsub -W depend=afterok:str(jobid)' +' '+runstring[i+j]

    


