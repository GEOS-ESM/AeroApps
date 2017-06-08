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

if __name__ == '__main__':
    
    #Defaults
    DT_hours = 24
    slurm    = 'run_lidar_sampler.j'
    tmp      = './tmp'

    parser = argparse.ArgumentParser()
    parser.add_argument("iso_t1",help='starting iso time')
    parser.add_argument("iso_t2",help='ending iso time')

    parser.add_argument('-D',"--DT_hours", default=DT_hours, type=int,
                        help="Timestep in hours for each file (default=%i)"%DT_hours)

    parser.add_argument('-s',"--slurm",default=slurm,
                        help="slurm script template")           

    parser.add_argument('-t','--tmp',default=tmp,
                        help="temp directory")

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
        newline = 'python -u run_lidar_sampler.py -v --nproc 6 --DT_hours 24 {} {} lidar.pcf >'.format(iso1,iso2) + ' slurm_${SLURM_JOBID}_py.out\n'
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

