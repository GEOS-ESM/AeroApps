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
    DT_jobs  = 1    #days
    slurm    = 'aeronet_g5nr_sampler.j'
    tmp      = './tmp'

    parser = argparse.ArgumentParser()
    parser.add_argument("iso_t1",help='starting iso time')
    parser.add_argument("iso_t2",help='ending iso time')
    parser.add_argument("prep_config",
                        help="prep config filename")

    parser.add_argument('-D',"--DT_hours", default=DT_hours, type=int,
                        help="Timestep in hours for each file (default=%i)"%DT_hours)

    parser.add_argument("--DT_jobs", default=DT_jobs, type=int,
                        help="Timestep in days for each job (default=%i)"%DT_jobs)    

    parser.add_argument('-s',"--slurm",default=slurm,
                        help="slurm script template")           

    parser.add_argument('-t','--tmp',default=tmp,
                        help="temp directory")

    args = parser.parse_args()

    Date    = isoparser(args.iso_t1)
    enddate = isoparser(args.iso_t2)
    Dt      = timedelta(hours=args.DT_hours)
    Djobs   = timedelta(days=args.DT_jobs)

    if not os.path.exists(args.tmp):
        os.makedirs(args.tmp)

    cwd = os.getcwd()
    while Date <= enddate:
        edate = Date + Djobs - Dt
        if edate > enddate: edate = enddate

        # copy template to temp
        outfile = '{}_{}.j'.format(args.slurm[:-2],Date.isoformat())
        outpath = '{}/{}'.format(args.tmp,outfile)

        shutil.copyfile(args.slurm,outpath)

        f = open(outpath)
        text = []
        for l in f:
            text.append(l)

        iso1 = Date.isoformat()
        iso2 = edate.isoformat()

        newline = "#SBATCH --output=slurm_{}_%j.out\n".format(iso1)
        text[11] = newline

        newline = "#SBATCH --error=slurm_{}_%j.err\n".format(iso1)
        text[12] = newline

        newline = 'python -u $BIN -v --nproc 6 --DT_hours {} {} {} {} >'.format(args.DT_hours,iso1,iso2,args.prep_config) + ' tmp/slurm_{}_'.format(iso1)+'${SLURM_JOBID}_py.out\n'
        text[-2] = newline
        f.close()

        f = open(outpath,'w')
        for l in text:
            f.write(l)
        f.close()


        os.chdir(args.tmp)
        subprocess.call('sbatch {}'.format(outfile),shell=True)
        os.chdir(cwd)

        Date += Djobs

