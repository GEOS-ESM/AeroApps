#!/bin/csh -f

#######################################################################
#                     Batch Parameters for Run Job
#######################################################################
#SBATCH --time=02:00:00
#SBATCH --ntasks=24 --cpus-per-task=1 --ntasks-per-node=24
#SBATCH --job-name=stn_sampler
#SBATCH -A s1412
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=patricia.castellanos@nasa.gov
#SBATCH --output=slurm_%j.out
#######################################################################
#                  System Environment Variables
#######################################################################

umask 022

limit stacksize unlimited

#######################################################################
#           Architecture Specific Environment Variables
#######################################################################
setenv G5DIR /discover/nobackup/pcastell/workspace/GAAS/src
setenv BIN /discover/nobackup/pcastell/workspace/GAAS/src/Components/stn_vlidort 

source $HOME/.cshrc
cd $G5DIR
source g5_modules

#######################################################################
#   Move to Run Directory
#######################################################################
cd $BIN

##################################################################
######
######         Do Sampling
######
##################################################################
python -u run_stn_sampler.py -v --nproc 24 --DT_hours 24 2006-01-01T00 2006-01-02T00 stn_sampler.pcf > run_stn_sampler_%j.out

