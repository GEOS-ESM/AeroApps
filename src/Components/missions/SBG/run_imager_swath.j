#!/bin/csh -f

#######################################################################
#                     Batch Parameters for Run Job
#######################################################################
#SBATCH --time=12:00:00
#SBATCH --ntasks=120 --cpus-per-task=1 --ntasks-per-node=120
#SBATCH --constraint=mil
#SBATCH --job-name=lidar_sampler
#SBATCH -A s2190
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=patricia.castellanos@nasa.gov
#SBATCH --output=slurm_%j.out
#SBATCH --error=slurm_%j.err
#######################################################################
#                  System Environment Variables
#######################################################################

umask 022

limit stacksize unlimited

#######################################################################
#           Architecture Specific Environment Variables
#######################################################################
source $HOME/.cshrc
source ./setup_env


##################################################################
######
######         Do Sampling
######
##################################################################
python -u imager_swath.py -v 2006-01-16T00 2006-01-17T00 -d 200 --dt_units milliseconds -D 5 --out_year 2006 imager_files.pcf ssd650.pcf sbg.pcf > slurm_${SLURM_JOBID}_py.out
