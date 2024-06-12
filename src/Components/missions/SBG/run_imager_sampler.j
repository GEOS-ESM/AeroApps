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
source ./setup_env

setenv TMPDIR $LOCAL_TMPDIR 

mkdir -p $LOCAL_TMPDIR/sampling
cp -r $AERODIR/src/Components/missions/SBG/sampling/* $LOCAL_TMPDIR/sampling
ln -s $LOCAL_TMPDIR/sampling ./sampling


##################################################################
######
######         Do Sampling
######
##################################################################
python -u run_imager_sampler.py -v --nproc 120 2006-01-15T00 2006-01-16T00 imager_files.pcf imager_sampler.pcf > slurm_${SLURM_JOBID}_py.out
rm -rf $LOCAL_TMPDIR/*
rm sampling
