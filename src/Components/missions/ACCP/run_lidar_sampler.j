#!/bin/csh -f

#######################################################################
#                     Batch Parameters for Run Job
#######################################################################
#SBATCH --time=12:00:00
#SBATCH --ntasks=12 --cpus-per-task=1 --ntasks-per-node=12
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
setenv G5DIR /discover/nobackup/pcastell/workspace/GAAS/src
setenv BIN $PWD 
setenv TMPDIR $LOCAL_TMPDIR 

mkdir -p $LOCAL_TMPDIR/sampling
cp -r /discover/nobackup/pcastell/workspace/GAAS/src/Components/missions/A-CCP/sampling/* $LOCAL_TMPDIR/sampling
ln -s $LOCAL_TMPDIR/sampling ./sampling

cp -r /discover/nobackup/pcastell/workspace/GAAS/src/Components/missions/A-CCP/tle/ $LOCAL_TMPDIR
ln -s $LOCAL_TMPDIR/tle ./tle

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
python -u run_lidar_sampler.py -v --nproc 12 --DT_hours 24 2006-01-01T00 2006-01-02T00 lidar.pcf > slurm_${SLURM_JOBID}_py.out
rm -rf $LOCAL_TMPDIR/*
rm sampling
rm tle
