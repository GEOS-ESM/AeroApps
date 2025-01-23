#!/bin/csh -f

#######################################################################
#                     Batch Parameters for Run Job
#######################################################################
#SBATCH --time=12:00:00
#SBATCH --ntasks=120 --cpus-per-task=1 --ntasks-per-node=120
#SBATCH --job-name=swath
#SBATCH --constraint=mil
#SBATCH -A s2190
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=patricia.castellanos@nasa.gov
#SBATCH --output=slurm_%j.out
#SBATCH --error=slurm_%j.err
#######################################################################
#
#                  System Environment Variables
#
#######################################################################

umask 022

limit stacksize unlimited

#######################################################################
#
#           Architecture Specific Environment Variables
#
#######################################################################
setenv BIN $PWD 
setenv TMPDIR $LOCAL_TMPDIR 
setenv RUNDIR /discover/nobackup/pcastell/workspace/aero_work/aist/aos-sky

mkdir -p ${LOCAL_TMPDIR}/tle
cp -r ${RUNDIR}/tle/ ${LOCAL_TMPDIR}
ln -s ${LOCAL_TMPDIR}/tle ./tle

source $HOME/.cshrc
source ./setup_env


##################################################################
######
######         Do Calculations
######
##################################################################
#python -u run_accp_polarimeter_swath.py -v --DT_hours 1 2006-01-01T00:00 2006-01-02T00:00 lidar_files.pcf gpm.pcf polar07.pcf  > slurm_${SLURM_JOBID}_py.out &
#rm -rf $LOCAL_TMPDIR/*
#rm tle
