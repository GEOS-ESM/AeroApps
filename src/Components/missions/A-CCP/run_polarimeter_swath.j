#!/bin/csh -f

#######################################################################
#                     Batch Parameters for Run Job
#######################################################################
#SBATCH --time=12:00:00
#SBATCH --ntasks=8 --cpus-per-task=1 --ntasks-per-node=8
#SBATCH --job-name=swath
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
python -u run_accp_polarimeter_swath.py -v --DT_hours 1 2006-01-01T00:00 2006-01-02T00:00 lidar_files.pcf gpm.pcf polar07.pcf  > slurm_${SLURM_JOBID}_py.out &
