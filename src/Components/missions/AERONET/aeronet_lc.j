#!/bin/csh -f

#######################################################################
#                     Batch Parameters for Run Job
#######################################################################
#SBATCH --time=03:00:00
#SBATCH --constraint=hasw
#SBATCH --ntasks=28 --cpus-per-task=1 --ntasks-per-node=28
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=patricia.castellanos@nasa.gov
#SBATCH --job-name=aeronet_vlidort
#SBATCH -A s1412
#SBATCH --output=slurm_%j.out
#SBATCH --error=slurm_%j.err

#######################################################################
#                    USER INPUTS
#######################################################################
setenv START 2006-01-01T00:00:00 
setenv END   2006-01-01T01:00:00
setenv OPTIONS "-D 1"
setenv AEROBIN /discover/nobackup/pcastell/workspace/GAAS/src/Components/missions/AERONET 


#######################################################################
#                  System Environment Variables
#######################################################################
umask 022

limit stacksize unlimited

#######################################################################
#           Architecture Specific Environment Variables
#######################################################################
setenv G5DIR /discover/nobackup/pcastell/workspace/GAAS/src
setenv RUN_CMD  "./aeronet_vlidort_lc.py"

source $HOME/.cshrc
cd $G5DIR
source g5_modules

#######################################################################
#   Move to Run Directory
#######################################################################
cd $AEROBIN

##################################################################
######
######         Perform single iteration of VLIDORT Run
######
##################################################################
$RUN_CMD $OPTIONS $START $END > slurm_${SLURM_JOBID}_py.out
