#!/bin/csh -f

#######################################################################
#                     Batch Parameters for Run Job
#######################################################################
#SBATCH --time=03:00:00
#SBATCH --constraint=hasw
#SBATCH --ntasks=28 --cpus-per-task=1 --ntasks-per-node=28
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=patricia.castellanos@nasa.gov
#SBATCH --job-name=pace_vlidort
#SBATCH -A s1412
#SBATCH --array=1-2
#SBATCH --output=slurm_%A_%a.out
#SBATCH --error=slurm_%A_%a.err

#######################################################################
#                    USER INPUTS
#######################################################################
setenv AEROBIN $PWD 


#######################################################################
#                  System Environment Variables
#######################################################################
umask 022

limit stacksize unlimited

#######################################################################
#           Architecture Specific Environment Variables
#######################################################################
setenv G5DIR /discover/nobackup/pcastell/workspace/GAAS/src
setenv RUN_CMD  "mpirun -np 28"

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
./clean_mem.sh
$RUN_CMD ./leo_vlidort_cloud.x leo_vlidort_cloud.rc ${SLURM_ARRAY_TASK_ID}
./clean_mem.sh

