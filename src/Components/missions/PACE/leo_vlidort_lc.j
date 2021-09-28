#!/bin/csh -f

#######################################################################
#                     Batch Parameters for Run Job
#######################################################################
#SBATCH --time=05:00:00
#SBATCH --time-min=05:00:00
#SBATCH --ntasks=28 --cpus-per-task=1 --ntasks-per-node=28
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=patricia.castellanos@nasa.gov
#SBATCH --job-name=pace_vlidort
#SBATCH -A g0620
#SBATCH --output=slurm_%A.out
#SBATCH --error=slurm_%A.err
#######################################################################
#                  System Environment Variables
#######################################################################

umask 022

limit stacksize unlimited

#######################################################################
#           Architecture Specific Environment Variables
#######################################################################
setenv G5DIR /discover/nobackup/pcastell/workspace/GAAS/src
setenv GEOBIN /discover/nobackup/pcastell/workspace/GAAS/src/Components/missions/PACE 
setenv RUN_CMD  "mpirun -np 28"

source $HOME/.cshrc
cd $G5DIR
source g5_modules

#######################################################################
#   Move to Run Directory
#######################################################################
cd $GEOBIN

##################################################################
######
######         Perform single iteration of VLIDORT Run
###### Make sure to clean up shared memory segements before and after the job
######
##################################################################
./clean_mem.sh
$RUN_CMD  ./leo_vlidort_cloud.x leo_vlidort_cloud_550-discover.rc 
./clean_mem.sh 

