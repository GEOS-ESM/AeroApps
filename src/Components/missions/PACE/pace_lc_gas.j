#!/bin/csh -f

#######################################################################
#                     Batch Parameters for Run Job
#######################################################################
#SBATCH --time=12:00:00
#SBATCH --constraint=sky
#SBATCH --ntasks=40 --cpus-per-task=1 --ntasks-per-node=40
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=patricia.castellanos@nasa.gov
#SBATCH --job-name=pace
#SBATCH --account=s2190
#SBATCH --output=slurm_%j.out
#SBATCH --error=slurm_%j.err

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
setenv RUN_CMD  "mpirun -np 40"

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
$RUN_CMD ./pace_vlidort_gasabs_multinode.x pace_vlidort.rc
./clean_mem.sh
