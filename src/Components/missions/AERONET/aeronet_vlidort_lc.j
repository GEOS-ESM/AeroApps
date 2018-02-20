#!/bin/csh -f

#######################################################################
#                     Batch Parameters for Run Job
#######################################################################
#SBATCH --time=02:00:00
#SBATCH --constraint=hasw
#SBATCH --ntasks=28 --cpus-per-task=1 --ntasks-per-node=28
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=patricia.castellanos@nasa.gov
#SBATCH --job-name=aeronet_vlidort
#SBATCH -A s1412
#SBATCH --output=slurm_%A.out
#SBATCH --error=slurm_%A.err

#######################################################################
#                    USER INPUTS
#######################################################################
setenv START 2006-01-01T00:00:00 
setenv END   2006-01-01T01:00:00


#######################################################################
#                  System Environment Variables
#######################################################################

umask 022

limit stacksize unlimited

#######################################################################
#           Architecture Specific Environment Variables
#######################################################################
setenv G5DIR /discover/nobackup/pcastell/workspace/GAAS/src
setenv AEROBIN /discover/nobackup/pcastell/workspace/GAAS/src/Components/missions/AERONET 
setenv RUN_CMD  "./do_aeronet.py -D 1"

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
$RUN_CMD  $START $END 
