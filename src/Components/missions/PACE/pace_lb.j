#!/bin/csh -f

#######################################################################
#                     Batch Parameters for Run Job
#######################################################################
#SBATCH --time=12:00:00
#SBATCH --constraint=sky
#SBATCH --ntasks=40 --cpus-per-task=1 --ntasks-per-node=40
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=patricia.castellanos@nasa.gov
#SBATCH --job-name=sample
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

source $HOME/.cshrc
cd $G5DIR
source g5_modules

#######################################################################
#   Move to Run Directory
#######################################################################
cd $AEROBIN

##################################################################
######
######         Run python code
######
##################################################################

./pace_lb.py pace_lb.20.pcf &
./pace_lb.py pace_lb.25.pcf &
./pace_lb.py pace_lb.30.pcf &
./pace_lb.py pace_lb.35.pcf &
./pace_lb.py pace_lb.40.pcf &
./pace_lb.py pace_lb.45.pcf &
./pace_lb.py pace_lb.50.pcf &
./pace_lb.py pace_lb.55.pcf &

wait

