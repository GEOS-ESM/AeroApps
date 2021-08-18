#!/bin/csh -f

#######################################################################
#                     Batch Parameters for Run Job
#######################################################################
#SBATCH --time=12:00:00
#SBATCH --constraint=sky
#SBATCH --ntasks=40 --cpus-per-task=1 --ntasks-per-node=40
#SBATCH --job-name=merge_random
#SBATCH --account=s2190
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
setenv G5BIN /discover/nobackup/pcastell/workspace/GAAS/Linux/bin
setenv BIN $PWD

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
python -u ./merge_random_mp.py  --ext all --fine 2006-01-01 2006-05-01 merge_files.pcf gpm.pcf polar07.pcf lidar_ext.pcf

#python -u ./merge_random_mp.py --merged all --levelb --surface --ext all --spc all --fine 2006-04-01 2006-05-01 merge_files.pcf gpm.pcf polar07.pcf lidar_ext.pcf
