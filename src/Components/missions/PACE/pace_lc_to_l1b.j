#!/bin/csh -f

#######################################################################
#                     Batch Parameters for Run Job
#######################################################################
#SBATCH --time=2:00:00
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

./pace_lc_to_l1b.py --do_single_xtrack --no_write_cld --no_write_aer --run_name morel_f0_1_nosleave_noaerosol --ALPHA_file alphaTable_v0/alpha_CK_Thuillier_o3_F0_1.nc4 --IRR_file irrTable_v0/F0_1_rsr_weighted_V0.nc  2006-03-24T00:50 

