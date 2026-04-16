#!/bin/csh -f

#######################################################################
#                     Batch Parameters for Run Job
#######################################################################
 
#SBATCH -J castnet_sampler
#SBATCH --nodes=1
#SBATCH --constraint=mil
#SBATCH --time=12:00:00
#SBATCH -A @GROUPID
#SBATCH -o output_castnet_sampler-%j.log
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL
##SBATCH --qos=debug
#######################################################################
#  Run sampler code for ASIA-AQ model runs at AMON surface stations
#######################################################################
#           Architecture Specific Environment Variables
#######################################################################

setenv SRC_DIR @SRCDIR
setenv PYTHONPATH ${SRC_DIR}/install/lib/Python

source $SRC_DIR/env@/g5_modules


#######################################################################
#          Run Sampler 
#######################################################################

python3 -u ./castnet_sampler.py sampling.yaml >& castnet_sampler-${SLURM_JOB_ID}.log
