#!/bin/csh -f

#######################################################################
#                     Batch Parameters for Run Job
#######################################################################
 
#SBATCH -J ctl2reference
#SBATCH --nodes=1
#SBATCH --constraint=mil
#SBATCH --time=12:00:00
#SBATCH -A @GROUPID
#SBATCH -o output_ctl2reference-%j.log
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL
##SBATCH --qos=debug
#######################################################################
#  Run sampler code for ASIA-AQ model runs at IMPROVE surface stations
#######################################################################
#           Architecture Specific Environment Variables
#######################################################################

setenv SRC_DIR @SRCDIR
setenv PYTHONPATH ${SRC_DIR}/install/lib/Python
setenv PATH ${PATH}:${SRC_DIR}/install/bin

source $SRC_DIR/env@/g5_modules


#######################################################################
#          Run CTL2KERCHUNK 
#######################################################################

python3 -u ./ctl2reference.py sampling.yaml >& ctl2reference-${SLURM_JOB_ID}.log
