#!/bin/csh -f

#######################################################################
#                     Batch Parameters for Run Job
#######################################################################
 
#SBATCH -J ctl2kerchunk
#SBATCH --nodes=1
#SBATCH --constraint=mil
#SBATCH --time=12:00:00
#SBATCH -A s2942
#SBATCH -o output_ctl2kerchunk-%j.log
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL
##SBATCH --qos=debug
#######################################################################
#  Run sampler code for ASIA-AQ model runs at IMPROVE surface stations
#######################################################################
#           Architecture Specific Environment Variables
#######################################################################

setenv SRC_DIR /gpfsm/dnb34/pcastell/workspace/AeroApps_asia_aq/AeroApps
setenv PYTHONPATH ${SRC_DIR}/install/lib/Python

source $SRC_DIR/env@/g5_modules


#######################################################################
#          Run CTL2KERCHUNK 
#######################################################################

python3 -u ./ctl2kerchunk.py sampling.yaml >& ctl2kerchunk-${SLURM_JOB_ID}.log
