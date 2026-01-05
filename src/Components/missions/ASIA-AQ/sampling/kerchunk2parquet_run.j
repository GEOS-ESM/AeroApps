#!/bin/csh -f

#######################################################################
#                     Batch Parameters for Run Job
#######################################################################
 
#SBATCH -J kerchunk2parquet
#SBATCH --nodes=1
#SBATCH --constraint=mil
#SBATCH --time=6:00:00
#SBATCH -A @GROUPID
#SBATCH -o output_kerchunk2parquet-%j.log
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

source $SRC_DIR/env@/g5_modules


#######################################################################
#          Run CTL2KERCHUNK 
#######################################################################

python3 -u ./kerchunk2parquet.py sampling.yaml >& kerchunk2parquet-${SLURM_JOB_ID}.log
