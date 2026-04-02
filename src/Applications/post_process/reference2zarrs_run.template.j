#!/bin/csh -f

#######################################################################
#                     Batch Parameters for Run Job
#######################################################################
#SBATCH --job-name=zarr_VARNAME
#SBATCH --time=WALLTIME
#SBATCH --constraint=mil
#SBATCH --nodes=1
#SBATCH -A @GROUPID
#SBATCH --output=LOG_DIR/zarr_VARNAME_%j.log
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL
##SBATCH --qos=debug
#######################################################################
#           Architecture Specific Environment Variables
#######################################################################

setenv SRC_DIR @SRCDIR
setenv PYTHONPATH ${SRC_DIR}/install/lib/Python

source $SRC_DIR/env@/g5_modules

#######################################################################
#          Run reference2zarrs
#######################################################################


python SCRIPT CONFIG \
    --varname VARNAME \
    --n_workers N_WORKERS \
    --chunk_time CHUNK_TIME \
    --log_file LOG_DIR/zarr_VARNAME_${SLURM_JOB_ID}.log


