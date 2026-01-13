#!/bin/csh -f

#######################################################################
#                     Batch Parameters for Run Job
#######################################################################
 
#SBATCH -J improve_sampler
#SBATCH --nodes=1
#SBATCH --constraint=mil
#SBATCH --time=12:00:00
#SBATCH -A @GROUPID
#SBATCH -o output_aaq_sampler-%j.log
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
#          Run Sampler 
#######################################################################
if (! -d ExtData) then
    ln -s /home/pcastell/opendap/dasilva_fvinput/ExtData/ .
endif

python3 -u ./improve_sampler.py sampling.yaml >& improve_sampler-${SLURM_JOB_ID}.log
python3 -u ./improve_derived.py sampling.yaml >& improve_derived-${SLURM_JOB_ID}.log
