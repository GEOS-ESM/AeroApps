#!/bin/csh -f

#######################################################################
#                     Batch Parameters for Run Job
#######################################################################
 
#SBATCH -J aaq_sampler
#SBATCH --nodes=1
#SBATCH --constraint=mil
#SBATCH --time=2:00:00
#SBATCH -A @GROUPID
#SBATCH -o output_aaq_sampler-%j.log
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL
##SBATCH --qos=debug
#######################################################################
#  Run sampler code for ASIA-AQ
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

python3 -u ./aaq_g3_hsrl_sampler.py sampling_m2_g3_hsrl.yaml >& aaq_sampler-${SLURM_JOB_ID}.log
#python3 -u ./aaq_derived_hsrl.py sampling_m2_g3_hsrl.yaml >& aaq_derived-${SLURM_JOB_ID}.log
