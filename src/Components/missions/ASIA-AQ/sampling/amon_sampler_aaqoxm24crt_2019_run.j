#!/bin/csh -f

#######################################################################
#                     Batch Parameters for Run Job
#######################################################################
 
#SBATCH -J amon_sampler
#SBATCH --nodes=1
#SBATCH --constraint=mil
#SBATCH --time=6:00:00
#SBATCH -A s2942
#SBATCH -o output_amon_sampler-%j.log
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL
##SBATCH --qos=debug
#######################################################################
#  Run sampler code for ASIA-AQ model runs at AMON surface stations
#######################################################################
#           Architecture Specific Environment Variables
#######################################################################

setenv SRC_DIR /gpfsm/dnb34/pcastell/workspace/AeroApps_asia_aq/AeroApps
setenv PYTHONPATH ${SRC_DIR}/install/lib/Python

source $SRC_DIR/env@/g5_modules


#######################################################################
#          Run Sampler 
#######################################################################
if (! -d ExtData) then
    ln -s /home/pcastell/opendap/dasilva_fvinput/ExtData/ .
endif

python3 -u ./amon_sampler.py sampling_aaqoxm24crt_2019.yaml >& amon_sampler-${SLURM_JOB_ID}.log
