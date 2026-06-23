#!/bin/csh -f

#######################################################################
#                     Batch Parameters for Run Job
#######################################################################
 
#SBATCH -J ict_sampler
#SBATCH --nodes=1
#SBATCH --constraint=mil
#SBATCH --time=2:00:00
#SBATCH -A s2191
#SBATCH -o output_ict_sampler-%j.log
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL
##SBATCH --qos=debug
#######################################################################
#  Run sampler code for Field Campaign ICARTT Files
#######################################################################
#           Architecture Specific Environment Variables
#######################################################################

setenv SRC_DIR /gpfsm/dnb34/acollow/AeroApps/AeroApps
setenv PYTHONPATH ${SRC_DIR}/install/lib/Python

source $SRC_DIR/env@/g5_modules


#######################################################################
#          Run Sampler 
#######################################################################
if (! -d ExtData) then
    ln -s /home/pcastell/opendap/dasilva_fvinput/ExtData/ .
endif

python3 -u ./samplecampaignforict.py sampling_firex.yaml #>& ict_sampler-${SLURM_JOB_ID}.log
python3 -u ./getderivedfields.py sampling_firex.yaml #>& ict_derived-${SLURM_JOB_ID}.log
