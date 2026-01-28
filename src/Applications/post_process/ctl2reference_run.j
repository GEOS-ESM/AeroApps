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

source $SRC_DIR/env@/g5_modules


#######################################################################
#          Run CTL2KERCHUNK 
#######################################################################

# Check if this is a continuation job (APPEND flag file exists)
set APPEND_FLAG = ""
if ( -f .ctl2reference_append_flag ) then
    set APPEND_FLAG = "--append"
endif

# Run the script
# modify max_batches to fit into the 12 hour run time. will depend on file sizes and disk
python3 -u ./ctl2reference.py ctl2reference.yaml --max_batches=13 ${APPEND_FLAG} >& ctl2reference-${SLURM_JOB_ID}.log

# Capture the exit status
set EXIT_STATUS = $status

# Check the exit status
if ( $EXIT_STATUS == 2 ) then
    # Work is complete, clean up flag file
    echo "Work complete. Exit status = 2"
    if ( -f .ctl2reference_append_flag ) then
        rm .ctl2reference_append_flag
    endif
else
    # Work not complete, create flag file and resubmit
    echo "Work incomplete. Exit status = ${EXIT_STATUS}. Resubmitting..."
    touch .ctl2reference_append_flag
    sbatch $0
endif
