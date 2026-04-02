#!/bin/tcsh
#
# reference2zarrs_run.csh
# Reads variables from parquet reference store, fills in slurm template,
# and submits one job per variable running up to MAX_JOBS concurrently
#

# Usage: reference2zarrs_run.csh <config.yaml>
#
# Reads variables from parquet reference store, fills in slurm template,
# and submits one job per variable running up to MAX_JOBS concurrently
#

# =====================================================================
# Check config file argument
# =====================================================================
if ($#argv != 1) then
    echo "Usage: $0 <config.yaml>"
    echo ""
    echo "  config.yaml  -- configuration yaml file"
    exit 1
endif

set CONFIG = $argv[1]

if (! -f $CONFIG) then
    echo "ERROR: Config file not found: $CONFIG"
    exit 1
endif

echo "Using config: $CONFIG"


# =====================================================================
# Configuration — edit these
# =====================================================================
set SCRIPT     = "reference2zarrs.py"
set TEMPLATE   = "reference2zarrs_run.template.j"
set LOG_DIR    = "logs/zarr_conversion"
set SLURM_DIR  = "slurm_scripts"
set MAX_JOBS   = 10
set N_WORKERS  = 20
set CHUNK_TIME = 1
set WALLTIME   = "3:00:00"
set PYTHON     = "python3"

# =====================================================================
# Check template exists
# =====================================================================
if (! -f $TEMPLATE) then
    echo "ERROR: SLURM template not found: $TEMPLATE"
    exit 1
endif

# =====================================================================
# Get list of variables from the parquet reference store
# =====================================================================

# write a temp python script to get variable list
set TMP_PY = "/tmp/get_vars_$$.py"
cat >! $TMP_PY << EOF
import yaml, fsspec, xarray as xr, sys

config = yaml.safe_load(open("$CONFIG"))
com_outdir = config['combined_references_dir']
ctl = config['reference_ctl']

ref_parq = f"{com_outdir}/reference_store_{ctl}.parq"
try:
    fs = fsspec.filesystem("reference", fo=ref_parq,
                           remote_protocol="file", lazy=True)
    ds = xr.open_dataset(fs.get_mapper(""), engine="zarr",
                         chunks={"time": 1}, consolidated=False)
    for vn in list(ds.data_vars)[0:1]:
        print(vn)
except Exception as e:
    print(f"ERROR: {e}", file=sys.stderr)
    sys.exit(1)
EOF

set VARIABLES = (`$PYTHON $TMP_PY`)
set exit_status = $status

# clean up temp file
rm -f $TMP_PY

if ($exit_status != 0) then
    echo "ERROR: Failed to get variable list from parquet reference store"
    exit 1
endif

if ($#VARIABLES == 0) then
    echo "ERROR: No variables found in parquet reference store"
    exit 1
endif

echo "Found $#VARIABLES variables:"
echo "  $VARIABLES"
echo ""

# =====================================================================
# Create log and slurm script directories
# =====================================================================
mkdir -p $LOG_DIR
mkdir -p $SLURM_DIR

# =====================================================================
# Submit jobs — max MAX_JOBS running concurrently
# =====================================================================
set n_submitted = 0
set n_total     = $#VARIABLES
set job_ids     = ()

echo "Submitting jobs (max $MAX_JOBS concurrent)..."
echo ""

foreach vn ($VARIABLES)

    # count currently running/pending jobs
    @ n_running = `squeue -u $USER -h -r | grep "zarr_" | wc -l`

    # wait if at max concurrent jobs
    while ($n_running >= $MAX_JOBS)
        echo "  $n_running jobs running -- waiting 30s for a slot..."
        sleep 30
        @ n_running = `squeue -u $USER -h -r | grep "zarr_" | wc -l`
    end

    # === Fill in the template for this variable ===
    set SLURM_SCRIPT = "$SLURM_DIR/zarr_${vn}.j"

    sed -e "s|VARNAME|${vn}|g" \
        -e "s|N_WORKERS|${N_WORKERS}|g" \
        -e "s|LOG_DIR|${LOG_DIR}|g" \
        -e "s|SCRIPT|${SCRIPT}|g" \
        -e "s|CONFIG|${CONFIG}|g" \
        -e "s|CHUNK_TIME|${CHUNK_TIME}|g" \
        -e "s|WALLTIME|${WALLTIME}|g" \
        $TEMPLATE > $SLURM_SCRIPT

    # submit the filled-in slurm script
    set JOB_ID = `sbatch --parsable $SLURM_SCRIPT`

    if ($status != 0) then
        echo "ERROR: Failed to submit job for $vn"
        exit 1
    endif

    set job_ids = ($job_ids $JOB_ID)
    @ n_submitted++

    echo "  [$n_submitted/$n_total] Submitted $vn -- job ID: $JOB_ID"
    echo "                          script: $SLURM_SCRIPT"

end

# =====================================================================
# Summary
# =====================================================================
echo ""
echo "======================================================="
echo "Submission Summary"
echo "======================================================="
echo "  Variables submitted: $n_submitted / $n_total"
echo "  Max parallel jobs:   $MAX_JOBS"
echo "  Log dir:             $LOG_DIR"
echo "  Slurm scripts:       $SLURM_DIR"
echo ""
echo "  Job IDs: $job_ids"
echo ""
echo "  Monitor with:"
echo "    squeue -u $USER"
echo ""
echo "  Check logs with:"
echo "    tail -f $LOG_DIR/zarr_AIRDENS.log"
echo ""
echo "  View generated slurm script with:"
echo "    cat $SLURM_DIR/zarr_AIRDENS.slurm"
echo ""
echo "  If some variables fail, rerun just those with:"
echo "    sbatch $SLURM_DIR/zarr_AIRDENS.slurm"
echo "======================================================="
