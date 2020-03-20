#!/bin/csh -x
#
#PBS -N PlumeRise
#PBS -W group_list=g0604
#PBS -l select=16
#PBS -l walltime=8:00:00
#PBS -j oe
#
#------------------------------------------------------------------------
#
# Runs the PlumeRise Model for a month
#
#------------------------------------------------------------------------

  set YYYYMM_END = 200409

  @ ncpus = 16 * 4

# Hardwired directories for now
# -----------------------------
  set ARCH = `uname -s`
  set exp_dir = /discover/nobackup/dasilva/pr_v1/
  set met_dir = /discover/nobackup/dasilva/ceres_b55/diag
  set bin_dir = /home/dasilva/workspace/prsa/$ARCH/bin 
  set path = ( $bin_dir $path )

# Get year/month for this segment
# -------------------------------
  cd $exp_dir/run
  set when = ( `cat date.rst` )
  set year  = $when[1]
  set month = $when[2]

# Start timer
# -----------
  /bin/rm -rf .zeit
  zeit_ci.x PlumeRise

# Stop here if end of experiment already
# --------------------------------------
  @ yyyymm = 100 * $year + $month
  if ( $yyyymm > $YYYYMM_END ) then
    echo "pr: end of experiment reached; nothing to do"
    exit 1
  endif

# Start server
# ------------
  set files = ( `find $met_dir/Y$year/M$month -name '*.diag.eta.*.nc' -print` )
  $bin_dir/prund.pl -d -H $HOSTNAME $files &

# Start clients under MPI
# -----------------------
  mpirun -np $ncpus $bin_dir/prund.pl -H $HOSTNAME $bin_dir/pr_instance.csh |& tee pr.log

# Next segment
# ------------
  @ month++
  if ( $month > 12 ) then
    @ month = 1
    @ year++
  endif

# Save "restart" date
# -------------------
  if ( $month < 10 ) then
    echo $year 0$month > date.rst
  else
    echo $year $month  > date.rst
  endif

  @ yyyymm = 100 * $year + $month
  if ( $yyyymm > $YYYYMM_END ) then
    echo "pr: end of experiment reached; no more jobs"
  else
    qsub -N PR_$yyyymm pr.j  # resubmit itself 
  endif

# Report timing
# -------------
  zeit_co.x PlumeRise
  zeit_pr.x 

# All done
# --------
  exit 0


