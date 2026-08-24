#!/bin/csh
#SBATCH --time=1:00:00
#SBATCH --nodes=1 --ntasks-per-node=125
#SBATCH --job-name=inspyre_sample_traj
#SBATCH --account=s3339

# Set base location
  cd /home/pcolarco/INSPYRE_25/pcolarco/sample

# setup environment
  setenv AEROAPPS $NOBACKUP/AeroApps
  source $AEROAPPS/env@/g5_modules
  setenv PYTHONPATH ./:$AEROAPPS/install/lib/Python/scat:$AEROAPPS/install/lib/Python/:$AEROAPPS/install/lib/Python/pyobs:$AEROAPPS/install/lib/Python/pyods
  set path = ( . $AEROAPPS/install/bin $path )
# ? Why do I need this next line?
  setenv LD_LIBRARY_PATH $AEROAPPS/install/lib/:$LD_LIBRARY_PATH

# run sampling scripts
  set tracks = `\ls -1 ../data/*RA.ict`
  echo $tracks
  exit
  prund.pl -H `hostname` -d `echo $tracks` &
  mpirun -np 48 prund.pl -H `hostname` ./sample_trajectory.py fp inst3_3d_aer_Nv %s
#  mpirun -np 48 prund.pl -H `hostname` ./sample_trajectory.py fp tavg1_2d_lfo_Nx %s
#  mpirun -np 48 prund.pl -H `hostname` ./sample_trajectory.py MERRA2 inst3_3d_aer_Nv %s

