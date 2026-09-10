#!/bin/csh
#SBATCH --time=1:00:00
#SBATCH --nodes=1 --ntasks-per-node=126
#SBATCH --job-name=cpexcv_sample_traj
#SBATCH --account=s3339

# setup environment
  setenv AEROAPPS /home/pcolarco/geos_aerosols/pcolarco/AeroApps
  source $AEROAPPS/env@/g5_modules
  setenv PYTHONPATH ./:$AEROAPPS/install/lib/Python/scat:$AEROAPPS/install/lib/Python/:$AEROAPPS/install/lib/Python/pyobs:$AEROAPPS/install/lib/Python/pyods
  set path = ( . $AEROAPPS/install/bin $path )
# ? Why do I need this next line?
  setenv LD_LIBRARY_PATH $AEROAPPS/install/lib/:$LD_LIBRARY_PATH


  ln -s /home/pcolarco/ExtData ExtData
  
# run sampling scripts
  set tracks = `\ls -1 /discover/nobackup/projects/gmao/iesa/aerosol/campaigns/CPEX-AW/*ict`

  prund.pl -H `hostname` -d `echo $tracks` &
#  mpirun -np 16 prund.pl -H `hostname` ./sample_trajectory.py m21c aer_inst_3hr_glo_Nv %s
#  mpirun -np 16 prund.pl -H `hostname` ./sample_trajectory.py MERRA2 inst3_3d_aer_Nv %s
  mpirun -np 16 prund.pl -H `hostname` ./sample_trajectory.py MERRA2 tavg1_2d_aer_Nx %s
#  mpirun -np 16 prund.pl -H `hostname` ./sample_trajectory.py fp inst3_3d_aer_Nv %s

