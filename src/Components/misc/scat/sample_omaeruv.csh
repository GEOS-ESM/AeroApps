#!/bin/csh
#SBATCH --time=2:00:00
#SBATCH --nodes=1 --ntasks-per-node=40
#SBATCH --constraint=sky
#SBATCH --job-name=sample_omaeruv
#SBATCH --account=s2323

# set OMI Level 2 information
  set inp = '/home/pcolarco/piesa/data/OMI/Level2/OMAERUV/'

# setup experiment
  set EXPID = c180R_v202_aura_schill
  cd $NOBACKUP/aura/${EXPID}

# setup environment
  setenv AEROAPPS $NOBACKUP/AeroApps/AeroApps
  source $AEROAPPS/env@/g5_modules
  setenv PYTHONPATH ./:$AEROAPPS/install/lib/Python/scat:$AEROAPPS/install/lib/Python/:$AEROAPPS/install/lib/Python/pyobs:$AEROAPPS/install/lib/Python/pyods
  set path = ( . $AEROAPPS/install/bin $path )
# ? Why do I need this next line?
  setenv LD_LIBRARY_PATH $AEROAPPS/install/lib/:$LD_LIBRARY_PATH

# setup outputs
  ln -s /gpfsm/dnb04/projects/p22/aerosol/experiments/colarco/$EXPID/inst3d_aer_v/ inst3d_aer_v

# prund try
#  setenv HYDRA_LAUNCHER_EXTRA_ARGS "--input none"
   set files = `ls -1 /home/pcolarco/piesa/data/OMI/Level2/OMAERUV/2016/09/??/OMI-Aura_L2-OMAERUV_2016m09*he5`
   prund.pl -H `hostname` -d `echo $files` &
   mpirun -np 20 prund.pl -H `hostname` ./sample_omaeruv.py -F --verbose -x $EXPID -d inst3d_aer_v %s inst3d_aer_v.ddf
