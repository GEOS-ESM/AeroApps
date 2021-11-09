#!/bin/csh
#SBATCH --time=12:00:00
#SBATCH --nodes=2 --ntasks-per-node=40
#SBATCH --job-name=omi_simul
#SBATCH --constraint=sky
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
  set OUTPATH = ./out/
  mkdir -p $OUTPATH

# run sampling scripts
   set files = `ls -1 /home/pcolarco/piesa/data/OMI/Level2/OMAERUV/2016/09/??/OMI-Aura_L2-OMAERUV_2016m09??t1[0,1,2,3]*he5`
   prund.pl -H `hostname` -d `echo $files` &
   mpirun -np 20 prund.pl -H `hostname` ./omi_simul.py %s

