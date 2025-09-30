#!/bin/csh
#SBATCH --time=00:40:00
#SBATCH --nodes=1 --ntasks-per-node=125
#SBATCH --job-name=omi_simul
#SBATCH --constraint=mil
#SBATCH --account=s1148

# Set base location
  cd /home/pcolarco/ARCSIX/magpie

# setup environment
  setenv AEROAPPS $NOBACKUP/AeroApps
  source $AEROAPPS/env@/g5_modules
  setenv PYTHONPATH ./:$AEROAPPS/install/lib/Python/scat:$AEROAPPS/install/lib/Python/:$AEROAPPS/install/lib/Python/pyobs:$AEROAPPS/install/lib/Python/pyods
  set path = ( . $AEROAPPS/install/bin $path )
# ? Why do I need this next line?
  setenv LD_LIBRARY_PATH $AEROAPPS/install/lib/:$LD_LIBRARY_PATH

  ln -s /discover/nobackup/projects/gmao/share/dasilva/fvInput ExtData

# run sampling scripts
  set dates = ""
  foreach MM ("05" "06" "07" "08")
  foreach DD ("01" "02" "03" "04" "05" "06" "07" "08" "09" "10" \
              "11" "12" "13" "14" "15" "16" "17" "18" "19" "20" \
              "21" "22" "23" "24" "25" "26" "27" "28" "29" "30" "31")
   set dates = `echo $dates`" 2024"$MM$DD
  end
  end
  echo $dates
  prund.pl -H `hostname` -d `echo $dates` &
#  mpirun -np 128 prund.pl -H `hostname` ./sample_stations.py c180R_arcsix inst3d_aer_v %s
#  mpirun -np 63 prund.pl -H `hostname` ./sample_stations.py fp inst3_3d_aer_Nv %s
  mpirun -np 63 prund.pl -H `hostname` ./sample_stations.py MERRA2 tavg1_2d_aer_Nx

