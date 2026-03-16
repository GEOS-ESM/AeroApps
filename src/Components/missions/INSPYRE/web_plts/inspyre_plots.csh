#!/bin/csh
#SBATCH --time=00:40:00
#SBATCH --nodes=1 --ntasks-per-node=40
#SBATCH --job-name=arcsix_plots
#SBATCH --account=s1148

  source ~/.cshrc

  setenv AEROAPPS $NOBACKUP/AeroApps/
  source $AEROAPPS/env@/g5_modules
  setenv PYTHONPATH ./:$AEROAPPS/install/lib/Python/scat:$AEROAPPS/install/lib/Python/:$AEROAPPS/install/lib/Python/pyobs:$AEROAPPS/install/lib/Python/pyods
  set path = ( . $AEROAPPS/install/bin $path )
# ? Why do I need this next line?
  setenv LD_LIBRARY_PATH $AEROAPPS/install/lib/:$LD_LIBRARY_PATH

#  module load python/GEOSpyD/Ana2019.10_py3.7

  which python
  
  set yy = YYYY
  set mm = MM
  set dd = DD
  set hh = HH

  set start = $yy$mm${dd}T$hh
  cd $NOBACKUP/AeroApps/src/Components/missions/INSPYRE/web_plts/$start

# run sampling scripts
   set files = `ls -1 /home/pcolarco/fp/forecast/Y$yy/M$mm/D$dd/H${hh}/GEOS.fp.fcst.inst1_2d_hwl_Nx.*.V01.nc4`
   prund.pl -H `hostname` -d `echo $files` &
   mpirun -np 40 prund.pl -H `hostname` simple_plot.n.py %s

# run sampling scripts
   set files = `ls -1 /home/pcolarco/fp/forecast/Y$yy/M$mm/D$dd/H${hh}/GEOS.fp.fcst.inst1_2d_hwl_Nx.*.V01.nc4`
   prund.pl -H `hostname` -d `echo $files` &
   mpirun -np 40 prund.pl -H `hostname` simple_plot.cloud.py %s

module load ffmpeg

# Post process
  foreach COLL (TOTEXTTAU OCEXTTAU SSEXTTAU DUEXTTAU SUEXTTAU \
                TOTEXTTAU_wcloud OCEXTTAU_wcloud SSEXTTAU_wcloud DUEXTTAU_wcloud SUEXTTAU_wcloud)
  set i = 0
  foreach file (`\ls -1 fp.$COLL.$yy$mm${dd}_$hh*png`)
   set iii = $i
   if($i < 10)  set iii = '0'$iii
   if($i < 100) set iii = '0'$iii
   ln -s $file $file:r:r.$iii.png
   @ i = $i + 1
  end
  ./mp4_encode -o fp.$COLL.$yy$mm${dd}_${hh}z.mp4 "fp.$COLL.???.png"
  rm -f fp.$COLL.???.png
  end

set web = "/discover/nobackup/projects/gmao/iesa/pub/aerosol/inspyre/$start"

\cp -f fp.*.$yy$mm${dd}_${hh}z.mp4 $web/
\cp -f fp.*.$yy$mm${dd}_${hh}+$yy$mm${dd}_${hh}00.png $web/
\cp -f fp.*.$yy$mm${dd}_${hh}+*_0000.png $web/plots
\cp -f fp.*.$yy$mm${dd}_${hh}+*_1200.png $web/plots
\rm -f *.png *.mp4
