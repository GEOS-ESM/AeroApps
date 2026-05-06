#!/bin/csh
#SBATCH --time=01:30:00
#SBATCH --nodes=1 --ntasks-per-node=42
#SBATCH --job-name=arcsix_plots
#SBATCH --account=s3339

  source ~/.cshrc

  setenv AEROAPPS /home/pcolarco/geos_aerosols/pcolarco/AeroApps/
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
  cd ROOTDIR/$start

# run sampling scripts
   set files = `ls -1 /home/pcolarco/fp/forecast/Y$yy/M$mm/D$dd/H${hh}/GEOS.fp.fcst.inst1_2d_hwl_Nx.*.V01.nc4`
   prund.pl -H `hostname` -d `echo $files` &
   mpirun -np 42 prund.pl -H `hostname` simple_plot.py %s

# run sampling scripts
   set files = `ls -1 /home/pcolarco/fp/forecast/Y$yy/M$mm/D$dd/H${hh}/GEOS.fp.fcst.inst1_2d_hwl_Nx.*.V01.nc4`
   prund.pl -H `hostname` -d `echo $files` &
   mpirun -np 42 prund.pl -H `hostname` simple_plot.cloud.py %s

# run sampling scripts
   set files = `ls -1 /home/pcolarco/fp/forecast/Y$yy/M$mm/D$dd/H${hh}/GEOS.fp.fcst.inst3_3d_aer_Nv.*.V01.nc4`
   prund.pl -H `hostname` -d `echo $files` &
   mpirun -np 42 prund.pl -H `hostname` sample_trajectory.py %s


   
module load ffmpeg

# Post process
  foreach COLL (TOTEXTTAU BREXTTAU OCEXTTAU SSEXTTAU DUEXTTAU SUEXTTAU \
                TOTEXTTAU_wcloud BREXTTAU_wcloud OCEXTTAU_wcloud SSEXTTAU_wcloud \
		DUEXTTAU_wcloud SUEXTTAU_wcloud)
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

  foreach COLL (lon80W lon100W lon110W lon120W lat040N lat050N)
  set i = 0
  foreach file (`\ls -1 *$COLL*png`)
   set iii = $i
   if($i < 10)  set iii = '0'$iii
   if($i < 100) set iii = '0'$iii
   ln -s $file fp.$COLL.$iii.png
   @ i = $i + 1
  end
  ./mp4_encode -o fp.$COLL.$yy$mm${dd}_${hh}z.mp4 "fp.$COLL.???.png"
  rm -f fp.$COLL.???.png
  end  
  
set web = "/discover/nobackup/projects/gmao/iesa/pub/aerosol/inspyre/$start"

\cp -f index.html $web
\cp -f fp.*.$yy$mm${dd}_${hh}z.mp4 $web/
\cp -f fp.*.$yy$mm${dd}_${hh}+$yy$mm${dd}_${hh}00.png $web/
\cp -f fp.*.$yy$mm${dd}_${hh}+*_0000.png $web/plots
\cp -f fp.*.$yy$mm${dd}_${hh}+*_1200.png $web/plots
mkdir -p ROOTDIR/samples/INSPYRE/sampled/GEOS-FP/$start
\mv -f samples/INSPYRE/sampled/GEOS-FP/$start/* ROOTDIR/samples/INSPYRE/sampled/GEOS-FP/$start
\rm -f *.png *.mp4
