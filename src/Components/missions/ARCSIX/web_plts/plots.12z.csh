#!/bin/tcsh

cd /discover/nobackup/pcolarco/AeroApps/src/Components/missions/ARCSIX/web_plts

# Get the current date
set YYYY  = `date '+20%y'` 
set MM    = `date '+%m'`
set DD    = `date '+%d'`
set HH    = "12"

# Target for a flight (edit)
set YYYT = $YYYY
set MT   = $MM
set DT   = $DD

# Get the end date
set YYYY_  = `date --date='+5 days' '+20%y'` 
set MM_    = `date --date='+5 days' '+%m'`
set DD_    = `date --date='+5 days' '+%d'`
set HH_    = "12"

echo $YYYY$MM$DD
echo $YYYY_$MM_$DD_

set start = $YYYY$MM${DD}T$HH
set end   = $YYYY_$MM_${DD_}T$HH_

set web = "/discover/nobackup/projects/gmao/iesa/pub/aerosol/arcsix/$start"
mkdir -p $web
mkdir -p $web/plots
mkdir -p $web/datagrams
chmod 775 $web $web/plots $web/datagrams
mkdir -p $start
mkdir -p $start/plots
mkdir -p $start/datagrams

# Make a web page
sed "s/YYYY/$YYYY/g;s/MM/$MM/g;s/DD/$DD/g;s/HH/$HH/g" index.tmp > $start/index.html

# ===============
# python plotting
# Update buoy locations
python buoys.py
\cp -f buoys.csv buoys.csv.$YYYY$MM${DD}_$HH
sed "s/YYYY/'$YYYY'/g;s/MM/'$MM'/g;s/DD/'$DD'/g;s/HH/'$HH'/g" arcsix_plots.csh > arcsix_plots.csh.$YYYY$MM${DD}_${HH}
qsub arcsix_plots.csh.$YYYY$MM${DD}_$HH
#rsync -auPv $start/ $web
#exit
# ===============





# Make the curtain plots
#module load python/GEOSpyD/Ana2019.10_py3.7
sed "s/YYYY/$YYYY/g;s/MM/$MM/g;s/DD/$DD/g;s/HH/$HH/g;s/YYYT/$YYYT/g;s/MT/$MT/g;s/DT/$DT/g" curtains.tmp > curtains.py
#sed "s/YYYY/$YYYY/g;s/MM/$MM/g;s/DD/$DD/g" curtains.12z.tmp > curtains.py
python curtains.py
\mv -f alert*png $start
\mv -f svalb*png $start
\mv -f green*png $start

#source ~jardizzo/src/FLUID/beta/utils/pyg_modules
source /discover/nobackup/jardizzo/software/GEOS-ESM/ARCSIX/WxMap/utils/pyg_modules
module load ffmpeg

# Get lists of what's plottable
#wxlist.py --theme chem2d_mission --theme ARCSIX
#wxlist.py --theme chem3d_mission --theme ARCSIX
#wxlist.py --theme custom_mission --theme ARCSIX | grep "du_"
#wxlist.py --theme weather_mission --theme ARCSIX

# Grab the current datagrams
# Don't know how to do just yet....

# Get the datagrams
set dir = /discover/nobackup/dao_ops/FLUID/datagrams/images/fp/Y$YYYY/M$MM/D$DD/H$HH
echo $dir
set locs = ("_74.7051239014_-94.9694061279" "_76.5161056519_-68.7689971924" \
            "_83.62_-31.2" "_80.05_-86.42" "_82.4508_-62.5072" \
            "_74.6_-67.4" "_87.0_87.0" "_78.92_0.0" "_87.0_28.0" \
            "_88.0_-60.0" "_81.6_-16.66" "_72.58_-38.46" "_78.9_11.88")
set site = ( Resolute_Bay Pituffik Kaffeklubben_Island Eureka Alert\
             Baffin_Bay Belle Fram_Strait Lisee Northern_Arctic_Ocean \
             Station_Nord Summit_Greenland Svalbard_Zeppelin)
foreach iloc (1 2 3 4 5 6 7 8 9 10 11 12 13 )
\cp -f $dir/meteo$locs[$iloc].png $start/datagrams/meteo_$site[$iloc].png
\cp -f $dir/oc$locs[$iloc].png $start/datagrams/oc_$site[$iloc].png
\cp -f $dir/ocmass$locs[$iloc].png $start/datagrams/ocmass_$site[$iloc].png
\cp -f $dir/co$locs[$iloc].png $start/datagrams/co_$site[$iloc].png
\cp -f $dir/total$locs[$iloc].png $start/datagrams/total_$site[$iloc].png
end

exit
# Below is FLUID plots
# Get the total aod
foreach field (totaot ssaot ocaot)
wxmap.py --theme chem2d_public --stream G5FPFC \
         --fcst_dt $start --start_dt $start --end_dt $end \
         --field $field --region nps --level 0 \
         --oname 'model.NASA_GEOS-5.%%Y%%m%%d%%H00.${tau}.$field.nps.png'
rm -f *.nav
./mp4_encode -o model.NASA_GEOS-5.$start.$field.nps.mp4 "model.NASA_GEOS-5.$YYYY$MM$DD${HH}00.???."${field}".nps.png"
\mv -f *png $start/plots
end

# Get the cloud and precip
wxmap.py --theme weather_mission --theme ARCSIX --stream G5FPFC \
         --fcst_dt $start --start_dt $start --end_dt $end \
         --field cldlow --region arcsix --level 0 \
         --oname 'model.NASA_GEOS-5.%%Y%%m%%d%%H00.${tau}.CloudLow.ARCSIX.png'
rm -f *.nav
./mp4_encode -o model.NASA_GEOS-5.$start.CloudLow.ARCSIX.mp4 "model.NASA_GEOS-5.$YYYY$MM$DD${HH}00.???.CloudLow.ARCSIX.png"
\mv -f *png $start/plots

# Get the PM by level
foreach lev (500 600 850)
wxmap.py --theme chem3d_public --stream G5FPFC \
         --fcst_dt $start --start_dt $start --end_dt $end \
         --field oc --region nps --level $lev \
         --oname 'model.NASA_GEOS-5.%%Y%%m%%d%%H00.${tau}.OC.'${lev}'hpa.nps.png'
rm -f *.nav
./mp4_encode -o model.NASA_GEOS-5.$start.OC.${lev}hpa.nps.mp4 "model.NASA_GEOS-5.$YYYY$MM$DD${HH}00.???.OC."${lev}"hpa.nps.png"
\mv -f *png $start/plots
end


# Get the CO by level
foreach lev (500 600 850)
wxmap.py --theme chem3d_public --stream G5FPFC \
         --fcst_dt $start --start_dt $start --end_dt $end \
         --field co --region nps --level $lev \
         --oname 'model.NASA_GEOS-5.%%Y%%m%%d%%H00.${tau}.CO.'${lev}'hpa.nps.png'
rm -f *.nav
./mp4_encode -o model.NASA_GEOS-5.$start.CO.${lev}hpa.nps.mp4 "model.NASA_GEOS-5.$YYYY$MM$DD${HH}00.???.CO."${lev}"hpa.nps.png"
\mv -f *png $start/plots
end

# Get the COBBNA by level
foreach lev (500 600 850)
wxmap.py --theme chem3d_public --stream G5FPFC \
         --fcst_dt $start --start_dt $start --end_dt $end \
         --field cobbna --region nps --level $lev \
         --oname 'model.NASA_GEOS-5.%%Y%%m%%d%%H00.${tau}.COBBNA.'${lev}'hpa.nps.png'
rm -f *.nav
./mp4_encode -o model.NASA_GEOS-5.$start.COBBNA.${lev}hpa.nps.mp4 "model.NASA_GEOS-5.$YYYY$MM$DD${HH}00.???.COBBNA."${lev}"hpa.nps.png"
\mv -f *png $start/plots
end

# Get the total aod
foreach field (totaot ssaot ocaot)
wxmap.py --theme chem2d_mission --theme ARCSIX --stream G5FPFC \
         --fcst_dt $start --start_dt $start --end_dt $end \
         --field $field --region arcsix --level 0 \
         --oname 'model.NASA_GEOS-5.%%Y%%m%%d%%H00.${tau}.$field.ARCSIX.png'
rm -f *.nav
./mp4_encode -o model.NASA_GEOS-5.$start.$field.ARCSIX.mp4 "model.NASA_GEOS-5.$YYYY$MM$DD${HH}00.???."${field}".ARCSIX.png"
\mv -f *png $start/plots
end

chmod g+w *mp4
\mv -f *mp4 $start
\cp -f python/fp.*.$YYYY$MM${DD}_${HH}z.mp4 $start/
\cp -f python/fp.*.$YYYY$MM${DD}_${HH}+$YYYY$MM${DD}_${HH}00.png $start/
rsync -auPv $start/ $web


exit

