#!/bin/sh

if [ $# -ne 4 ]; then
  echo "Usage: $0 [model] [sccyymmddhh] [eccyymmddhh] [hinc]"
  exit 1
else
  model=$1
  sccyymmddhh=$2
  eccyymmddhh=$3
  hinc=$4
fi

sdate=`echo $sccyymmddhh | cut -c1-8`
stime=`echo $sccyymmddhh | cut -c9,10`
stime=`expr $stime \* 10000`
edate=`echo $eccyymmddhh | cut -c1-8`
etime=`echo $eccyymmddhh | cut -c9,10`
etime=`expr $etime \* 10000`

dir=`dirname $0`

portal_dir=`cat $dir/wxmap.rc | grep '^portal_dir:' | cut -d' ' -f2`
local_dir=`cat $dir/wxmap.rc | grep '^local_dir:' | cut -d' ' -f2`

fcst_dt=`timetag $sdate $stime %Y%m%dT%H%M%S`
start_dt=$fcst_dt
end_dt=`timetag $edate $etime %Y%m%dT%H%M%S`
fnode=`timetag $sdate $stime %Y%m%dT%H%M%S`

# Panels
# ======

if [ "$model" == "ALL" ]; then

  path=$portal_dir/FIREX-AQ/%%Y%%m%%dT%%H%%M%%S/'ALL/$field/$level'
  name='model.ALL.%%Y%%m%%dT%%H%%M%%S.${tau}.$field.$level.$region.png'
  pathname=$path/$name

  wxd.py --theme custom_mission --theme FIREX-AQ --fcst_dt $fcst_dt --start_dt $fcst_dt --end_dt $end_dt --t_deltat $hinc --stream GEOS --oname $pathname --navigate 'off'

  exit 0
fi

# Datagrams
# =========

dest=$portal_dir/FIREX-AQ/$fnode/$model/datagrams
mkdir -p $dest

wxd.py --theme custom_mission --theme FX_$model --fcst_dt $fcst_dt --time_dt $fcst_dt --start_dt $fcst_dt --end_dt $fcst_dt --stream $model --region fx-central --oname $dest/'$field_name/model.'$model.$fnode'.${tau}.$field_name.0.$tag_name.png' --fields _gram

# Regional Plots
# ==============

path=$portal_dir/FIREX-AQ/%%Y%%m%%dT%%H%%M%%S/'$stream/$field/$level'
name='model.$stream.%%Y%%m%%dT%%H%%M%%S.${tau}.$field.$level.$region.png'
pathname=$path/$name

wxd.py --theme custom_mission --theme FX_$model --fcst_dt $fcst_dt --start_dt $fcst_dt --end_dt $end_dt --t_deltat $hinc --stream $model --oname $pathname --navigate 'off'
