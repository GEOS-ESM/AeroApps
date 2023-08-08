#!/bin/sh

if [ $# -ne 2 ]; then
  echo "Usage: $0 [ccyymmdd] [hhmmss]"
  exit 1
else
  idate=$1
  itime=$2
fi

exec >/dev/null 2>&1

SENTINEL=/discover/nobackup/projects/gmao/nca/pub/firexaq/opendap/GEOSFP/inst1_2d_hwl_Nx.%Y%m%d_%H

while [ 1 ]; do

  hour=`expr $itime / 10000`
  if [ $hour -eq 0 ]; then
    flen=05
  else
    flen=05
  fi

  edattim=`timetag $idate $itime {%Y%m%d_%H%M%S}%+d$flen`
  edate=`echo $edattim | cut -d'_' -f1`
  etime=`echo $edattim | cut -d'_' -f2`
  sentinel=`timetag $idate $itime $SENTINEL`
  sentinel=`timetag $edate $etime $sentinel`

  if [ -f $sentinel ]; then

#   sleep 300

    cat <<EOF | sed -n '1,$s/^ *//p' > firex_catalog.j
    #!/bin/csh -fx

    #SBATCH --job-name=firex_catalog
    #SBATCH --account=s1321
    #SBATCH --time=0:40:00
    #SBATCH --qos=daohi
    #SBATCH --ntasks=28
    #SBATCH --export=NONE
    #SBATCH --constraint=hasw
    #SBATCH --output=/discover/nobackup/dao_ops/jardizzo/FLUID/firex_catalog_${idate}_$itime.log

    limit stacksize unlimited
    source /home/dao_ops/jardizzo/FLUID/firex-aq/utils/pyg_modules
    cd /home/dao_ops/jardizzo/FLUID/firex-aq/utils

    firex_catalog.sh $idate $itime
EOF

#   qsub firex_catalog.j
    exit 0

    dattim=`timetag $idate $itime {%Y%m%d_%H%M%S}%+H12`
    idate=`echo $dattim | cut -d'_' -f1`
    itime=`echo $dattim | cut -d'_' -f2`

  else

    sleep 300

  fi

done

exit 0
