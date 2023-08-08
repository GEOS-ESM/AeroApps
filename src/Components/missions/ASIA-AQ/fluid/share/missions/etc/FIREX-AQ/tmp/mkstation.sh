#!/bin/sh

#URL=https://fluid.nccs.nasa.gov/gram/total/LATxLON/?region=camp2ex
URL=https://fluid.nccs.nasa.gov/gram/static/plots/total_LAT_LON.png

while read station; do

  name=`echo $station | cut -d',' -f1 | tr '_' ' '`

  if [ "$name" == "name" ]; then continue; fi

# lon=`echo $station | cut -d',' -f2 | sed s'/0*$//'`
  lon=`echo $station | cut -d',' -f2`
  lat=`echo $station | cut -d',' -f3`

  ilon=`echo $lon | cut -d'.' -f1`
  if [ $ilon -lt 0 ]; then
    lon=`python -c "print $lon + 360.0"`
  fi

  url=`echo $URL | sed -n s/LAT/$lat/p | sed -n s/LON/$lon/p`
  echo "  - $lon $lat $name $url"

# ilon=`echo $lon | cut -d'.' -f1`
# if [ $ilon -lt 0 ]; then
#   lon=`python -c "print $lon + 360.0"`
#   echo "  - $lon $lat $name $url"
# elif [ $ilon -gt 180 ]; then
#   lon=`python -c "print $lon - 360.0"`
#   echo "  - $lon $lat $name $url"
# fi

done

exit 0
