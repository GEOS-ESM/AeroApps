#!/bin/sh

#URL=https://fluid.nccs.nasa.gov/gram/total/LATxLON/?region=camp2ex
URL=https://fluid.nccs.nasa.gov/gram/static/plots/total_LAT_LON.png

flag=$1

while read station; do

  name=`echo $station | cut -d',' -f1 | tr '_' ' '`

  if [ "$name" == "name" ]; then continue; fi

# lon=`echo $station | cut -d',' -f2 | sed s'/0*$//'`
  lon=`echo $station | cut -d',' -f2`
  ilon=`echo $lon | cut -d'.' -f1`
  lat=`echo $station | cut -d',' -f3 | tr ' ' '\0'`
  url=`echo $URL | sed -n s/LAT/$lat/p | sed -n s/LON/$lon/p`

  lon360=$lon
  if [ $ilon -lt 0 ]; then
    lon360=`python -c "print $lon + 360.0"`
  fi

  lon180=$lon
  if [ $ilon -gt 180 ]; then
    lon180=`python -c "print $lon - 360.0"`
  fi

  if [ $flag -eq 1 ]; then
      echo "  - $lon360 $lat $name $url"
  elif [ $flag -eq -1 ]; then
      echo "  - $lon180 $lat $name $url"
  elif [ $flag -eq 2 ]; then
      echo "  - $lon360 $lat $name $url"
      echo "  - $lon180 $lat $name $url"
  fi

done

exit 0
