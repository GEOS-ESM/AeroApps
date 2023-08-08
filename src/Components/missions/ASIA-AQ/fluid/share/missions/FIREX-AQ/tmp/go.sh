#!/bin/sh

while read field; do

  wxmap.py --theme custom_mission --theme FIREX-AQ --fcst_dt 20190715 --time_dt 20190715 --region fx-se --stream GEOS --level 0 --field $field --oname '$field.png'
  echo "Finished $field"
  sleep 10
done

exit 0
