#!/bin/sh

echo '<html>'
echo '<body bgcolor="#ffffff" text="black">'

#echo '<table style="width:100%">'
echo '<table>'

while read iname; do

  field=`echo $iname | cut -d'.' -f1`
  level=`echo $iname | cut -d'.' -f2`
  bname=`basename $iname`
  region=`basename $iname ".png" | tr '[a-z]' '[A-Z]'`
  thumb=`basename $iname ".png"`.thumb.png
# imtag='<img src="'$bname'" style="border:1px solid white;" />'
  imtag='<img src="'$bname'" />'

  echo '<tr>'

  echo '<td align="center" valign="center">'
  echo $imtag
  echo '</td>'

  echo '<td align="center" valign="center">'
  echo '<font size="+2">'$region'</font>'
  echo '</td>'

  echo '</tr>'

  echo '<tr> <td> <br> </td> </tr>'

done

echo '</body>'
echo '</html>'


exit 0
