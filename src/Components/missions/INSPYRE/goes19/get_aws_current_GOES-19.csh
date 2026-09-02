#!/bin/csh 

    set USER = "Govindaraj5195"
    set TAG  = "Ushag%25"
    set fixgrp = 'chmod -R o+r,o+X,g+r,g+w,g+X \!*; chgrp -R g5extdata \!*'
    set DAY_TABLE = ( 31 29 31 30 31 30 31 31 30 31 30 31 )
    set INST_DIR = GOES-19
    set HOST_DIR = "noaa-goes19"
#    set STORAGE_DIR = /css/geostationary/NonOptimized/L1/${INST_DIR}-ABI-L1B-FULLD
    set STORAGE_DIR = ./${INST_DIR}-ABI-L1B-FULLD
#   set STORAGE_DIR = /discover/nobackup/rgovinda/NonOptimized/L1/${INST_DIR}-ABI-L1B-FULLD
#   set STORAGE_DIR = /css/geostationary/BackStage/${INST_DIR}-ABI-L1B-FULLD
#   set STORAGE_DIR = /discover/nobackup/projects/eis_fire/data/${INST_DIR}-ABI-L1B-FULLD

    set MONTH_TABLE = ( 01 02 03 04 05 06 07 08 09 10 11 12 )

      set nymd_beg = `date '+%Y%m%d'`
      set YYYY = `echo $nymd_beg | cut -c1-4`
      set   MM = `echo $nymd_beg | cut -c5-6`
      set   DD = `echo $nymd_beg | cut -c7-8`

      set JUL_DAYs = `./get_julday $nymd_beg`
      set JDAY = $JUL_DAYs[2]
      echo $JDAY

#          set SYNOP_TABLE =  ( 00 01 02 03 04 05 06 07 08 09 10 11 12 13 14 15 16 17 18 19 20 21 22 23 )

           set SYNOP_TABLE =  `./get_available_GOES-19_synoptic $HOST_DIR $YYYY $JDAY`
           foreach SYNOP ( `echo $SYNOP_TABLE` )
             mkdir -p $STORAGE_DIR/$YYYY/$JDAY/$SYNOP
             cd $STORAGE_DIR/$YYYY/$JDAY/$SYNOP
            nohup aws s3 ls --no-sign-request s3://$HOST_DIR/ABI-L1b-RadF/$YYYY/$JDAY/$SYNOP/ >& file.list

             set kount = 0
             set lknt  = 0
             set fknt = `cat file.list | awk '{print $4}' | wc -l`
             foreach FILE ( `cat file.list | awk '{print $4}'` )
             @ kount = $kount + 1
             @ lknt = $lknt + 1
               nohup  aws s3 cp --no-sign-request s3://$HOST_DIR/ABI-L1b-RadF/$YYYY/$JDAY/$SYNOP/$FILE . &
 
                 if ( $kount == 4 || $lknt == $fknt ) then
                   wait
                   set kount = 0
                endif
            end


           /bin/rm file.list
         end
         $fixgrp $STORAGE_DIR/$YYYY/$JDAY
         chgrp g5extdata $STORAGE_DIR/$YYYY/$JDAY
         @ DD = $DD + 1
