#!/bin/csh

#############################################################################
# get_aero_obs.csh - retrieve observations for aerosol analysis
#
# !REVISION HISTORY:
#
#  23Jun2016  Todling   Created from Joe Stassi original GEOSdas.csm code
#
#############################################################################

setenv MYNAME get_aero_obs.csh

if ( $#argv < 3 ) then
   echo " "
   echo " \\begin{verbatim} "
   echo " "
   echo " NAME"
   echo " "
   echo "  $MYNAME  - retrieve observations for AEROSOL analysis "
   echo " "
   echo " SYNOPSIS"
   echo " "
   echo "  $MYNAME  nymd nhms nstep"
   echo " "
   echo " where"
   echo "   nymd   -  initial date of forecast, as in YYYYMMDD "
   echo "   nhms   -  initial time of forecast, as HHMMSS"
   echo "   nstep  -  number of times to  in bkg (im) (for history)"
   echo " "
   echo " DESCRIPTION "
   echo " "
   echo " "
   echo " Example of valid command line:"
   echo " $MYNAME  20151118 000000 1 "
   echo " "
   echo " REQUIRED RESOURCE FILES"
   echo " "
   echo "   obsys-gaas.rc    - database of aerosols observations"
   echo " "
   echo " REQUIRED ENVIRONMENT VARIABLES"
   echo " "
   echo "    ATMENSETC     - location of resource files        "
   echo "    ATMENSLOC     - location of current ensemble      "
   echo "    ASYNBKG       - frequency of background (minutes) "
   echo "    FVBCS         - location of fvInput               "
   echo "    FVHOME        - location of experiment            "
   echo "    FVROOT        - location of DAS build             "
   echo "    FVWORK        - location of work directory        "
   echo "    GID           - group ID to run job under         "
   echo "    TIMEINC       - analysis frequency (minutes)      "
   echo " "
   echo " OPTIONAL ENVIRONMENT VARIABLES"
   echo " "
   echo " "
   echo " REMARKS"
   echo " "
   echo " SEE ALSO"
   echo " "
   echo " AUTHOR"
   echo "   Ricardo Todling (Ricardo.Todling@nasa.gov), NASA/GMAO "
   echo "     Last modified: 24Jun2016      by: R. Todling"
   echo " \\end{verbatim} "
   echo " \\clearpage "
   exit(0)
endif
set bnymd  = $1
set bnhms  = $2
set nstep  = $3

set hh =  `echo $bnhms | cut -c1-2`
set yyyymmddhh = ${bnymd}${hh}


setenv FAILED 0

if ( !($?FVHOME)        ) setenv FAILED 1
if ( !($?FVROOT)        ) setenv FAILED 1
if ( !($?FVWORK)        ) setenv FAILED 1
if ( !($?GID)           ) setenv FAILED 1
if ( !($?group_list)    ) setenv FAILED 1

if ( -e $FVWORK/.DONE_${MYNAME}.$yyyymmddhh ) then
   echo "${MYNAME}: Already done."
   exit(0)
endif

if ( $FAILED ) then
  echo " ${MYNAME}: not all required env vars defined"
  exit 1
endif

if ( !($?AERO_OBSDBRC) )  setenv AERO_OBSDBRC  obsys-gaas.rc
if ( !($?DATA_QUEUE)   )  setenv DATA_QUEUE    datamove
if ( !($?DO_DMGET)     )  setenv DO_DMGET      1
if ( !($?IGNORE_0)     )  setenv IGNORE_0      0
if ( !($?MODIS_L2_HDF) )  setenv MODIS_L2_HDF  0

if ( $?FVSPOOL ) then
   set spool = "-s $FVSPOOL "
else
   set diren = `dirname $FVHOME`
   set spool = "-s $diren/spool "
endif

set path = ( . $FVHOME/run $FVROOT/bin $path )

  alias fname1 'echo \!* >! $fname'   # write first line of $fname
  alias fname2 'echo \!* >> $fname'   # append to $fname


             set fname = "acqaerobs1.pbs"
             if ( ! -e $fname ) /bin/rm $fname
             set acqdate = ${bnymd}_`echo $bnhms | cut -c1-2`z
             set acqlog = $FVHOME/run/acqaerobs1.log.$acqdate.txt
             touch ${FVWORK}/acquire.FAILED
             fname1 "#\!/bin/csh -xvf"
             fname2 "#$group_list"
             fname2 "#PBS -N acqaerobs1"
             fname2 "#PBS -o acqaerobs1.log.$acqdate.txt"
             fname2 "#PBS -l nodes=1:ppn=1"
             fname2 "#PBS -l walltime=2:00:00"
             fname2 "#PBS -q $DATA_QUEUE"
             fname2 "#PBS -S /bin/csh"
             fname2 "#PBS -j eo"
             fname2 ""
             fname2 "setenv IGNORE_0  $IGNORE_0"
             fname2 "setenv FVWORK $FVWORK"
             fname2 "setenv DO_DMGET $DO_DMGET"
             fname2 "set path = ( $path )"
             fname2 "cd $FVWORK"

                set flag = "-rc $AERO_OBSDBRC"
                @ numhrs = $nstep * 6

                # check for AVHRR data availability
                #----------------------------------
                obsys_check.pl $flag patmosx_asc $bnymd $bnhms $numhrs
                setenv PATMOSX $status

                # check whether user wants to access NRT MODIS HDF data 
                #------------------------------------------------------
                if ($MODIS_L2_HDF) then
                   setenv MOD04_NNR 0
                   setenv MYD04_NNR 0
                else

                   # check for MODIS Terra data availability
                   #----------------------------------------
                   obsys_check.pl $flag mod04_land_nnr $bnymd $bnhms $numhrs
                   setenv MOD04_NNR $status

                   # check for MODIS Aqua data availability
                   #---------------------------------------
                   obsys_check.pl $flag myd04_land_nnr $bnymd $bnhms $numhrs
                   setenv MYD04_NNR $status

                endif

                # check for misr data availability
                #---------------------------------
                obsys_check.pl $flag misr_F12_bright $bnymd $bnhms $numhrs
                setenv MISR_BRIGHT $status

                # check for AERONET data availability
                #------------------------------------
                obsys_check.pl $flag aeronet_obs $bnymd $bnhms $numhrs
                setenv AERONET $status

                # acquire available data
                #-----------------------
                if ($PATMOSX) then
                   @ mstep = $nstep  #$nstep * 2
                   fname2 "acquire_obsys -v -d $FVWORK $spool -ssh \"
                   fname2 "      -strict $bnymd $bnhms 030000 $mstep \"
                   fname2 "      patmosx_asc,patmosx_des -drc obsys-gaas.rc"
                   fname2 "setup_gaas_obs.pl $FVWORK -v -avhrr"
                endif
                if ($MOD04_NNR) then
                   @ mstep = $nstep #$nstep * 2
                   fname2 "acquire_obsys -v -d $FVWORK $spool -ssh \"
                   fname2 "      -strict $bnymd $bnhms 030000 $mstep \"
                   fname2 "      mod04_land_nnr,mod04_ocean_nnr -drc obsys-gaas.rc"
                endif
                if ($MYD04_NNR) then
                   @ mstep = $nstep #$nstep * 2
                   fname2 "acquire_obsys -v -d $FVWORK $spool -ssh \"
                   fname2 "      -strict $bnymd $bnhms 030000 $mstep \"
                   fname2 "      myd04_land_nnr,myd04_ocean_nnr -drc obsys-gaas.rc"
                endif
                #-----------------------------------------------------------
                # NOTE: ana_aod.j script accesses MODIS L2 HDF data directly
                #       from intermediate directory (see modis_l2.pcf)
                #-----------------------------------------------------------
                #--if ($MODIS_L2_HDF) then
                #--   set bdtgaas = (`tick $bnymd $bnhms -5400`)
                #--   @ mstep = $nstep * 72
                #--   fname2 "acquire_obsys -v -d $FVWORK $spool -ssh \"
                #--   fname2 "      $bdtgaas[1] $bdtgaas[2] 000500 $mstep \"
                #--   fname2 "      mod04_051_flk,myd04_051_flk"
                #--   fname2 "setup_gaas_obs.pl $FVWORK -v -modis"
                #--endif
                #-----------------------------------------------------------
                if ($MISR_BRIGHT) then
                   @ mstep = 1
                   fname2 "acquire_obsys -v -d $FVWORK $spool -ssh \"
                   fname2 "      -strict $bnymd $bnhms 240000 $mstep \"
                   fname2 "      misr_F12_bright -drc obsys-gaas.rc"
                endif
                if ($AERONET) then
                   @ mstep = 1
                   fname2 "acquire_obsys -v -d $FVWORK $spool -ssh \"
                   fname2 "      -strict $bnymd $bnhms 240000 $mstep \"
                   fname2 "      aeronet_obs -drc obsys-gaas.rc"
                endif

                fname2 ' if ( ! $status ) /bin/rm -f ${FVWORK}/acquire.FAILED'
                fname2 "exit"

                qsub -W block=true -o $acqlog $fname  # acquire observations; ignore status return

# if made it here, should be ok
# -----------------------------
touch $FVWORK/.DONE_${MYNAME}.$yyyymmddhh
echo " ${MYNAME}: Complete "
exit(0)
