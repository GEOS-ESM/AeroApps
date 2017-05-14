#!/bin/csh

#############################################################################
# run_gaas_ana.csh - retrieve observations for aerosol analysis
#
# !REVISION HISTORY:
#
#  23Jun2016  Todling   Created from Joe Stassi original GEOSdas.csm code
#  03Mar2017  Todling   Consistent w/ GEOSdas.csm in GEOSadas-5_16_5
#
#############################################################################

setenv MYNAME run_gaas_ana.csh

if ( $#argv < 6 ) then
   echo " "
   echo " \\begin{verbatim} "
   echo " "
   echo " NAME"
   echo " "
   echo "  $MYNAME  - run GAAS analysis"
   echo " "
   echo " SYNOPSIS"
   echo " "
   echo "  $MYNAME expid nymd nhms nstep"
   echo " "
   echo " where"
   echo "   expid   -  experiment name"
   echo "   nymd    -  initial date of forecast, as in YYYYMMDD "
   echo "   nhms    -  initial time of forecast, as HHMMSS"
   echo "   nstep   -  THIS WILL BE MADE OBSOLETE"
   echo "   workdir -  directory where work takes place"
   echo "   egressfn-  egress file name"
   echo " "
   echo " DESCRIPTION "
   echo " "
   echo " "
   echo " Example of valid command line:"
   echo " $MYNAME d512a 20151118 000000 1 workdir egressfn"
   echo " "
   echo " REQUIRED RESOURCE FILES"
   echo " "
   echo " "
   echo " REQUIRED ENVIRONMENT VARIABLES"
   echo " "
   echo "    ATMENSETC     - location of resource files        "
   echo "    ATMENSLOC     - location of current ensemble      "
   echo "    ASYNBKG       - frequency of background (minutes) "
   echo "    FVBCS         - location of fvInput               "
   echo "    FVHOME        - location of experiment            "
   echo "    FVROOT        - location of DAS build             "
   echo "    GID           - group ID to run job under         "
   echo "    TIMEINC       - analysis frequency (minutes)      "
   echo " "
   echo " OPTIONAL ENVIRONMENT VARIABLES"
   echo " "
   echo "   AODBLOCKJOB        - allows blocking batch jobs and halt "
   echo "                        before proceeding                   "
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
set expid = $1
set nymd  = $2
set nhms  = $3
set nstep = $4
set workdir  = $5
set egressfn = $6

set hh =  `echo $nhms | cut -c1-2`
set yyyymmddhh = ${nymd}${hh}

setenv FAILED 0

if ( !($?FVHOME)        ) setenv FAILED 1
if ( !($?FVROOT)        ) setenv FAILED 1
if ( !($?GID)           ) setenv FAILED 1
if ( !($?MPIRUN_AOD)    ) setenv FAILED 1

setenv FVWORK $workdir
if ( -e $FVWORK/.DONE_${MYNAME}.$yyyymmddhh ) then
   echo "${MYNAME}: Already done."
   exit(0)
endif

if ( $FAILED ) then
  echo " ${MYNAME}: not all required env vars defined"
  exit 1
endif

if ( !($?AODBLOCKJOB)  )  setenv AODBLOCKJOB   0
if ( !($?AERO_OBSDBRC) )  setenv AERO_OBSDBRC  obsys-gaas.rc
if ( !($?DATA_QUEUE)   )  setenv DATA_QUEUE    datamove
if ( !($?DO_DMGET)     )  setenv DO_DMGET      1
if ( !($?IGNORE_0)     )  setenv IGNORE_0      0
if ( !($?MODIS_L2_HDF) )  setenv MODIS_L2_HDF  0
if ( !($?NCPUS_AOD)    )  setenv NCPUS_AOD     1
if ( !($?PBS_BIN)      )  setenv PBS_BIN       "/usr/slurm/bin"
if ( !($?group_list)   )  setenv group_list    "SBATCH -A $GID"

if ( $?FVSPOOL ) then
   set spool = "-s $FVSPOOL "
else
   set diren = `dirname $FVHOME`
   set spool = "-s $diren/spool "
endif

set path = ( . $FVHOME/run $FVROOT/bin $path )

cd $FVWORK

  set aod_date = $nymd
  set aod_time = $nhms

  set flag = "-rc $AERO_OBSDBRC"
  @ numhrs = $nstep * 6

  # check for AVHRR data availability
  #----------------------------------
  obsys_check.pl $flag patmosx_asc $nymd $nhms $numhrs
  setenv PATMOSX $status

  # check whether user wants to access NRT MODIS HDF data 
  #------------------------------------------------------
  if ($MODIS_L2_HDF) then
     setenv MOD04_NNR 0
     setenv MYD04_NNR 0
  else

     # check for MODIS Terra data availability
     #----------------------------------------
     obsys_check.pl $flag mod04_land_nnr $nymd $nhms $numhrs
     setenv MOD04_NNR $status

     # check for MODIS Aqua data availability
     #---------------------------------------
     obsys_check.pl $flag myd04_land_nnr $nymd $nhms $numhrs
     setenv MYD04_NNR $status

  endif

  # check for misr data availability
  #---------------------------------
  obsys_check.pl $flag misr_F12_bright $nymd $nhms $numhrs
  setenv MISR_BRIGHT $status

  # check for AERONET data availability
  #------------------------------------
  obsys_check.pl $flag aeronet_obs $nymd $nhms $numhrs
  setenv AERONET $status


  # set up separate work directory
  #-------------------------------  
  setenv AODWORK $FVWORK/aod.$aod_date.$aod_time
  mkdir $AODWORK
  touch $AODWORK/.no_archiving

  # Copy over pcf files for NNR
  # ---------------------------
  /bin/cp $FVHOME/run/gaas/avhrr_l2a.pcf $AODWORK
  /bin/cp $FVHOME/run/gaas/modis_l2a.pcf $AODWORK
 
  # Prepare ana.rc from template
  # ----------------------------
  set ana_rc_tmpl = $FVHOME/run/gaas/ana.rc.tmpl
  set ana_rc = $AODWORK/ana.rc
  vED -env $ana_rc_tmpl -o $ana_rc

  if ($PATMOSX)      sed -i "s/#___AVHRR___//" $ana_rc
  if ($MODIS_L2_HDF) then
      sed -i "s/#___TERRA_NRT___//" $ana_rc
      sed -i "s/#___AQUA_NRT___//" $ana_rc
  endif
  if ($MOD04_NNR)    sed -i "s/#___TERRA___//" $ana_rc
  if ($MYD04_NNR)    sed -i "s/#___AQUA___//" $ana_rc
  if ($MISR_BRIGHT)  sed -i "s/#___MISR___//" $ana_rc
  if ($AERONET)      sed -i "s/#___AERONET___//" $ana_rc

  echo "cat $ana_rc"
  cat $ana_rc

  # copy other resource files and aerosol bkgs to work directory
  # ------------------------------------------------------------
  \cp $FVHOME/run/gaas/*.rc $AODWORK

  foreach aod_time_offset ( 000000 030000 )
     set gaastime = (`tick $aod_date $aod_time 0 $aod_time_offset`)
     set ymd = $gaastime[1]
     set hms = $gaastime[2]

     set flg1 = "-template $EXPID $ymd $hms"
     set flg2 = "-rc $AODWORK/gaas.rc gaas_bkg_filename"
     set aer_f = `echorc.x $flg1 $flg2`

     echo $aer_f
     if($aer_f == "") then  # this test is to catch funk behavior under MPT where sometime env not correct
        echo "file $aer_f not found in $AODWORK, aborting ... "
        exit 1
     endif
     /bin/cp $aer_f $AODWORK
  end

# write job script to AODWORK
# ---------------------------
  set tmpl = "$FVHOME/run/gaas/ana_aod.j.tmpl"
  set jobf = "$AODWORK/ana_aod.j"
  set aod_label = `echo ${aod_date}_${aod_time} | cut -c1-11`
  set gaasLOG  = $FVWORK/${EXPID}.ana_aod.log.${aod_label}z.txt
  set gaasLOGx = $FVHOME/run/ana_aod.abnormal.log.${aod_label}z.txt

  set siteID = `$FVROOT/bin/get_siteID.pl`
  set nodeflg = ""
  if ($siteID == "nas") then
     if ($NCPUS_PER_NODE == 12) set nodeflg = ":model=wes"
     @ num_aod_nodes = $NCPUS_AOD / $NCPUS_PER_NODE
     if ( $num_aod_nodes * $NCPUS_PER_NODE < $NCPUS_AOD ) @ num_aod_nodes++
     set pbsflg1 = "select=$num_aod_nodes"
     set pbsflg2 = "ncpus=$NCPUS_PER_NODE"
     set pbsflg3 = "mpiprocs=$NCPUS_PER_NODE"
     set line5 = "#PBS -l ${pbsflg1}:${pbsflg2}:$pbsflg3$nodeflg"
  else
     set line5 = "#PBS -l ncpus=$NCPUS_AOD"
  endif

  setenv GAASFAIL $AODWORK/ana_aod.FAILED.$$
  touch $GAASFAIL

  sed -e 5i"$line5" < $tmpl >! $jobf
  sed -i "s|>>>PATMOSX<<<|${PATMOSX}|" $jobf
  sed -i "s|>>>MODIS_L2_HDF<<<|${MODIS_L2_HDF}|" $jobf

  sed -i "s|>>>AODWORK<<<|${AODWORK}|" $jobf
  sed -i "s|>>>FVHOME<<<|${FVHOME}|" $jobf
  sed -i "s|>>>FVWORK<<<|${FVWORK}|" $jobf
  sed -i "s|>>>GAASFAIL<<<|${GAASFAIL}|" $jobf
  sed -i "s|>>>MPIRUN_AOD<<<|${MPIRUN_AOD}|" $jobf
  sed -i "s|>>>NCPUS_AOD<<<|${NCPUS_AOD}|" $jobf
  sed -i "s|>>>NYMD<<<|${aod_date}|" $jobf
  sed -i "s|>>>NHMS<<<|${aod_time}|" $jobf
 
  echo "cat $jobf"
  if ($status) then
     exit 1
  endif
  #cat $jobf

  setenv aod_parallel_flag 0 
  if ($aod_parallel_flag) then
     /bin/cp $FVWORK/AOD_list $AODWORK

     # run at command line
     #--------------------
     chmod 744 $jobf
     $jobf $aod_parallel_flag >&! $gaasLOG &

  else

     # submit job and save job ID
     #---------------------------
     if ( $AODBLOCKJOB ) then
        set jobID = `$PBS_BIN/qsub -W block=true -V -o $gaasLOG $jobf`
     else
        set jobID = `$PBS_BIN/qsub -V -o $gaasLOG $jobf`
     endif

     if ($?aodJobIDs) then
        setenv aodJobIDs "${aodJobIDs}:$jobID"
     else
        setenv aodJobIDs $jobID
     endif

  endif

# if made it here, all good
# -------------------------
if ( $egressfn == "DEFAULT" ) then
  echo " ${MYNAME}: Complete "
  touch $FVWORK/.DONE_${MYNAME}.$yyyymmddhh
else
  if ( -e $AODWORK/ANAAOD_EGRESS ) then
     echo " ${MYNAME}: Complete "
     touch $egressfn
  else
     exit (1)
  endif
endif
exit(0)

