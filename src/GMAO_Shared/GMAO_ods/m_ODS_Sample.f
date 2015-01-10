!====================================================================
!-----------------------------------------------------------------------
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!-----------------------------------------------------------------------
!BOP
! !MODULE: m_ODS_Sample -- Utility routines for the ODS sample program
!
! !DESCRIPTION:
!     This module contains the utility routines, global declarations
!     and hardwired parameters for the program ODS sample.
!
! !INTERFACE:
!
      module      m_ODS_Sample
      implicit    NONE
      private	! except
      public ::   Init
      public ::   Inquire_ODS
      public ::   Set_SynTimes
      public ::   NewTime
      public ::   Bin_Obs
      public ::   ItoA
!
! !REVISION HISTORY:
!     19Oct2000  C. Redder  Original code
!     15Jun2002  da Silva   Split from ods_sample.f to conform
!                           to module file name conventions, thus
!                           allowing automatic dependency generation.
!
!EOP
!-----------------------------------------------------------------

!     Special characters
!     ------------------
      character (  len  =   1 ), parameter ::
     .   BLK    = achar (  32 )             ! blank
      character (  len  =   1 ), parameter ::
     .   EOL    = achar (  10 )             ! end of line  mark

!     Error status and handling
!     -------------------------
      integer, parameter   ::
     .   NLines_Max   = 40
      character ( len = 100),        dimension ( NLines_Max ) ::
     .   ErrorMessage              ! Output message
      integer, parameter   ::
     .   MesBufSz     = 2000       ! Maximum length of ...
      character ( len = MesBufSz ) ::
     .   MessageBuffer             ! buffer for output message
      character ( len = 100 )  ::
     .   Fmt                       ! Format of output message
      character (len=*), parameter :: MyModule = 'ODS_Sample'

      contains
*...................................................................

!-------------------------------------------------------------------------
!         NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!-------------------------------------------------------------------------
!BOP
!
! !ROUTINE:  Init() --- Initializer
!
! !DESCRIPTION: This propgram reads command line arguments to come
!               up with several operating parameters. The command
!               line syntax is:
! \label{Sample:Init}
! \begin{verbatim}
!     ods_sample.x  [-o tplate] [-t f1:f2] nsyn infile(s)
! \end{verbatim}
!  {\bf Required parameters:}
! \begin{verbatim}
!        nsyn      -  number of synoptic times per day in the output
!                     file(s). Valid numbers are 1,2,3,4,6,8,12 and 24.
!        infile(s) -  a list input ODS files
! \end{verbatim}
!  {\bf Optional parameters:}
! \begin{verbatim}

!        -o tplate -  Specify GrADS-like template to tplate for output
!                     file name generation. The following descriptors
!                     are allowed: %s, %y4, %m2, %m3, %d2, and %h2. See
!                     the module, m_StrTemplate, in the mpeu library for
!                     a complete list of valid descriptors and their
!                     meanings. Default: %s.resampled.obs.%y4%m2%d2.ods.
!                     Note: the descriptor, %s, corresponds to the
!                     experiment ID as derived from the first input ODS
!                     file.
!                     Examples: The template,
!                          expid.preanal.obs.%y4%m2%d2.ods,
!
!                     will create daily output ODS files of the form,
!                          expid.preanal.obs.YYYYMMDD.ods
!
!                     and the template,
!                          expid.preanal.obs.%y4%m2%.ods
!
!                     will create monthly output ODS files of the form,
!                          expid.preanal.obs.YYYYMM.ods
!
!        -t f1:f2  -  Define the time window around the main synoptic
!                     time, t0, as [t0+f1*dt,t0+f2*dt) where dt =
!                     24/nsyn and f1 and f2 must be defined so that 0
!                     <= f1,f2 <= 1 and f1 + f2 = 1 (i.e. no gaps
!                     between windows). Default: f1 = f2 = 0.5.
!                     Examples:
!
!                     -t -0.5:+0.5  window centered around main
!                                   synoptic time
!                     -t  0.0:+1.0  forwardly  uncentered window
!                     -t -1.0:+0.0  backwardly uncentered window
!
! !INTERFACE: 
!
      subroutine Init ( expid, nsyn, f, tplate, nfiles, infiles ) 
!
! !OUTPUT PARAMETERS:
!
      use         m_die, only : Warn
      implicit    NONE
      character (len=*), intent (out) ::
     .   expid                 ! Experiment ID Tag
      integer,           intent (out) ::
     .   nsyn                  ! The number of synoptic times per day.
      real,              intent (out), dimension (2) ::
     .   f                     ! The time interval boundaries (t1,t2)
                               !    where t = 0Z.
      character (len=*), intent (out) ::
     .   tplate                ! The template for the output filenames.
      integer,           intent (out) ::
     .   nfiles                ! The number of input files.
      character (len=*), intent (out), dimension (:) ::
     .   infiles               ! The input filenames.
!
! !BUGS:
!
! !SEE ALSO: 
!
!     I2k_Label()   selects a label (key)
!
! !FILES USED:  
!
!     File name supplied on input. The file is opened, read and then
!     and then closed.
!
! !REVISION HISTORY: 
!
!     20Oct2000  C. Redder  Original code.
!
!EOP
!-------------------------------------------------------------------------

      character (len=*), parameter ::
     .           MyName = MyModule // '::Init'

!     Size of buffer and file names
!     -----------------------------
      integer     BufSize
      parameter ( BufSize = 255 )

!     Functions referenced for ...
!     ----------------------------
      integer     IArgC     ! ... extracting the number of command
                            !     arguments

*     Error handling
*     --------------
      character                * ( BufSize )
     .            ErrorMessage   ( 10 )

*     Extracting command-line input
*     -----------------------------
      integer     NArg, iArg, iArg_Basic
      integer     iArg_nsyn
      integer     iArg_TPlate
      integer     iArg_TRange
      character   Buffer       * ( BufSize )
      character   Buffer2      * ( BufSize )
      integer     len_buf

*     Setting parameters
*     ------------------
      integer     iFile, iColon, len_expid
      integer     nfiles_max
      logical     bad_f1f2, bad_nsyn
      real        tsum

*...........................................

*     Number of arguments in command-line input
*     -----------------------------------------
      NArg          = IArgC ()
      if ( NArg .eq. 0 ) call Usage ()

*     Initialize
*     ----------
      iArg_TPlate   =  0
      iArg_TRange   =  0
      iArg          =  0
      iArg_Basic    =  0
      nfiles        =  0
      nfiles_max    =  size ( infiles )

*     Scan for arguments
*     ------------------
      do while ( iArg .lt. NArg )
         iArg = iArg + 1
         call GetArg ( iArg, Buffer )

*        Option for ...
*        --------------

*        ... setting the template for the output filenames
*        -------------------------------------------------
         if      ( Buffer ( : 2 ) .eq. '-o' ) then
            iArg = iArg  + 1
            iArg_TPlate  = iArg

*        ... setting the time intervals
*        ------------------------------
         else if ( Buffer ( : 2 ) .eq. '-t' ) then
            iArg = iArg  + 1
            iArg_TRange  = iArg

*        Any other option produces an error
*        ----------------------------------
         else if ( Buffer ( : 1 ) .eq. '-'  ) then
            len_buf = Len_Trim ( Buffer )
            write ( ErrorMessage, 901 ) Buffer ( : len_buf )
            call Warn  ( MyName, ErrorMessage ( : 1 ) )
            call Usage ()

*        Non-optional arguments concering ...
*        ------------------------------------
         else
            iArg_Basic   =  iArg_Basic + 1

*           ... the number of the synoptic time per day
*           ---------------------
            if      ( iArg_Basic .eq. 1 ) then
               iArg_nsyn =  iArg

*           Any other argument is assumed to be an input filename
*           -----------------------------------------------------
            else
               nfiles    =  nfiles + 1
               if ( nfiles .le. nfiles_max )
     .            infiles ( nfiles ) = Buffer

            end if
         end if
      end do

*     Check to ensure that sufficient number of arguments exist
*     ---------------------------------------------------------
      if       ( nfiles .lt. 1 ) then
         write ( ErrorMessage, 902 )
         call Warn    ( MyName, ErrorMessage ( : 1 ) )
         call Usage   ()

*     Check to ensure that there are not too many input ODS files
*     -----------------------------------------------------------
      else if ( nfiles .gt. nfiles_max ) then
         write ( ErrorMessage, 903 ) nfiles, nfiles_max
         call Warn    ( MyName, ErrorMessage ( : 3 ) )
         call Usage   ()

      end if

*     Check to determine if the command line input ends
*     with an option expecting another argument
*     -------------------------------------------------
      if     ( max ( iArg_TPlate,
     .               iArg_TRange ) .gt. NArg ) then
         write ( ErrorMessage, 904 )
         call Warn    ( MyName, ErrorMessage ( : 2 ) )
         call Usage   ()

      end if

*     Extract the experiment ID tag ...
*     ---------------------------------
      len_expid = scan ( infiles ( 1 ), '.' ) - 1
c      len_expid = scan ( infiles ( 1 ), '._' ) - 1
      if ( len_expid .eq. 0 ) len_expid = Len_Trim ( infiles ( 1 ) )
      expid     = infiles ( 1 ) ( : len_expid )

*     ... the number of synoptic times per day
*     ----------------------------------------
      call GetArg ( iArg_nsyn, Buffer )
      len_buf  =  Len_Trim ( Buffer )
      bad_nsyn = .true.
      read ( Buffer, '(bn,i10)', err = 5 ) nsyn
      if ( nsyn             .le. 0 ) go to 5
      if ( mod ( 24, nsyn ) .ne. 0 ) go to 5  ! Need to look at this 
                                              !   more closely.
      bad_nsyn = .false.
 5    continue

      if ( bad_nsyn ) then
         write ( ErrorMessage, 905 ) Buffer ( : len_buf )
         call Warn    ( MyName, ErrorMessage ( : 1 ) )
         call Usage   ()

      end if

*     ... the output file tamplate
*     ----------------------------
      tplate = '%s.resampled.obs.%y4%m2%d2.ods'
      if ( iArg_TPlate .gt. 0 ) 
     .   call GetArg ( iArg_TPlate, tplate )

*     ... and the time range for each synoptic
*         time in the output files
*     ----------------------------------------
      Buffer  = '.5:.5'
      if ( iArg_TRange .gt. 0 ) 
     .   call GetArg ( iArg_TRange, buffer )
      len_buf = Len_Trim ( Buffer )
      iColon  = index ( Buffer, ':' )
      if ( iColon .eq. 0 ) iColon = len_buf + 1

      bad_f1f2 = .true.
      if ( iColon .le. 1       )    go to 10
      Buffer2   =  Buffer ( : iColon - 1   )
      read ( Buffer2, '(bn,f15.0)', err = 10 ) f ( 1 )
      if ( f ( 1 ) .lt. -1.0e-6 .or.
     .     f ( 1 ) .gt.  1.000001 ) go to 10

      if ( iColon .ge. len_buf )    go to 10
      Buffer2   =  Buffer (   iColon + 1 : )
      read ( Buffer2, '(bn,f15.0)', err = 10 ) f ( 2 )
      if ( f ( 2 ) .lt. -1.0e-6 .or.
     .     f ( 2 ) .gt.  1.000001 ) go to 10

      tsum = f ( 1 ) + f ( 2 )
      if ( tsum .gt. 1.000001 .or.
     .     tsum .lt. 0.999999 ) go to 10
      bad_f1f2 = .false.
 10   continue

      if ( bad_f1f2 ) then
         write ( ErrorMessage, 906 ) Buffer ( : len_buf )
         call Warn    ( MyName, ErrorMessage ( : 1 ) )
         call Usage   ()

      end if
      return
 901  format ( ' Illegal option, ', a )
 902  format ( ' Insufficient number of arguments. ' )
 903  format ( ' Too many input O-F file names. ',
     .         /, '   Number of filenames    = ', i5
     .         /, '   Maximum number allowed = ', i5 )
 904  format ( ' Command-line input ends with an ',   
     .         /, '   option expecting another argument. ' )
 905  format ( ' Invalid value for the number of synoptic times, ', a )
 906  format ( ' Invalid time window, ', a )

      end subroutine Init

*....................................................................

!-------------------------------------------------------------------------
!         Nasa/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!-------------------------------------------------------------------------
! BOP
!
! !ROUTINE:  Usage () --- Writes usage
! 
! !DESCRIPTION: This routine prints out a summary of the
!               interface for command-line input. The routine
!               prints to standard error.
!
! !INTERFACE:
!
      subroutine Usage ()
!
! !INPUT PARAMETERS:
!
      use          m_die, only : die
      implicit     NONE
!
! !SEE ALSO:
!
! !REVISION HISTORY: 
!
!  12Oct2000  C. Redder   Origional code
!
! EOP
!-------------------------------------------------------------------------

      character (len=*), parameter ::
     .           MyName = 'Usage:'
      integer     iLine, NLines

      iLine  = 0
      iLine  = min ( iLine + 1, NLines_Max )
      ErrorMessage ( iLine )
     .       =  ' '
      iLine  = min ( iLine + 1, NLines_Max )
      ErrorMessage ( iLine )
     .       =  ' ods_sample.x  [-o tplate] [-t f1:f2] '
     .       //                 'nsyn infile(s)'
      iLine  = min ( iLine + 1, NLines_Max )
      ErrorMessage ( iLine )
     .       =  ' '
      iLine  = min ( iLine + 1, NLines_Max )
      ErrorMessage ( iLine )
     .       =  ' Required parameters: '
      iLine  = min ( iLine + 1, NLines_Max )
      ErrorMessage ( iLine )
     .       =  '    nsyn      -  number of synoptic times per '
     .       //                  'day in the output '
      iLine  = min ( iLine + 1, NLines_Max )
      ErrorMessage ( iLine )
     .       =  '                 file(s). Valid numbers are 1,2,3,4,6,'
     .       //                  '8,12 and 24.'
      iLine  = min ( iLine + 1, NLines_Max )
      ErrorMessage ( iLine )
     .       =  '    infile(s) -  a list input ODS files '
      iLine  = min ( iLine + 1, NLines_Max )
      ErrorMessage ( iLine )
     .       =  ' '
      iLine  = min ( iLine + 1, NLines_Max )
      ErrorMessage ( iLine )
     .       =  ' Optional parameters: '
      iLine  = min ( iLine + 1, NLines_Max )
      ErrorMessage ( iLine )
     .       =  '    -o tplate -  Specify GrADS-like template to '
     .       //                  'tplate for output'
      iLine  = min ( iLine + 1, NLines_Max )
      ErrorMessage ( iLine )
     .       =  '                 file name generation. The following '
     .       //                  'descriptors'
      iLine  = min ( iLine + 1, NLines_Max )
      ErrorMessage ( iLine )
     .       =  '                 are allowed: %s, %y4, %m2, %m3, %d2, '
     .       //                  'and %h2. See'
      iLine  = min ( iLine + 1, NLines_Max )
      ErrorMessage ( iLine )
     .       =  '                 the module, m_StrTemplate, in the '
     .       //                  'mpeu library for'
      iLine  = min ( iLine + 1, NLines_Max )
      ErrorMessage ( iLine )
     .       =  '                 a complete list of valid descriptors '
     .       //                  'and their'
      iLine  = min ( iLine + 1, NLines_Max )
      ErrorMessage ( iLine )
     .       =  '                 meanings. Default: '
     .       //                  '%s.resampled.obs.%y4%m2%d2.ods.'
      iLine  = min ( iLine + 1, NLines_Max )
      ErrorMessage ( iLine )
     .       =  '                 Note: the descriptor, %s, '
     .       //                  'corresponds to the '
      iLine  = min ( iLine + 1, NLines_Max )
      ErrorMessage ( iLine )
     .       =  '                 experiment ID as derived from the '
     .       //                  'first input ODS '
      iLine  = min ( iLine + 1, NLines_Max )
      ErrorMessage ( iLine )
     .       =  '                 file.'
      iLine  = min ( iLine + 1, NLines_Max )
      ErrorMessage ( iLine )
     .       =  '    -t f1:f2  -  Define the time window around the '
     .       //                  'main synoptic'
      iLine  = min ( iLine + 1, NLines_Max )
      ErrorMessage ( iLine )
     .       =  '                 time, t0, as [t0+f1*dt,t0+f2*dt) '
     .       //                  'where dt = '
      iLine  = min ( iLine + 1, NLines_Max )
      ErrorMessage ( iLine )
     .       =  '                 24/nsyn and f1 and f2 must be '
     .       //                  'defined so that 0 '
      iLine  = min ( iLine + 1, NLines_Max )
      ErrorMessage ( iLine )
     .       =  '                 <= f1,f2 <= 1 and f1 + f2 = 1 (i.e. '
     .       //                  'no gaps'
      iLine  = min ( iLine + 1, NLines_Max )
      ErrorMessage ( iLine )
     .       =  '                 between windows). Default: f1 = '
     .       //                  'f2 = 0.5. '

*     Write message
*     -------------
      NLines = iLine
      call die ( MyName, ErrorMessage ( : NLines ) )

      return
      end subroutine Usage 

*....................................................................
!-------------------------------------------------------------------------
!         NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!-------------------------------------------------------------------------
! BOP
!
! !ROUTINE:  Info() --- Print message 
! 
! !DESCRIPTION: Prints a message to stdout.
!
! !INTERFACE:
!
      subroutine Info ( Message )
!
! !INPUT PARAMETERS: 
!
      use        m_stdio, only : stdout
      implicit   NONE
      character ( len = * ), intent (in), dimension (:) ::  Message

!
! !REVISION HISTORY: 
!
!     12Oct2000  C. Redder   Initial code
!
! EOP
!-------------------------------------------------------------------------

      integer :: LLine, iLine, LMes, LMes_New, NLines

*     For each line ...
*     -----------------
      NLines         = size ( Message )
      LMes           = 0
      do iLine = 1, NLines
         LLine       = Len_Trim ( Message ( iLine ) )

*        ... append to buffer if the line is ...
*        ---------------------------------------
         if ( LLine .eq. 0 ) then
            LMes_New = min ( LMes + 2,         MesBufSz )
            MessageBuffer ( LMes + 1 : LMes_New )
     .               = BLK // EOL                 ! ... blank
                                                  ! -------------
         else
            LMes_New = min ( LMes + LLine + 1, MesBufSz )
            MessageBuffer ( LMes + 1 : LMes_New )
     .               = Message ( iLine ) ( : LLine ) // EOL
     .                                            ! ... non-blank
                                                  ! -------------
         end if
         LMes = LMes_New

      end do
      write ( stdout, '( a )' ) MessageBuffer ( : LMes )

      return
      end subroutine Info

*....................................................................
!-------------------------------------------------------------------------
!         NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!-------------------------------------------------------------------------
! BOP
!
! !ROUTINE:  ItoA() --- Convert integer to string with left justification
! 
! !DESCRIPTION: Convert integer to string with left justification
!
! !INTERFACE:
!
      character (len=20) function ItoA ( Num )
!
! !INPUT PARAMETERS: 
!
      implicit NONE
      integer, intent (in)  ::  Num
!
! !REVISION HISTORY: 
!
!     24Oct2000  C. Redder   Initial code
!
! EOP
!-------------------------------------------------------------------------

      integer                :: iNum
      character ( len = 20 ) :: str

      write ( str, '(i20)' ) Num
      iNum = verify   ( str, BLK )
      ItoA = str      ( iNum : )

      return
      end function ItoA

!.................................................................

!-------------------------------------------------------------------------
!         Nasa/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!-------------------------------------------------------------------------
! BOP
!
! !ROUTINE:  Inquire_ODS () --- Gets ODS file parameters
! 
! !DESCRIPTION: The routine gets first and last valid synoptic
!               Julian hour as well as the number of synoptic 
!               hours on file.
!
! !INTERFACE:
!
      subroutine Inquire_ODS ( File, NSyn, JHrBeg, JHrEnd )
!
! !INPUT PARAMETERS:
!
      use        m_die, only : die
      implicit   NONE
      character (len=*), intent (in)  ::
     .   File            ! Name of ODS file
!
! !OUTPUT PARAMETERS:
!
      integer,           intent (out) ::
     .   NSyn            ! Number of synoptic times per day
      integer,           intent (out) ::
     .   JHrBeg          ! Earliest ...
      integer,           intent (out) ::
     .   JHrEnd          ! ... latest synoptic Julian hour on file
!
! !SEE ALSO:
!
! !REVISION HISTORY: 
!
!     19Oct2000  C. Redder   Origional code
!
! EOP
!-------------------------------------------------------------------------
!
!     Find beginning and ending Julian days on file.
!
      integer id, ier, JDayBeg, JDayEnd, HrEnd
      character (len=*), parameter :: MyName =   MyModule // '::'
     .                                       // 'Inquire_ODS'

*     Open the input file ...
*     -----------------------
      call ODS_Open ( id, trim ( File ), 'r', ier ) ! open the file
      if ( ier .ne. 0 ) call die ( MyName, 'Could not open ods file '
     .                               // trim ( File ) )

*     ... and extract the number of synoptic times per day
*     ----------------------------------------------------
      call ODS_IGet ( id, 'nsyn', nsyn, ier )
      if ( ier .ne. 0 ) call die ( MyName, 'could not read ods file '
     .                               // trim ( File ) )

*     ... the first Julian day on file
*     --------------------------------
      call ODS_IGet ( id, 'syn_beg:first_julian_day',
     .                JDayBeg, ier )
      if ( ier .ne. 0 ) call die ( MyName, 'could not read ods file '
     .                               // trim ( File ) )

*     ... the last Julian day and ...
*     --------------------------------
      call ODS_IGet ( id, 'syn_beg:latest_julian_day',
     .                JDayEnd, ier )
      if ( ier .ne. 0 ) call die ( MyName, 'could not read ods file '
     .                               // trim ( File ) )
*     ... synoptic hour on file.
*     --------------------------
      call ODS_IGet ( id, 'syn_beg:latest_synoptic_hour',
     .                HrEnd, ier )
      if ( ier .ne. 0 ) call die ( MyName, 'could not read ods file '
     .                               // trim ( File ) )

*     Determine the first and last Julian hour on file
*     ------------------------------------------------
      JHrBeg = ( JDayBeg - 1 ) * 24
      JHrEnd = ( JDayEnd - 1 ) * 24 + HrEnd

*     Close input file
*     ----------------
      call ODS_close ( id, MyName, ier )

      end Subroutine Inquire_ODS

!.................................................................

!-------------------------------------------------------------------------
!         Nasa/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!-------------------------------------------------------------------------
! BOP
!
! !ROUTINE:  Set_SynTimes () --- Sets synoptic times for output file
! 
! !DESCRIPTION: The routine sets the first and last synoptic times
!               for the output file for a given synoptic time for
!               the input file.
!
! !INTERFACE:
!
      subroutine Set_SynTimes ( JHr1,    SamTimes, f, NSyn2, 
     .                          JHrBeg2, JHrEnd2 )
!
! !INPUT PARAMETERS:
!
      implicit   NONE
      integer,           intent (in)  ::
     .   JHr1            ! Synoptic Julian hour assigned during data
                         !   aquisition from the input file.
      integer,           intent (in), dimension (:) ::
     .   SamTimes        ! Sample times of each observation (in minutes
                         !   since 0:00 GMT) 
      integer,           intent (in)  ::
     .   NSyn2           ! Number of synoptic times per day in the
                         !   output file.
      real,              intent (in), dimension (2) ::
     .   f               ! Factors for defining the boundaries of the
                         !   time window for each synoptic time in the
                         !   output file
!
! !OUTPUT PARAMETERS:
!
      integer,           intent (out) ::
     .   JHrBeg2         ! Earliest ...
      integer,           intent (out) ::
     .   JHrEnd2         ! ... latest synoptic Julian hour
!
! !SEE ALSO:
!
! !REVISION HISTORY: 
!
!     19Oct2000  C. Redder   Origional code
!
! EOP
!-------------------------------------------------------------------------

      integer :: SynTime1, TimeMin1, TimeMax1
      integer ::           TimeMin2, TimeMax2

*     Set time window for input file
*     ------------------------------
      SynTime1 = mod    ( JHr1, 24 ) * 60
      TimeMin1 = minval ( SamTimes ) + SynTime1
      TimeMax1 = maxval ( SamTimes ) + SynTime1

*     Set first ...
*     -------------
      JHrBeg2  = JHr1 - mod ( JHr1, 24 )
      TimeMin2 = nint ( -f ( 1 ) * real ( 60 * 24 / NSyn2 ))
      do while ( TimeMin2 .gt. TimeMin1 )
         TimeMin2 = TimeMin2 - 60 * 24 / NSyn2
         JHrBeg2  = JHrBeg2  -      24 / NSyn2

      end do

      TimeMin2 = TimeMin2 + 60 * 24 / NSyn2
      do while ( TimeMin2 .lt. TimeMin1 )
         TimeMin2 = TimeMin2 + 60 * 24 / NSyn2
         JHrBeg2  = JHrBeg2  +      24 / NSyn2

      end do

*     last synoptic times for the output file.
*     ----------------------------------------
      JHrEnd2  = JHr1 - mod ( JHr1, 24 )
      TimeMax2 = nint ( f ( 2 ) * real ( 60 * 24 / NSyn2 ))
      do while ( TimeMax2 .le. TimeMax1 )
         TimeMax2 = TimeMax2 + 60 * 24 / NSyn2
         JHrEnd2  = JHrEnd2  +      24 / NSyn2

      end do

      TimeMax2 = TimeMax2 - 60 * 24 / NSyn2
      do while ( TimeMax2 .gt. TimeMax1 )
         TimeMax2 = TimeMax2 - 60 * 24 / NSyn2
         JHrEnd2  = JHrEnd2  -      24 / NSyn2

      end do

      return
      end subroutine Set_SynTimes

!..............................................................


!-------------------------------------------------------------------------
!         NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!-------------------------------------------------------------------------
! BOP
!
! !ROUTINE:      NewTime - Adjust sampling time for new reference Julian hour
!
! !INTERFACE:
      subroutine NewTime ( RefJHr1, RefJHr2, Time )
!
! !INPUT PARAMETERS:
      implicit NONE
      integer, intent (in) ::  RefJHr1  ! Old and ...
      integer, intent (in) ::  RefJHr2  ! ... new reference Julian hour
!
! !OUTPUT PARAMETERS:
      integer, intent (inout), dimension (:) ::
     .                         Time     ! Time (in minutes) since the
                                        !   reference Julian hour
!
! !DESCRIPTION:
! \label{MODS:NewRefJHr}
!     Adjusts the time to account for the change in the reference
!     Julian hour.
!
! !REVISION HISTORY: 
!     10Nov2000   Redder   Origional version, adapted from the obs_io
!                          library
!
! EOP
!-------------------------------------------------------------------------

      integer  iVal, NVal, AdjTime

      AdjTime = 60 * ( RefJHr1 - RefJHr2 )

      NVal    = size ( Time )
      do iVal = 1, NVal
         Time ( iVal ) = Time ( iVal ) + AdjTime

      end do

      return
      end subroutine NewTime

!.................................................................

!-------------------------------------------------------------------------
!         Nasa/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!-------------------------------------------------------------------------
! BOP
!
! !ROUTINE:  Bin_Obs () --- Bin observations according synoptic time window.
! 
! !DESCRIPTION: The routine bins the observations according to the
!               synoptic time windows as defined by command-line input.
!
! !INTERFACE:
!
      subroutine Bin_Obs ( JHr1, f,   NSyn2, JHrBeg2,
     .                     Type, Obs, Binned_Obs )
!
! !INPUT PARAMETERS:
!
      use m_ods, only : ods_vect
      use m_die, only : die

      implicit   NONE
      integer,           intent (in) ::
     .   JHr1            ! Synoptic Julian hour assigned during data
                         !   aquisition from the input file.
      real,              intent (in), dimension (2) ::
     .   f               ! Factors for defining the boundaries of the
                         !   time window for each synoptic time in the
                         !   output file
      integer,           intent (in) ::
     .   NSyn2           ! Number of synoptic times per day in the
                         !   output file.
      integer,           intent (in) ::
     .   JHrBeg2         ! The Julian hour of the first synoptic time
                         !   window.
      character (len=*), intent (in) ::
     .   Type            ! ODS file type: 'pre_anal' or 'post_anal'
!
! !INPUT/OUTPUT PARAMETERS:
!
      type ( ods_vect ), intent (inout) ::
     .   Obs             ! Observation vectors obtained from the input
                         !   file.  On output, the observations are
                         !   sorted according to the assigned synoptic
                         !   times for the output files.
      type ( ods_vect ), intent (inout), dimension (:) ::
     .   Binned_Obs      ! Binned observations.  The number of synoptic
                         !   time windows is determined by the array
                         !   size.
!
! !SEE ALSO:
!
! !REVISION HISTORY: 
!
!     20Oct2000  C. Redder   Origional code
!     14Jun2002  Todling     Added Xvec
!
! EOP
!-------------------------------------------------------------------------

      character (len=*), parameter ::
     .           MyName = MyModule // '::Init'

      integer :: iOb, NObs
      integer :: NAllSyn, iSyn
      integer :: iSegBeg, iSegEnd, LSeg
      integer :: JHrEnd2
      integer, dimension (:), allocatable ::
     .           SegLoc, SegLen, SynNum
      integer :: TimeMin2, Time2
      integer :: iret

      NAllSyn = size ( Binned_Obs )
      JHrEnd2 = JHrBeg2 + ( NAllSyn - 1 ) * NSyn2 / 24

!     Allocate space for the segment locations ...
!     --------------------------------------------
      allocate ( SegLoc ( NAllSyn ), stat = iret )
      if ( iret .ne. 0 ) then
         write ( ErrorMessage, 901 ) 'SegLoc', 
     .                                trim ( ItoA ( NAllSyn )),
     .                                trim ( ItoA ( iret    ))
         call die ( MyName, ErrorMessage ( :1 ) )

      end if

!     ... and location of each bin
!     ----------------------------
      allocate ( SegLen ( NAllSyn ), stat = iret )
      if ( iret .ne. 0 ) then
         write ( ErrorMessage, 901 ) 'SegLen', 
     .                                trim ( ItoA ( NAllSyn )),
     .                                trim ( ItoA ( iret    ))
         call die ( MyName, ErrorMessage ( :1 ) )

      end if
      SegLen = 0

!     ... synoptic time index for all observations
!     --------------------------------------------
      NObs = Obs % Data % NObs
      allocate ( SynNum ( NObs ), stat = iret )
      if ( iret .ne. 0 ) then
         write ( ErrorMessage, 901 ) 'SynLen', 
     .                                trim ( ItoA ( NObs    )),
     .                                trim ( ItoA ( iret    ))
         call die ( MyName, ErrorMessage ( :1 ) )

      end if

!     Determine the synoptic time index for each observation
!     ------------------------------------------------------
      TimeMin2 = nint ( -f ( 1 ) * real ( 60 * 24 / NSyn2 ))
     .         - ( JHr1 - JHrBeg2 ) * 60
      do iOb = 1, NObs
         Time2          =   Obs % Data % Time ( iOb )
         SynNum ( iOb ) = ( Time2 - TimeMin2 )
     .                  / ( 60 * 24 / NSyn2 ) + 1

!        ... and the size and ...
!        ------------------------
         SegLen ( SynNum ( iOb ) ) = SegLen ( SynNum ( iOb ) ) + 1

      end do

!     ... and the location of each segment
!     ------------------------------------
      iSegBeg = 1
      do iSyn = 1, NAllSyn
         SegLoc ( iSyn ) = iSegBeg
         iSegBeg         = iSegBeg + SegLen ( iSyn )

      end do

!     Sort the observations according to the bin numbers
!     --------------------------------------------------
      do iOb = 1, NObs
         iSyn = SynNum ( iOb  )
         Obs % Data % kid ( SegLoc ( iSyn ) ) = iOb
         SegLoc ( iSyn ) = SegLoc ( iSyn ) + 1

      end do

      Obs % Data % Lat      ( : NObs )
     .    = Obs % Data % Lat    ( Obs % Data % kid ( : NObs ) )
      Obs % Data % Lon      ( : NObs )
     .    = Obs % Data % Lon    ( Obs % Data % kid ( : NObs ) )
      Obs % Data % Lev      ( : NObs )
     .    = Obs % Data % Lev    ( Obs % Data % kid ( : NObs ) )
      Obs % Data % kx       ( : NObs )
     .    = Obs % Data % kx     ( Obs % Data % kid ( : NObs ) )
      Obs % Data % kt       ( : NObs )
     .    = Obs % Data % kt     ( Obs % Data % kid ( : NObs ) )
      Obs % Data % ks       ( : NObs )
     .    = Obs % Data % ks     ( Obs % Data % kid ( : NObs ) )
      Obs % Data % xm       ( : NObs )
     .    = Obs % Data % xm     ( Obs % Data % kid ( : NObs ) )
      Obs % Data % Time     ( : NObs )
     .    = Obs % Data % Time   ( Obs % Data % kid ( : NObs ) )
      Obs % Data % Obs      ( : NObs )
     .    = Obs % Data % Obs    ( Obs % Data % kid ( : NObs ) )
      Obs % Data % qcExcl   ( : NObs )
     .    = Obs % Data % qcExcl ( Obs % Data % kid ( : NObs ) )
      Obs % Data % qcHist   ( : NObs )
     .    = Obs % Data % qcHist ( Obs % Data % kid ( : NObs ) )


      if ( type .eq. 'post_analysis' ) then
         Obs % Data % OmF   ( : NObs )
     .       = Obs % Data % OmF    ( Obs % Data % kid ( : NObs ) )
         Obs % Data % OmA   ( : NObs )
     .       = Obs % Data % OmA    ( Obs % Data % kid ( : NObs ) )
         Obs % Data % Xvec  ( : NObs )
     .       = Obs % Data % Xvec   ( Obs % Data % kid ( : NObs ) )

      end if

!     Reset the segment locations
!     ---------------------------
      iSegBeg = 1
      do iSyn = 1, NAllSyn
         SegLoc ( iSyn ) = iSegBeg
         iSegBeg         = iSegBeg + SegLen ( iSyn )

      end do

!     For each synoptic time, define ...
!     ----------------------------------
      do iSyn = 1, NAllSyn
         LSeg    =  SegLen ( iSyn )
         iSegBeg =  SegLoc ( iSyn )
         iSegEnd = iSegBeg + LSeg - 1

!        ... the size and
!        ----------------
         Binned_Obs     ( iSyn ) % Meta % nkt
     .         =             Obs % Meta % nkt
         Binned_Obs     ( iSyn ) % Meta % nkx
     .         =             Obs % Meta % nkx
         Binned_Obs     ( iSyn ) % Meta % nqc
     .         =             Obs % Meta % nqc
         Binned_Obs     ( iSyn ) % Meta % ncr
     .         =             Obs % Meta % ncr
         Binned_Obs     ( iSyn ) % Meta % nsyn
     .         =             NSyn2
!ams         Binned_Obs     ( iSyn ) % Data % nsyn
!ams     .         =             NSyn2

!        ... contents of the meta data structure
!        ---------------------------------------
         Binned_Obs ( iSyn ) % Meta % kt_names
     .         =>        Obs % Meta % kt_names
         Binned_Obs ( iSyn ) % Meta % kt_units
     .         =>        Obs % Meta % kt_units
         Binned_Obs ( iSyn ) % Meta % kx_names
     .         =>        Obs % Meta % kx_names
         Binned_Obs ( iSyn ) % Meta % kx_meta
     .         =>        Obs % Meta % kx_meta
         Binned_Obs ( iSyn ) % Meta % qcx_names
     .         =>        Obs % Meta % qcx_names

!        ... the number of 
!        ----------------
         Binned_Obs ( iSyn ) % Data % NObs = LSeg
         Binned_Obs ( iSyn ) % Data % NVct = LSeg

!        ... and the content of the observation vectors
!        ----------------------------------------------
         Binned_Obs ( iSyn ) % Data % kid
     .         => Obs % Data % kid    ( iSegBeg : iSegEnd )
         Binned_Obs ( iSyn ) % Data % Lat
     .         => Obs % Data % Lat    ( iSegBeg : iSegEnd )
         Binned_Obs ( iSyn ) % Data % Lon
     .         => Obs % Data % Lon    ( iSegBeg : iSegEnd )
         Binned_Obs ( iSyn ) % Data % Lev
     .         => Obs % Data % Lev    ( iSegBeg : iSegEnd )
         Binned_Obs ( iSyn ) % Data % kx 
     .         => Obs % Data % kx     ( iSegBeg : iSegEnd )
         Binned_Obs ( iSyn ) % Data % kt 
     .         => Obs % Data % kt     ( iSegBeg : iSegEnd )
         Binned_Obs ( iSyn ) % Data % ks
     .         => Obs % Data % ks     ( iSegBeg : iSegEnd )
         Binned_Obs ( iSyn ) % Data % xm
     .         => Obs % Data % xm     ( iSegBeg : iSegEnd )
         Binned_Obs ( iSyn ) % Data % Time
     .         => Obs % Data % Time   ( iSegBeg : iSegEnd )
         Binned_Obs ( iSyn ) % Data % Obs
     .         => Obs % Data % Obs    ( iSegBeg : iSegEnd )
         Binned_Obs ( iSyn ) % Data % qcExcl
     .         => Obs % Data % qcExcl ( iSegBeg : iSegEnd )
         Binned_Obs ( iSyn ) % Data % qcHist
     .         => Obs % Data % qcHist ( iSegBeg : iSegEnd )


         if ( type .eq. 'post_analysis' ) then
            Binned_Obs ( iSyn ) % Data % OmF
     .         => Obs % Data % OmF    ( iSegBeg : iSegEnd )
            Binned_Obs ( iSyn ) % Data % OmA
     .         => Obs % Data % OmA    ( iSegBeg : iSegEnd )
            Binned_Obs ( iSyn ) % Data % Xvec
     .         => Obs % Data % Xvec   ( iSegBeg : iSegEnd )

         end if

         Binned_Obs ( iSyn ) % Data % kid ( : LSeg )
     .         =  ( / ( iOb, iOb = 1, LSeg ) / )

      end do

!     Clean up
!     --------
      deallocate ( SegLoc )
      deallocate ( SegLen )
      deallocate ( SynNum )

      return
 901  format ( ' Error in allocating the array ', a, '(', a,
     .          '); stat = ', a )

      end subroutine Bin_Obs
!....................................................................

      end module m_ODS_Sample
