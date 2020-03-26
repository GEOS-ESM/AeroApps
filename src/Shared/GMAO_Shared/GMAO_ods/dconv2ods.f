!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!-----------------------------------------------------------------------
!BOP
!
! !ROUTINE: dconv2ods:  converts diag_conv files to ODS files
!
! !INTERFACE:

      program dconv2ods

! !USES:

      use m_ods

      implicit none 
 
! !USAGE: see the routine usage() below
!
! !DESCRIPTION: Converts GSI diag\_conv-type files to ODS.
!               This is a preliminary version and not all
!               entries in the diag\_conv file will appear in
!               the corresponding ODS file.
!
! !REVISION HISTORY:
!
!	29Apr2004  Todling  - Initial code.
!
!EOP
      integer, parameter :: NFILES_MAX = 100 ! Max number of input files

      character(len=255) :: fname

      integer  ios, ier, lu, i, k, ik, id, is(1), iks, ikx, ns
      integer  iu, iv, iknow, nfiles, ifile
      integer  minloc, mobs, kobs
      integer  iargc, argc
      integer  nymd, nhms
      logical  newprof

      character(len=255) infile (NFILES_MAX) ! input filenames
      character(len=255) prefix              ! output filename prefix
      character(len=255) outfile             ! output filename

      logical append_mode
      logical daily
      integer synhour

      type(ods_vect) ods
      integer           :: nobs_ods
      character(len=80) :: ftype
      
!     Parse command line
!     ------------------
      call init ( infile, nfiles_max, nfiles,
     .            synhour, append_mode, daily, prefix )


!     Loop over input files
!     ---------------------
      do ifile = 1, nfiles
          fname = trim(infile(ifile))

          nymd = -1            ! get data for next synoptic time on file
          nhms =  0
          call ODSNxTime ( trim(infile(ifile)), nymd, nhms )
          if ( nymd .eq. -1 ) then 
               print *, 'End-Of-File'
               exit     ! skip to next file
          end if

          print *, 'calling ODS_Get'

!         Read diag_conv file from GSI
!         ------------------
          call ODS_Get ( fname, nymd, nhms, ftype, ods, ier, ncf=.true. )
            if ( ier/=0 ) then
                 print *, 'Error read file ', trim(fname)
                 exit
            endif

          if ( daily ) then
               write(outfile,'(2a,i8.8)') trim(prefix), '.ods.t', nymd
          else
               write(outfile,'(a,''.ods'')') trim(prefix)
          end if

!         Write out ODS file
!         ------------------          
          call ODS_Put ( trim(outfile), ftype, nymd, nhms, ods, ier )
          print *, ' number of obs written to ODS: ', ods%data%nobs

      end do ! < ifile >

      end

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!BOP
! !ROUTINE: init: initialize odsselect
!
! !DESCRIPTION:
!
! !INTERFACE:
!
      subroutine init ( infile, nfiles_max, nfiles,
     .                  synhour,
     .                  append_mode, daily, prefix )

! !INPUT PARAMETERS:
!
      implicit NONE
      integer,       intent(in)  :: nfiles_max

! !INPUT/OUTPUT PARAMETERS:

      character*255, intent(out) :: prefix
      logical,       intent(out) :: daily
      logical,       intent(out) :: append_mode
      integer,       intent(out) :: synhour
!
! !OUTPUT PARAMETERS:

      character*255, intent(out) :: infile(nfiles_max)
      integer,       intent(out) :: nfiles


! !REVISION HISTORY:
!	19Apr2004 - R. Todling - Initial code.
!
!EOP
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

      character*4, parameter :: myname = 'init'

      integer       iret, i, lv, iarg, argc, iargc
      integer        k, k1, k2, ic, lt
      character*255 argv
      character*10  SS

!     Set defaults
!     ------------
      prefix      = 'CONV'   ! DEFAULT: output file names prefixed with 'CONV'
      daily       = .FALSE.  ! DEFAULT: single output file
      append_mode = .FALSE.  ! DEFAULT: no append
      synhour     = -1       ! DEFAULT: process all synoptic times

!     Parse command line
!     ------------------
      argc =  iargc()
      if ( argc .lt. 1 ) call usage()
      nfiles = 0
      iarg = 0
      do i = 1, 32767
         iarg = iarg + 1
         if ( iarg .gt. argc ) go to 111
         call GetArg ( iArg, argv )
         if (index(argv,'-daily' ) .gt. 0 ) then
             daily = .TRUE.
         elseif (index(argv,'-append') .gt. 0 ) then
            append_mode = .TRUE.
         elseif (index(argv,'-synhour') .gt. 0 ) then
            if ( iarg+1 .gt. argc ) call usage()
            iarg = iarg + 1
            call GetArg ( iArg, SS )
            read(SS,*) synhour
         elseif (index(argv,'-o' ) .gt. 0 ) then
            if ( iarg+1 .gt. argc ) call usage()
            iarg = iarg + 1
            call GetArg ( iArg, prefix )
         else
            nfiles = nfiles + 1
            if ( nfiles .gt. nfiles_max ) then
               print *, 'Maximum number of input files = ', nfiles_max
               stop
            end if
            infile(nfiles) = argv
         end if
      end do
 111  continue
      if ( nfiles .lt. 1 ) call usage()

      print *, ' Will process the following files: '
      do i = 1, nfiles
         print *, trim(infile(i)) 
      enddo

      return
      end subroutine init

      subroutine usage
      print *
      print *, 'Usage:'
      print *
      print *, 'dconv2ods.x [-daily] [-append] [-synhour HH] [-o ID] diag_conv_file(s)'
      print *
      print *, 'where'
      print *
      print *,'-daily            create daily output files'
      print *,'                    (default: all output in single file)'
      print *,'-append           append to same synoptic time: '
      print *,'                    CAUTION: must specify -synhour'
      print *,'                    (default: no append)'
      print *,'-synhour HH       process synoptic hour HH only'
      print *,'                    (default: all)'
      print *,'-o ID             use ID for naming output files'
      print *,'                    (default: CONV)'
      print *,' diag_conv(s)     diag_conv file(s) from GSI output'
      print *
      print *
      print *, 'NOTE: diag_conv files assumed to contain obs for a single syn time'
      print *
      stop

      end subroutine usage
