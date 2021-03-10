!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 610.1, GEOS/DAS      !
!-----------------------------------------------------------------------
!BOP
!
! !ROUTINE: odsshuffle:  Shuffles ODS more or less randomly.
!
! !INTERFACE:

      program odsshuffle

!
! !USAGE:   see the routine usage() below

! !USES:

      use m_MergeSorts
      use m_ods
      use m_random
      implicit NONE

! !DESCRIPTION: Reads an existing ODS file; shuffles its contents
!               somewhat, and output another ODS file.
!
! !TODO:
!
!     implement a true random shuffle, using random number generator
!     and sort algorithm.
!
! !REVISION HISTORY:
!
!     21Sep98 (Dee) - original code
!     28Sep98 (Dee) - changed Read_ODS interface
!     04Jan99 (Dee) - updated for latest OBS_IO
!     05Feb99 (Dee) - truly random
!     09Dec99 (Todling) - increased nobs_max;
!                       - no halt when attribute not found
!     30Oct00 (Todling) - increased max kx
!     21Feb01 (Todling) - modified to function w/ m_ods.
!     06Jun02 (Todling) - added ref to m_random
!
!EOP
!BOC

      character*10, parameter :: myname = 'odsshuffle'

      integer, parameter :: NFILES_MAX = 100 ! Max number of input files

!     Local variables
!     ---------------
      integer nfiles, ifile, lf, isyn, nymd, nhms
      integer ierr, nobs, synhour, seed, i
      logical eof, verbose, post_anal, append_mode, daily

      character*255 infile (NFILES_MAX) ! input filenames
      character*255 prefix              ! output filename prefix
      character*255 outfile             ! output filename
      character*80 ftype

      integer, allocatable :: indx (:)  ! permutation index for shuffle
      real*8, allocatable :: rndm (:)  ! random numbers

!     storage for ODS:
!     ---------------
      type ( ods_vect ) ods

      real   , pointer :: lat   (:)
      real   , pointer :: lon   (:)
      real   , pointer :: lev   (:)
      integer, pointer :: time  (:)
      integer, pointer :: kt    (:)
      integer, pointer :: kx    (:)
      integer, pointer :: ks    (:)
      real   , pointer :: xm    (:)
      integer, pointer :: qcexcl(:)
      integer, pointer :: qchist(:)
      real   , pointer :: obs   (:)
      real   , pointer :: omf   (:)
      real   , pointer :: oma   (:)
      real   , pointer :: xvec  (:)

!     Option flags:
!     ------------
      synhour = -1         ! DEFAULT: process all synoptic times
      prefix  = 'SHUFFLE'  ! DEFAULT: output file names prefixed with 'SHUFFLE'
      seed    = 0          ! DEFAULT: initial seed for random number generator
      daily   = .TRUE.

!     Parse command line and load resources
!     -------------------------------------
      call init ( infile, NFILES_MAX, nfiles, prefix, synhour, seed )

      call zufalli ( seed )       ! initialize random number generator

!     Loop over input files
!     ---------------------
      do ifile = 1, nfiles

!       Loop over all synoptic times on this file
!       -----------------------------------------
        do isyn = 1, 32767

!         Read all data for this synoptic time
!         ------------------------------------

          post_anal = .TRUE.   ! can handle post-analysis ODS
	  nymd = -1            ! get data for next synoptic time on file
	  nhms =  0
          call ODSNxTime ( trim(infile(ifile)), nymd, nhms )
          if ( nymd .eq. -1 ) then
               print *, 'End-Of-File'
               exit     ! skip to next file
          end if

          print *, 'calling ODS_Get'

          call ODS_Get ( trim(infile(ifile)), nymd, nhms, ftype, ods, ierr )

          print *, 'completed ODS_Get'

!         Set number of observations found in the file
!         --------------------------------------------
          nobs   =  ods%data%nobs

          if ( ierr .gt. 0 ) then
               print *, 'ODS_Get error: ierr = ', ierr
               exit     ! skip to next file
          end if

          if ( nobs .eq. 0 ) then
               print *, 'No data for this synoptic time'
               cycle    ! skip to next synoptic time
          end if

          if ( synhour .ge. 0 .and. nhms .ne. synhour*10000 ) then
               print *, 'Skipping this synoptic time'
               cycle    ! skip to next synoptic time
          end if

          if( ftype(1:3)=='pre' ) post_anal = .FALSE.

!         Set pointers
!         ------------
          lat    => ods%data%lat
          lon    => ods%data%lon
          lev    => ods%data%lev
          time   => ods%data%time
          kt     => ods%data%kt
          kx     => ods%data%kx
          ks     => ods%data%ks
          xm     => ods%data%xm
          qchist => ods%data%qchist
          qcexcl => ods%data%qcexcl
          obs    => ods%data%obs
          omf    => ods%data%omf
          oma    => ods%data%oma
          xvec   => ods%data%xvec

          allocate ( rndm(nobs), indx(nobs), stat=ierr )
          if ( ierr .ne. 0 ) then
               print *, 'Error, Alloc(rndn,indx), ierr = ', ierr
               call exit (2)
          end if

          print *, 'number of obs read = ', nobs
          print *, 'date = ', nymd
          print *, 'time = ', nhms

!         reorder data
!         ------------

          print *, 'shuffling input data'
          call normalen ( nobs, rndm )

          call IndexSet  ( nobs, indx )
          call IndexSort ( nobs, indx, rndm(1:nobs), descend=.false. )

          kt    (1:nobs) = kt    ( (/ (indx(i), i=1,nobs) /) )
          kx    (1:nobs) = kx    ( (/ (indx(i), i=1,nobs) /) )
          ks    (1:nobs) = ks    ( (/ (indx(i), i=1,nobs) /) )
          lon   (1:nobs) = lon   ( (/ (indx(i), i=1,nobs) /) )
          lat   (1:nobs) = lat   ( (/ (indx(i), i=1,nobs) /) )
          lev   (1:nobs) = lev   ( (/ (indx(i), i=1,nobs) /) )
          time  (1:nobs) = time  ( (/ (indx(i), i=1,nobs) /) )
          obs   (1:nobs) = obs   ( (/ (indx(i), i=1,nobs) /) )
          xm    (1:nobs) = xm    ( (/ (indx(i), i=1,nobs) /) )
          qcexcl(1:nobs) = qcexcl( (/ (indx(i), i=1,nobs) /) )
          qchist(1:nobs) = qchist( (/ (indx(i), i=1,nobs) /) )

          if ( post_anal ) then
             OmF   (1:nobs) = OmF   ( (/ (indx(i), i=1,nobs) /) )
             OmA   (1:nobs) = OmA   ( (/ (indx(i), i=1,nobs) /) )
             Xvec  (1:nobs) = Xvec  ( (/ (indx(i), i=1,nobs) /) )
          end if

!         Write to ods file
!         -----------------
          print *, 'calling ODS_Put'

          append_mode = .FALSE.
          if ( daily ) then
               write(outfile,'(2a,i8.8)') trim(prefix), '.ods.t', nymd
          else
               write(outfile,'(a,''.ods'')') trim(prefix)
          end if

          ods%data%nobs = nobs   ! write out only select part
          call ODS_Put ( outfile, ftype, nymd, nhms, ods, ierr,
     .                   append=append_mode )
          print *, 'completed ODS_Put'

          if ( ierr .ne. 0 ) then
               print *, 'ODS_Put error: ierr = ', ierr
          end if

!         Nullify pointers
!         -----------------
          if(associated(lat)   ) nullify(lat)
          if(associated(lon)   ) nullify(lon)
          if(associated(lev)   ) nullify(lev)
          if(associated(time)  ) nullify(time)
          if(associated(kt)    ) nullify(kt)
          if(associated(kx)    ) nullify(kx)
          if(associated(ks)    ) nullify(ks)
          if(associated(xm)    ) nullify(xm)
          if(associated(qchist)) nullify(qchist)
          if(associated(qcexcl)) nullify(qcexcl)
          if(associated(obs)   ) nullify(obs)
          if(associated(omf)   ) nullify(omf)
          if(associated(oma)   ) nullify(oma)
          if(associated(xvec)  ) nullify(xvec)

          deallocate ( rndm, indx )

        end do  ! loop over synoptic times

      end do  ! loop over files

      stop

      end ! program odsshuffle

!EOC

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 610.1, GEOS/DAS      !
!-----------------------------------------------------------------------
!BOP
! !ROUTINE: init: initialize odsshuffle
!
! !DESCRIPTION:
!
! !INTERFACE:
!
      subroutine init ( infile, NFILES_MAX, nfiles, prefix, synhour, seed )

! !INPUT PARAMETERS:
!
      implicit NONE
      integer,       intent(in)  :: nfiles_max
!
! !OUTPUT PARAMETERS:

      character*255, intent(out) :: infile(nfiles_max)
      integer,       intent(out) :: nfiles

! !INPUT/OUTPUT PARAMETERS:

      integer,       intent(inout) :: synhour
      integer,       intent(inout) :: seed
      character*255, intent(inout) :: prefix
!
!
! !REVISION HISTORY:
!
!       21Sep98 - D.Dee   Initial code
!
!EOP
!BOC

      character*4, parameter :: myname = 'init'

      integer iret, i, lv, iarg, argc, iargc
      character*255 argv
      character*2 HH
      character*10 SS

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
         if (index(argv,'-prefix' ) .gt. 0 ) then
            if ( iarg+1 .gt. argc ) call usage()
            iarg = iarg + 1
            call GetArg ( iArg, prefix )
         elseif (index(argv,'-seed') .gt. 0 ) then
            if ( iarg+1 .gt. argc ) call usage()
            iarg = iarg + 1
            call GetArg ( iArg, SS )
            read(SS,*) seed
         elseif (index(argv,'-synhour') .gt. 0 ) then
            if ( iarg+1 .gt. argc ) call usage()
            iarg = iarg + 1
            call GetArg ( iArg, HH )
            read(HH,*) synhour
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

!     Echo the parameters
!     -------------------
      print *
      if (synhour .ge. 0) print *, 'Will process synoptic hour ',
     .                             synhour, ' only'
      print *
      print *, 'Initial seed for random number generator is ', seed
      print *
      print *, 'Input files: ', nfiles
      print *
      do i = 1, nfiles
         lv = len_trim(infile(i))
         print *, ' o ', infile(i)(1:lv)
      end do
      print *
      lv = len_trim(prefix)
      print *, 'Output filename prefix: ', prefix(1:lv)

      return

      end ! subroutine init

!EOC

      subroutine usage()
      print *
      print *, 'Usage:  odsshuffle [-prefix ID] [-synhour HH] [-seed SS] ',
     .                               'odsfile(s)'
      print *
      print *, 'where'
      print *
      print *,   '-prefix ID     use ID as a prefix for output file names'
      print *,   '              (default: SHUFFLE)'
      print *,   '-synhour HH    process synoptic hour HH only'
      print *,   '              (default: all)'
      print *,   '-seed SS       use SS as seed for random number generator'
      print *,   '              (default: 0)'
      print *,   ' odsfile(s)    ODS file(s)'
      print *
      stop
      end

