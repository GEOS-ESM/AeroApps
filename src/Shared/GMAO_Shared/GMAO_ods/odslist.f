
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!-----------------------------------------------------------------------
!BOP
!
! !ROUTINE: odslist:  Prints contents of ods files to stdout.
!
! !INTERFACE:

      program odslist

!
! !USAGE: see the routine usage() below
!
! !DESCRIPTION: Lists contents of an ODS file in ASC format.
!
! !TODO:
!
!  option for selecting which attributes to print

! !USES:

      use m_MergeSorts
      use m_ods

      implicit NONE

!
! !REVISION HISTORY:
!
!     21Sep98 (Dee) - original code
!     28Sep98 (Dee) - changed Read_ODS interface
!     10Nov98 (Dee) - only 1 synoptic time, 1 file
!     04Jan99 (Dee) - updated for latest OBS_IO
!     09Dec99 (Todling) - increased nobs_max;
!                       - no halt when attribute not found
!     14Apr00 (Dee) - added 1p to print formats
!     16Feb01 (Todling) - upgraded to function w/ m_ods
!     03Oct01 (Todling) - added xvec.
!     15Jun04 (Todling) - added -ncf to handle GSI files.
!
!EOP
!BOC

      character*7, parameter :: myname = 'odslist'

      integer, parameter :: NFILES_MAX =   1 ! Max number of input files

!     Local variables
!     ---------------
      integer nfiles, ifile, isyn, ksyn, nymd, nhms
      integer i, ierr, nobs, synhour
      logical eof, verbose, post_anal, nonames, nosort, ncf

      character*255 infile (NFILES_MAX) ! input filenames
      character*255 outfile             ! output filename
      character*80  ftype               ! file type: pre- or post-anal

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
     
      integer, allocatable :: indx(:)

!     Option flags:
!     ------------
      synhour = -1         ! DEFAULT: process first synoptic time
      nonames = .FALSE.    ! DEFAULT: do include attribute names
      nosort  = .FALSE.    ! DEFAULT: do sort
      ncf     = .FALSE.    ! DEFAULT: non-complaint format, i.e., non-ODS files
      outfile = '@DEFAULT' ! DEFAULT: output file name (depends on nonames)

!     Parse command line and load resources
!     -------------------------------------
      call init ( infile, NFILES_MAX, nfiles,
     .                            synhour, nonames, nosort, outfile, ncf )

!     Figure out range of synoptic time loop
!     --------------------------------------
      ksyn = 32767
      if(ncf) ksyn = 1

!     Loop over input files
!     ---------------------
      do ifile = 1, nfiles

!       Loop over first synoptic times on this file
!       -------------------------------------
        do isyn = 1, ksyn

!         Read all data for this synoptic time
!         ------------------------------------

          post_anal = .TRUE.  ! can handle post-analysis ODS
	  nymd = -1           ! get data for next synoptic time on file
	  nhms =  0
          call ODSNxTime ( trim(infile(ifile)), nymd, nhms )
          if ( nymd .eq. -1 ) then
               print *, 'End-Of-File'
               exit     ! skip to next file
          end if
          call ODS_Get ( trim(infile(ifile)), nymd, nhms, ftype, ods, ierr, ncf=ncf )

          if ( ierr .gt. 0 ) then
               print *, 'ODS_Get error: ierr = ', ierr
               exit     ! skip to next file
          end if

          if( ftype(1:3)=='pre' ) post_anal = .FALSE.

!         Set number of observations found in the file
!         --------------------------------------------
          nobs   =  ods%data%nobs

          if ( nobs .eq. 0 ) then
               print *, 'No data for this synoptic time'
               cycle    ! skip to next synoptic time
          end if

          if ( synhour .ge. 0 .and. nhms .ne. synhour*10000 ) then
               print *, 'Skipping this synoptic time'
               cycle    ! skip to next synoptic time
          end if

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

          print *, 'number of obs read = ', nobs
          print *, 'date = ', nymd
          print *, 'time = ', nhms
	  
	  if (.not. nosort ) then

             allocate ( indx(nobs) )
	  
	     call IndexSet  ( nobs, indx )
	     call IndexSort ( nobs, indx, qcexcl(1:nobs), descend=.false. )
	     call IndexSort ( nobs, indx,    lat(1:nobs), descend=.false. )
	     call IndexSort ( nobs, indx,    lon(1:nobs), descend=.false. )
	     call IndexSort ( nobs, indx,    lev(1:nobs), descend=.false. )
	     call IndexSort ( nobs, indx,   time(1:nobs), descend=.false. )
	     call IndexSort ( nobs, indx,     kt(1:nobs), descend=.false. )
	     call IndexSort ( nobs, indx,     ks(1:nobs), descend=.false. )
	     call IndexSort ( nobs, indx,     kx(1:nobs), descend=.false. )
  
	     kt    (1:nobs) = kt    ( (/ (indx(i), i=1,nobs) /) )
             kx    (1:nobs) = kx    ( (/ (indx(i), i=1,nobs) /) )
             ks    (1:nobs) = ks    ( (/ (indx(i), i=1,nobs) /) )
             lon   (1:nobs) = lon   ( (/ (indx(i), i=1,nobs) /) )
             lat   (1:nobs) = lat   ( (/ (indx(i), i=1,nobs) /) )
             lev   (1:nobs) = lev   ( (/ (indx(i), i=1,nobs) /) )
             time  (1:nobs) = time  ( (/ (indx(i), i=1,nobs) /) )
             obs   (1:nobs) = obs   ( (/ (indx(i), i=1,nobs) /) )
             OmF   (1:nobs) = OmF   ( (/ (indx(i), i=1,nobs) /) )
             OmA   (1:nobs) = OmA   ( (/ (indx(i), i=1,nobs) /) )
             xm    (1:nobs) = xm    ( (/ (indx(i), i=1,nobs) /) )
             qcexcl(1:nobs) = qcexcl( (/ (indx(i), i=1,nobs) /) )
             qchist(1:nobs) = qchist( (/ (indx(i), i=1,nobs) /) )
             Xvec  (1:nobs) = Xvec  ( (/ (indx(i), i=1,nobs) /) )

             deallocate ( indx )
	  
	  end if

          call ODS_list ( outfile, nobs, nonames, post_anal,
     .                    lat, lon, lev, time,
     .                    kt, kx, ks, xm, qcexcl, qchist,
     .                    obs, omf, oma, xvec )

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

          exit    ! OK, done.

        end do  ! loop over synoptic times

      end do  ! loop over files

      stop

      end ! program odslist

!EOC

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!-----------------------------------------------------------------------
!BOP
! !ROUTINE: init: initialize odslist
!
! !DESCRIPTION:
!
! !INTERFACE:
!
      subroutine init ( infile, NFILES_MAX, nfiles,
     .                          synhour, nonames, nosort, outfile, ncf )

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
      logical,       intent(inout) :: nonames
      logical,       intent(inout) :: nosort
      character*255, intent(inout) :: outfile
      logical,       intent(inout) :: ncf
!
!
! !REVISION HISTORY:
!
!       21Sep98 - D.Dee   Initial code
!
!EOP
!BOC

      character*14, parameter :: myname = 'init'

      integer iret, i, lv, iarg, argc, iargc
      character*255 argv
      character*2 HH

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
         if (index(argv,'-nonames') .gt. 0 ) then ! non-documented
            nonames = .TRUE.
         else if (index(argv,'-odt') .gt. 0 ) then
            nonames = .TRUE.
         else if (index(argv,'-strip') .gt. 0 ) then
            nonames = .TRUE.
         elseif (index(argv,'-nosort') .gt. 0 ) then
            nosort = .TRUE.
         elseif (index(argv,'-ncf') .gt. 0 ) then
            ncf = .TRUE.
         elseif (index(argv,'-synhour') .gt. 0 ) then
            if ( iarg+1 .gt. argc ) call usage()
            iarg = iarg + 1
            call GetArg ( iArg, HH )
            read(HH,*) synhour
         elseif (index(argv,'-o') .gt. 0 ) then
            if ( iarg+1 .gt. argc ) call usage()
            iarg = iarg + 1
            call GetArg ( iArg, outfile )
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

      if ( trim(outfile) .eq. '@DEFAULT' ) then
           if ( nonames ) then
                outfile = 'odslist.odt'
           else
                outfile = 'odslist.txt'
           end if
      end if

!     Echo the parameters
!     -------------------
      print *
      if (synhour .ge. 0) print *, 'Will process synoptic hour ',
     .                             synhour, ' only'
      print *
      print *, 'Input files: ', nfiles
      print *
      do i = 1, nfiles
         lv = len_trim(infile(i))
         print *, ' o ', infile(i)(1:lv)
      end do
      print *
      print *, 'Output filename: ', outfile

      return

      end ! subroutine init

!EOC

!-------------------------------------------------------------------------
!         NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!-------------------------------------------------------------------------
!BOP
!
! !ROUTINE: ODS_list --- Prints a listing of ODS attributes to stdout
!
! !INTERFACE:
!
      subroutine ODS_list(outfile, nobs, nonames, postanal,
     .                          lat, lon, lev, time,
     .                          kt, kx, ks, xm, qcx, qch,
     .                          obs, omf, oma, xvec )

! !DESCRIPTION:
!
! !REVISION HISTORY:
!
!EOP
!-------------------------------------------------------------------------
!BOC

      character*255 outfile
      integer nobs
      logical nonames
      logical postanal
      real    lat(nobs)
      real    lon(nobs)
      real    lev(nobs)
      integer time(nobs)
      integer kt(nobs)
      integer kx(nobs)
      integer ks(nobs)
      real    xm(nobs)
      integer qcx(nobs)
      integer qch(nobs)
      real    obs(nobs)
      real    omf(nobs)
      real    oma(nobs)
      real    xvec(nobs)

      integer i, ios

      open(99, file=outfile,
     .         form='formatted',
     .         iostat=ios)
      if ( ios .ne. 0 ) then
           print *, 'Cannot open file ', outfile
           return
      else
           print *, 'Writing ODS to file ', outfile
      end if

      if ( postanal) then

         if ( .not. nonames ) then

            do i = 1, nobs
               write(99,'("kx=",i3," ks=",i9," kt=",i2,
     .         " time=",i5," lev=",f6.1," lon=",f7.2," lat=",f7.2,
     .         " obs=",1p,e9.2," omf=",e9.2," oma=",e9.2,
     .         " qcx=",i3," qch=",i3," xm=",e9.2," xvec=",e9.2)')
     .              kx(i), ks(i), kt(i),
     .              time(i), lev(i), lon(i), lat(i),
     .              obs(i), omf(i), oma(i),
     .              qcx(i), qch(i), xm(i), xvec(i)
            end do

         else

            do i = 1, nobs
               write(99,'(i3," ",i9," ",i2,
     .         " ",i5," ",f6.1," ",f7.2," ",f7.2,
     .         " ",1p,e9.2," ",e9.2," ",e9.2,
     .         " ",i3," ",i3," ",e9.2," ",e9.2)')
     .              kx(i), ks(i), kt(i),
     .              time(i), lev(i), lon(i), lat(i),
     .              obs(i), omf(i), oma(i),
     .              qcx(i), qch(i), xm(i), xvec(i)
            end do

         end if

      else

         if ( .not. nonames ) then

            do i = 1, nobs
               write(99,'("kx=",i3," ks=",i9," kt=",i2,
     .         " time=",i5," lev=",f6.1," lon=",f7.2," lat=",f7.2,
     .         " obs=",1p,e9.2,
     .         " qcx=",i3," qch=",i3," xm=",e9.2)')
     .              kx(i), ks(i), kt(i),
     .              time(i), lev(i), lon(i), lat(i),
     .              obs(i),
     .              qcx(i), qch(i), xm(i)
            end do

         else

            do i = 1, nobs
               write(99,'(i3," ",i9," ",i2,
     .         " ",i5," ",f6.1," ",f7.2," ",f7.2,
     .         " ",1p,e9.2,
     .         " ",i3," ",i3," ",e9.2)')
     .              kx(i), ks(i), kt(i),
     .              time(i), lev(i), lon(i), lat(i),
     .              obs(i),
     .              qcx(i), qch(i), xm(i)
            end do

         end if

      end if

      return
      end
!EOC

      subroutine usage()
      print *
      print *, 'odslist - Create ASCII Listing from ODS File'
      print *
      print *, 'Usage:  odslist [-strip] [-nosort] [-ncf]',
     .                     '[-synhour HH] [-o FNAME] odsfile'
      print *
      print *, 'where'
      print *
      print *, '-strip          strip attribute names from output file,'
      print *, '                creating a standard ODT file for use'
      print *, '                with ods_maker.x'
      print *, '                (default: attribute names are included)'
      print *, '-nosort          do not sort'
      print *, '                (default: do)'
      print *, '-ncf             specify when input files are non-ODS files'
      print *, '                (default: ignore)'
      print *, '-synhour HH      process synoptic hour HH'
      print *, '                (default: first on file)'
      print *, '-o FNAME         output file name'
      print *, '                (default: odslist.txt or odslist.odt)'
      print *, ' -odt           same as -strip'
      print *, ' odsfile         ODS file'
      print *
      stop
      end

