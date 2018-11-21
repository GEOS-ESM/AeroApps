!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!-----------------------------------------------------------------------
!BOP
!
! !ROUTINE: odsemean:  Operates on specific attributes in ODS
!
! !INTERFACE:

      program odsemean

!
! !DESCRIPTION:
!           Reads ODS files, selects data according to criteria specified
!           via a resource file, operates on chosen attribute and summarizes
!           results by writing them to an ASC file.
!

      use m_ods
      use m_inpak90
      use m_chars,   only: lowercase
      use m_stdio,   only: stdout, stderr

      implicit NONE

!
! !REVISION HISTORY:
!
!     17Aug2017 Todling  Initial code
!
!EOP
!BOC

      character(len=*), parameter :: myname = 'odsemean'

      character(len=*), parameter :: egress = 'ODSEMEAN_EGRESS'
      integer, parameter :: NFILES_MAX = 500 ! Max number of input files
      integer, parameter ::  NTMAX = 12*31*4 ! 4 syn times a day, for a year (roughly)

      logical, parameter :: debug = .false.
      logical            :: verb  = .false.

!     Local variables
!     ---------------
      integer nfiles, ifile, lf, isyn, ksyn, nymd, nhms, nkt, nkx, nlv, ncfound
      integer i, j, ic, nc, nt, ierr, nobs, nsel, synhour, nop, nops

      character*255 infile (NFILES_MAX) ! input filenames
      character*255 fileout, outfile, outodsfn   ! output filename
      character*80  ftype
      logical       anotherclass, allkxs, allevs, lstdv

!     Attributes to accummulate
!     -------------------------
      real,    allocatable :: hxf(:)
      integer, allocatable :: qcx(:)

!     storage for ODS:
!     ---------------
      type ( ods_vect ) ods
      type ( ods_vect ) odss

!     Option flags:
!     ------------
      synhour = -1       ! DEFAULT: process all synoptic times

!     Parse in command line
!     ---------------------
      call init ( infile, nfiles_max, nfiles, nymd, nhms, outfile, verb )

!     Figure out range of synoptic time loop
!     --------------------------------------
      ksyn = 32767
      nt   = 0

!     Loop over input files
!     ---------------------
      do ifile = 1, nfiles

!         Read all data for this synoptic time
!         ------------------------------------
          if(verb) print *, 'calling ODS_Get'

          call ODS_Get ( trim(infile(ifile)), nymd, nhms, ftype, ods, ierr )

          if(verb) print *, 'completed ODS_Get'

          if ( ierr .gt. 0 ) then
               print *, 'ODS_Get error: ierr = ', ierr
               exit     ! skip to next file
          end if

!         Set number of observations found in the file
!         --------------------------------------------
          nobs   =  ods%data%nobs

          if ( nobs .eq. 0 ) then
               print *, 'No data for this synoptic time'
               cycle    ! skip to next synoptic time
          end if

          if ( synhour .ge. 0 .and. nhms .ne. synhour*10000 ) then
               if(verb) print *, 'Skipping this synoptic time'
               cycle    ! skip to next synoptic time
          end if

          nt = nt + 1
          if ( nt>ntmax) then
              print *, 'Max number of syn hours reached'
              print *, 'Need to lower number if processed files'
              stop
          endif

          if (verb) then
              print *, 'number of obs read = ', nobs
              print *, 'date = ', nymd
              print *, 'time = ', nhms
          endif

          if (ifile==1) then
             call ODS_Select ( ods, nobs, nsel, ierr, odss=odss )
             if ( nsel /= nobs ) then
                  print *, 'ODS_Select error: no match in obs count ', nobs, nsel
             endif
             allocate(hxf(nobs))
             allocate(qcx(nobs))
             hxf=0.0
             qcx=ods%data%qcexcl
          endif

!         Accummulate hxf
!         ---------------
          where (qcx==0)  ! avoid operating on UNDEF (missing values)
             hxf = hxf + ods%data%obs-ods%data%omf
          endwhere

          call ODS_Clean ( ods, ierr )
             if ( ierr .ne. 0 ) then
                  print *, 'ODS_Clean error: ierr = ', ierr
             endif

      end do  ! loop over files

!     Wherever makes sense ...
      where (qcx==0)

!        Calculate mean hxf
!        ------------------
         hxf = hxf/real(nfiles)

!        Recalculate omf with mean hxf
!        -----------------------------
         odss%data%omf = odss%data%obs - hxf

      elsewhere
         odss%data%omf = obs_missing
      endwhere

      if (trim(outfile) == 'odsemean') then
         write(outodsfn,'(2a,i8.8,a,i2.2,a)') trim(outfile), '.obs.', nymd, '_', nhms/10000, 'z.ods'
      else
         outodsfn=trim(outfile)
      endif
      print *, 'Writing ods file for debug purposes: ', trim(outodsfn)
      call ODS_Put ( trim(outodsfn), ftype, nymd, nhms, odss, ierr )

      if ( ierr==0 ) then
         close(999)
         open (999,file=trim(egress),form='formatted')
         close(999)
      endif

!     Clean up
!     --------
      deallocate(hxf)
      deallocate(qcx)
      end program odsemean ! program odsemean

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!-----------------------------------------------------------------------
!BOP
! !ROUTINE: init: initialize odstats
!
! !DESCRIPTION:
!
! !INTERFACE:
!
      subroutine init ( infile, nfiles_max, nfiles, nymd, nhms, 
     .                  outfile, verb )

! !USES:

! !INPUT PARAMETERS:
!
      implicit NONE
      integer,       intent(in)  :: nfiles_max

! !OUTPUT PARAMETERS:

      character(len=*), intent(out) :: infile(nfiles_max)
      integer,          intent(out) :: nfiles
      integer,          intent(out) :: nymd, nhms
      character(len=*), intent(out) :: outfile
      logical,          intent(out) :: verb
!
!
! !REVISION HISTORY:
!     17Aug2017 Todling - Initial code (stripped off odsstats)
!
!EOP
!BOC

      character*4, parameter :: myname_ = 'init'

      integer iret, i, ic, lt, lv, iarg, argc, iargc
      real swap
      character*255 argv
      character*10 SS

      outfile = 'odsemean'
      verb = .false.
      nymd = -1
      nhms = -1

!     Parse command line
!     ------------------

      argc =  iargc()
      if ( argc .lt. 1 ) call usage()
      nfiles = 0
      iarg = 0
      ic = 0
      do i = 1, 32767
         iarg = iarg + 1
         if ( iarg .gt. argc ) go to 111
         call GetArg ( iArg, argv )
         if (index(argv,'-o' ) .gt. 0 ) then
            if ( iarg+1 .gt. argc ) call usage()
            iarg = iarg + 1
            call GetArg ( iArg, outfile )
         else if (index(argv,'-verbose' ) .gt. 0 ) then
            verb = .true.
         else
            ic = ic + 1
            if(ic==1) read(argv,*) nymd
            if(ic==2) read(argv,*) nhms
            if(ic>2) then
               nfiles = nfiles + 1
               if ( nfiles .gt. nfiles_max ) then
                  print *, 'Maximum number of input files = ', nfiles_max
                  stop
               end if
               infile(nfiles) = argv
            endif
         end if
      end do
 111  continue
      if ( nfiles .lt. 1 ) call usage()

      print *
      write(6,'(a,i8.8,a,i6.6)') 'On date: ', nymd,' and time ', nhms
      print *, 'Input files: ', nfiles
      print *
      do i = 1, nfiles
         lv = len_trim(infile(i))
         print *, ' o ', infile(i)(1:lv)
      end do
      print *
      print *, 'Output filename: ', trim(outfile)

      return

      end subroutine init! subroutine init

!EOC

!-------------------------------------------------------------------------

      subroutine usage()
      print *
      print *, 'Usage:'
      print *
      print *, 'obsemean [-o ID] nymd nhms odsfile(s)'
      print *
      print *, 'where'
      print *
      print *,'-o  ID         use ID for naming output files'
      print *,'                (default: obsemean.ods)'
      print *,'-verbose       sets verbose on (default: off)'
      print *
      print *
      print *,' odsfile(s)    ODS file(s)'
      print *
      print *
      stop
      end
