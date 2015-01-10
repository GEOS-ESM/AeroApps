!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!-----------------------------------------------------------------------
!BOP
!
! !ROUTINE: odsselect:  Selects from ODS files.
!
! !INTERFACE:

      program odsselect

!
! !DESCRIPTION:
!           Reads ODS files, selects data according to criteria specified
!           on the command line, and writes new ODS files containing the
!           selected data only.
!
! !USAGE:   see the routine usage() below
!

      use m_ods

      implicit NONE

!
! !REVISION HISTORY:
!
!     23Sep98 (Dee) - original code
!     28Sep98 (Dee) - changed Read_ODS interface
!     08Dec98 (Dee) - added 'daily' option
!     04Jan99 (Dee) - updated for latest OBS_IO
!     09Dec99 (Todling) - increased nobs_max;
!                       - no halt when attribute not found
!     10Apr00 (Dee) - added 'lon' and 'lat' options,
!                         modified kx, kt, o options.
!     30Oct00 (Todling) - increased max kx.
!     17Feb01 (Todling) - modified to function w/ m_ods.
!     22Mar01 (Todling) - Added capability to handle multiple, 
!                         non-consecutive levels
!     03May01 (Todling) - Added opt to replace qcexcl flag
!     30Jan02 (Todling) - Added opt to select by obs time.
!     15Jun04 (Todling) - added -ncf to handle GSI files.
!     07Dec04 (Dee)     - added -rmdupl to remove duplicates
!     20Jul05 (Todling) - added ods_clean to solve mem leak issue
!
!EOP
!BOC

      character*9, parameter :: myname = 'odsselect'

      integer, parameter :: NFILES_MAX = 512 ! Max number of input files
      integer, parameter :: SMAX     = 100   ! Max number of arguments
      integer, parameter :: MAXLEV   = 30    ! Max number of levs to select

!     Local variables
!     ---------------
      integer nfiles, ifile, lf, isyn, ksyn, nymd, nhms, nkt, nkx, nqcx
      integer i, j, ierr, ierr_att(13), nobs, synhour, is
      integer miter, jiter
      logical post_anal, append_mode, eof, verbose, daily
      logical kxselect, ktselect, levselect, qcselect
      logical lonselect, latselect, timeselect, ncf, rmdupl
      integer nkxs, nkts
      integer iqcx
      integer kxs(SMAX), kts(SMAX)
      integer nlev
      real    levs(MAXLEV)
      real    lev1, lev2, lat1, lat2, lon1, lon2 
      integer t1, t2
      character*4 qcs
      character*80 ftype 

      character*255 infile (NFILES_MAX) ! input filenames
      character*255 prefix              ! output filename prefix
      character*255 outfile             ! output filename

      integer, allocatable :: iselect(:) ! indices of selected data

!     storage for ODS:
!     ---------------
      type ( ods_vect ) ods

      real   , pointer :: lat (:)
      real   , pointer :: lon (:)
      real   , pointer :: lev (:)
      integer, pointer :: time(:)
      integer, pointer :: kt  (:)
      integer, pointer :: kx  (:)
      integer, pointer :: ks  (:)
      real   , pointer :: xm  (:)
      integer, pointer :: qcx (:)
      integer, pointer :: qch (:)
      real   , pointer :: obs (:)
      real   , pointer :: omf (:)
      real   , pointer :: oma (:)
      real   , pointer :: xvec(:)

!     Option flags:
!     ------------
      nkxs = 0            ! DEFAULT: all kx
      nkts = 0            ! DEFAULT: all kt
      lat1 = 0.0
      lat2 = -1.0         ! DEFAULT: all latitudes
      lon1 = 0.0
      lon2 = -1.0         ! DEFAULT: all longitudes
      lev1 = 0.0  
      lev2 = -1.0         ! DEFAULT: all levels
      nlev = 0
      iqcx = -1           ! DEFAULT: leave qcx flag alone
      t1  =  999	  ! DEFAULT: all times
      t2  = -999
      qcs = 'GOOD'        ! DEFAULT: only data that passed QC
      prefix  = 'SELECT'  ! DEFAULT: output file names prefixed with 'SELECT'
      synhour = -1        ! DEFAULT: process all synoptic times
      daily   = .FALSE.   ! DEFAULT: single output file
      append_mode = .FALSE.! DEFAULT: no append
      ncf     = .FALSE.    ! DEFAULT: non-complaint format, i.e., non-ODS files
      rmdupl  = .FALSE.    ! DEFAULT: do not eliminate duplicates
      miter = -1           ! DEFAULT: not expecting sensitivities in diag files
      jiter = -1           ! DEFAULT: not expecting sensitivities in diag files

!     Parse command line and load resources
!     -------------------------------------
      call init ( infile, NFILES_MAX, nfiles,
     .            SMAX,   kxs, nkxs, kts, nkts, lev1, lev2, nlev, levs, qcs,
     .                    lat1, lat2, lon1, lon2, t1, t2, 
     .                    append_mode, daily, prefix, synhour,
     .                    iqcx, miter, jiter, ncf, rmdupl )

!     Figure out range of synoptic time loop
!     --------------------------------------
      ksyn = 32767
      if(ncf) ksyn = 1

!     Loop over input files
!     ---------------------
      do ifile = 1, nfiles

!       Loop over all synoptic times on this file
!       -----------------------------------------
        do isyn = 1, ksyn

!         Read all data for this synoptic time
!         ------------------------------------

          post_anal = .TRUE.   ! can handle post_analysis ODS
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
               print *, 'Skipping this synoptic time'
               cycle    ! skip to next synoptic time
          end if

          if( ftype(1:3)=='pre' ) post_anal = .FALSE.

!         Set pointers
!         ------------
          lat  => ods%data%lat
          lon  => ods%data%lon
          lev  => ods%data%lev
          time => ods%data%time
          kt   => ods%data%kt
          kx   => ods%data%kx
          ks   => ods%data%ks
          xm   => ods%data%xm
          qch  => ods%data%qchist
          qcx  => ods%data%qcexcl
          obs  => ods%data%obs
          omf  => ods%data%omf
          oma  => ods%data%oma
          xvec => ods%data%xvec

          print *, 'number of obs read = ', nobs
          print *, 'date = ', nymd
          print *, 'time = ', nhms


!         Select data
!         -----------
          allocate ( iselect(nobs) )

!         Create a list of indices of selected data

          is = 0
          do i = 1, nobs

             if (nkxs .eq. 0) then
                 kxselect = .TRUE.
             else
                 kxselect = .FALSE.
                 do j = 1, nkxs
                    if (kx(i) .eq. kxs(j)) kxselect = .TRUE.
                 end do
             end if

             if (nkts .eq. 0) then
                 ktselect = .TRUE.
             else
                 ktselect = .FALSE.
                 do j = 1, nkts
                    if (kt(i) .eq. kts(j)) ktselect = .TRUE.
                 end do
             end if

             if (lat2<lat1) then
                 latselect = .TRUE.
             elseif (lat1<=lat(i) .and. lat(i)<=lat2) then
                 latselect = .TRUE.
             else
                 latselect = .FALSE.
             end if
                 
             if (lon2<lon1) then
                 lonselect = .TRUE.
             elseif (lon1<=lon(i) .and. lon(i)<=lon2) then
                 lonselect = .TRUE.
             else
                 lonselect = .FALSE.
             end if

             if (t2<t1) then
                 timeselect = .TRUE.
             elseif (t1<=time(i) .and. time(i)<=t2) then
                 timeselect = .TRUE.
             else
                 timeselect = .FALSE.
             end if
                 
             if (lev2<lev1) then
                 levselect = .TRUE.
             elseif (lev1<=lev(i) .and. lev(i)<=lev2) then
                 levselect = .TRUE.
             else
                 levselect = .FALSE.
             end if
             if ( count(lev(i)==levs(1:nlev)) .ne. 0 ) levselect = .TRUE.

             qcselect = .TRUE.
             if (qcs .eq. 'GOOD') then
                 if (qcx(i) .ne. 0) qcselect = .FALSE.
             elseif (qcs .eq. 'BAD') then
                 if (qcx(i) .eq. 0) qcselect = .FALSE.
             end if

             if ( kxselect .and. ktselect .and.
     .           lonselect .and. latselect .and.
     .           timeselect.and.
     .           levselect .and. qcselect      ) then
                 is = is + 1
                 iselect(is) = i
             end if
          end do

          nobs = is    ! number of selected data

          if (nobs .eq. 0) then
              print *, 'No data satisfy selection criteria.'
              deallocate ( iselect )
              cycle    ! skip to next synoptic time
          else
              print *, 'number of obs selected = ', nobs
          end if

!         move selected data to the front (iselect is monotone)

          do is = 1, nobs
             i  = iselect(is)
             lat(is)  = lat(i)
             lon(is)  = lon(i)
             lev(is)  = lev(i)
             time(is) = time(i)
             kx(is)   = kx(i)
             kt(is)   = kt(i)
             ks(is)   = ks(i)
             xm(is)   = xm(i)
             qcx(is)  = qcx(i)
             qch(is)  = qch(i)
             obs(is)  = obs(i)
          end do
          if( iqcx>=0 ) where(qcx(1:nobs)==iqcx) qcx(1:nobs) = 0

          if ( post_anal ) then
               do is = 1, nobs
                  i  = iselect(is)
                  omf(is)  = omf(i)
                  oma(is)  = oma(i)
                  xvec(is) = xvec(i)
               end do
          end if
	  
          ods%data%nobs = nobs
	  
!         Eliminate duplicates, based on kx,kt,time,lat,lon,lev,obs
!         ---------------------------------------------------------
          if (rmdupl) then
	  
	    call ods_rmdupl ( ods )
	    
            print *, 'eliminated ', nobs-ods%data%nobs, ' duplicate observations'
	    nobs = ods%data%nobs
            print *, 'number of obs remaining = ', nobs
	    
	  end if

!         Write to ods file
!         -----------------
          print *, 'calling ODS_Put'

	  if ( daily ) then
	       write(outfile,'(2a,i8.8)') trim(prefix), '.ods.t', nymd
	  else
	       write(outfile,'(a,''.ods'')') trim(prefix)
	  end if
	  
          call ODS_Put ( outfile, ftype, nymd, nhms, ods, ierr, 
     .                   append=append_mode )
          print *, 'completed ODS_Put'

          if ( ierr .ne. 0 ) then
               print *, 'ODS_Put error: ierr = ', ierr
          end if

!         Nullify pointers
!         -----------------
          if(associated(lat) ) nullify(lat)
          if(associated(lon) ) nullify(lon)
          if(associated(lev) ) nullify(lev)
          if(associated(time)) nullify(time)
          if(associated(kt)  ) nullify(kt)
          if(associated(kx)  ) nullify(kx)
          if(associated(ks)  ) nullify(ks)
          if(associated(xm)  ) nullify(xm)
          if(associated(qch) ) nullify(qch)
          if(associated(qcx) ) nullify(qcx)
          if(associated(obs) ) nullify(obs)
          if(associated(omf) ) nullify(omf)
          if(associated(oma) ) nullify(oma)
          if(associated(xvec)) nullify(xvec)

          deallocate ( iselect )

          call ODS_Clean ( ods, ierr )
             if ( ierr .ne. 0 ) then
                  print *, 'ODS_Clean error: ierr = ', ierr
             endif
        end do  ! loop over synoptic times

      end do  ! loop over files

      stop

      end ! program odsselect

!EOC

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!-----------------------------------------------------------------------
!BOP
! !ROUTINE: init: initialize odsselect
!
! !DESCRIPTION:
!
! !INTERFACE:
!
      subroutine init ( infile, nfiles_max, nfiles,
     .                  smax, kxs, nkxs, kts, nkts, lev1, lev2, nlev, levs,
     .                  qcs, lat1, lat2, lon1, lon2, t1, t2,
     .                  append_mode, daily, prefix, synhour,
     .                  iqcx, miter, jiter, ncf, rmdupl )

! !USES:

      use m_ods_obsdiags, only : ods_obsdiags_init
      use m_ods_obsdiags, only : ods_obsdiags_setparam

! !INPUT PARAMETERS:
!
      implicit NONE
      integer,       intent(in)  :: nfiles_max
      integer,       intent(in)  :: smax

! !INPUT/OUTPUT PARAMETERS:

      character*255, intent(inout) :: prefix
      integer,       intent(inout) :: synhour
      integer,       intent(inout) :: nkxs
      integer,       intent(inout) :: nkts
      integer,       intent(inout) :: nlev
      character*4,   intent(inout) :: qcs
      logical,       intent(inout) :: daily
      logical,       intent(inout) :: append_mode
      real,          intent(inout) :: lat1
      real,          intent(inout) :: lat2
      real,          intent(inout) :: lon1
      real,          intent(inout) :: lon2
      real,          intent(inout) :: lev1
      real,          intent(inout) :: lev2
      real,          intent(inout) :: levs(*)
      integer,       intent(inout) :: iqcx
      integer,       intent(inout) :: t1
      integer,       intent(inout) :: t2
      integer,       intent(inout) :: miter
      integer,       intent(inout) :: jiter
      logical,       intent(inout) :: ncf
      logical,       intent(inout) :: rmdupl
!
! !OUTPUT PARAMETERS:

      character*255, intent(out) :: infile(nfiles_max)
      integer,       intent(out) :: nfiles
      integer,       intent(out) :: kxs(smax)
      integer,       intent(out) :: kts(smax)
!
!
! !REVISION HISTORY:
!       20Sep98 - D.Dee   Initial code
!       08Dec98 - D.Dee   Added -daily option
!       10Apr00 - D.Dee   Added -lat, -lon options;
!                         modified -kx, -kt, -lev, -o options.
!       22Mar01 - Todling Added capability to handle multiple, 
!                         non-consecutive levels
!       03May01 - Todling Added opt to allow reset of qcexcl flag.
!       30Jan02 - Todling Added opt to select by obs time.
!       13Feb09 - Todling Add -sens option
!     21Jan2014 - Todling Add -prepsigo option
!
!EOP
!BOC

      character*4, parameter :: myname = 'init'

      integer iret, i, lv, iarg, argc, iargc
      integer k, k1, k2, ic, ic1, ic2, lt, ii, ib
      real swap
      character*255 argv
      character*255 SS
      logical osens,adjsigo

!     Parse command line
!     ------------------

      osens = .false.
      adjsigo = .true.
      argc =  iargc()
      if ( argc .lt. 1 ) call usage()
      nfiles = 0
      iarg = 0
      do i = 1, 32767
         iarg = iarg + 1
         if ( iarg .gt. argc ) go to 111
         call GetArg ( iArg, argv )
         if (index(argv,'-o' ) .gt. 0 ) then
            if ( iarg+1 .gt. argc ) call usage()
            iarg = iarg + 1
            call GetArg ( iArg, prefix )
         elseif (index(argv,'-synhour') .gt. 0 ) then
            if ( iarg+1 .gt. argc ) call usage()
            iarg = iarg + 1
            call GetArg ( iArg, SS )
            read(SS,*) synhour
         elseif (index(argv,'-resetqcx') .gt. 0 ) then
            if ( iarg+1 .gt. argc ) call usage()
            iarg = iarg + 1
            call GetArg ( iArg, SS )
            read(SS,*) iqcx
         elseif (index(argv,'-kx') .gt. 0 ) then
            if ( iarg+1 .gt. argc ) call usage()
            iarg = iarg + 1
            call GetArg ( iArg, SS )	    
            ic1 = index(SS,':')
            ic2 = index(SS,',')
            if(ic1/=0.and.ic2/=0) then
                print *, myname, 'invalid kx option kx format'
                stop
            endif
            ic = ic1 + ic2         ! string is k1 or k1:k2 or k1,k2,...
            lt = len_trim(SS)
            if (ic .eq. 0) then    ! no colon, therefore k1
                read(SS,*) k1
                nkxs      = 1
                kxs(nkxs) = k1
            else if(ic1/=0) then   ! colon, therefore k1:k2
                read(SS(1:ic-1) ,*) k1
                read(SS(ic+1:lt),*) k2
                if (nkxs+k2-k1+1 .gt. smax) then
                    print *, myname, 'Too many kx values.'
                    stop
                else
                    kxs(nkxs+1:nkxs+k2-k1+1) = (/ (k, k=k1,k2) /)
                    nkxs = nkxs+k2-k1+1
                endif
            else if(ic2/=0) then   ! commas, therefore k1,k2,...
                ib = 1
                do ii = 1, lt
                   if(SS(ii:ii).eq.',') then
                      nkxs=nkxs+1
                      if (nkxs .gt. smax) then
                          print *, myname, 'Too many kx values.'
                          stop
                      endif
                      read(SS(ib:ii-1),*) k1
                      kxs(nkxs) = k1
                      ib=ii+1
                   endif
                enddo
                nkxs=nkxs+1
                if (nkxs .gt. smax) then
                    print *, myname, 'Too many kx values.'
                    stop
                endif
                read(SS(ib:lt),*) k1
                kxs(nkxs) = k1
            end if
         elseif (index(argv,'-kt') .gt. 0 ) then
            if ( iarg+1 .gt. argc ) call usage()
            iarg = iarg + 1
            call GetArg ( iArg, SS )	    
	    ic = index(SS,':')     ! string is k1 or k1:k2
            lt = len_trim(SS)
            if (ic .eq. 0) then    ! no colon, therefore k1
                read(SS,*) k1
		k2 = k1
            else                   ! colon, therefore k1:k2
                read(SS(1:ic-1) ,*) k1
                read(SS(ic+1:lt),*) k2
            end if	    
            if (nkts+k2-k1+1 .gt. smax) then
                print *, myname, 'Too many kt values.'
                stop
	    else
	        kts(nkts+1:nkts+k2-k1+1) = (/ (k, k=k1,k2) /)
                nkts = nkts+k2-k1+1
            end if
         elseif (index(argv,'-lev') .gt. 0 ) then
            if ( iarg+1 .gt. argc ) call usage()
            iarg = iarg + 1
            call GetArg ( iArg, SS )	    
	    ic = index(SS,':')     ! string is lev1 or lev1:lev2
            lt = len_trim(SS)
            if (ic .eq. 0) then    ! no colon, therefore lev1
                nlev = nlev + 1
                read(SS,*) levs(nlev)
            else                   ! colon, therefore lev1:lev2
                read(SS(1:ic-1) ,*) lev1
                read(SS(ic+1:lt),*) lev2
            end if
	    if (lev2<lev1) then    ! let's be nice to the user..
	        swap = lev2
		lev2 = lev1
		lev1 = swap
	    end if	    
         elseif (index(argv,'-lat') .gt. 0 ) then
            if ( iarg+1 .gt. argc ) call usage()
            iarg = iarg + 1
            call GetArg ( iArg, SS )	    
	    ic = index(SS,':')     ! string is lat1 or lat1:lat2
            lt = len_trim(SS)
            if (ic .eq. 0) then    ! no colon, therefore lat1
                read(SS,*) lat1
		lat2 = lat1
            else                   ! colon, therefore lat1:lat2
                read(SS(1:ic-1) ,*) lat1
                read(SS(ic+1:lt),*) lat2
            end if
	    if (lat2<lat1) then    ! let's be nice to the user..
	        swap = lat2
		lat2 = lat1
		lat1 = swap
	    end if	    
         elseif (index(argv,'-lon') .gt. 0 ) then
            if ( iarg+1 .gt. argc ) call usage()
            iarg = iarg + 1
            call GetArg ( iArg, SS )	    
	    ic = index(SS,':')     ! string is lon1 or lon1:lon2
            lt = len_trim(SS)
            if (ic .eq. 0) then    ! no colon, therefore lon1
                read(SS,*) lon1
		lon2 = lon1
            else                   ! colon, therefore lon1:lon2
                read(SS(1:ic-1) ,*) lon1
                read(SS(ic+1:lt),*) lon2
            end if
	    if (lon2<lon1) then    ! let's be nice to the user..
	        swap = lon2
		lon2 = lon1
		lon1 = swap
	    end if	    
         elseif (index(argv,'-time') .gt. 0 ) then
            if ( iarg+1 .gt. argc ) call usage()
            iarg = iarg + 1
            call GetArg ( iArg, SS )
            ic = index(SS,':')     ! string is t1 or t1:t2
            lt = len_trim(SS)
            if (ic .eq. 0) then    ! no colon, therefore t1
                read(SS,*) t1
                t2 = t1
            else                   ! colon, therefore t1:t2
                read(SS(1:ic-1) ,*) t1
                read(SS(ic+1:lt),*) t2
            end if
            if (t2<t1) then    ! let's be nice to the user..
                swap = t2
                t2   = t1
                t1   = swap
            end if
         elseif (index(argv,'-qc') .gt. 0 ) then
            if ( iarg+1 .gt. argc ) call usage()
            iarg = iarg + 1
            call GetArg ( iArg, SS )
            read(SS,*) qcs
         elseif (index(argv,'-append') .gt. 0 ) then
	    append_mode = .TRUE.
         elseif (index(argv,'-ncf') .gt. 0 ) then
	    ncf = .TRUE.
         elseif (index(argv,'-miter') .gt. 0 ) then   ! can only be invoked w/ -ncf
            if ( iarg+1 .gt. argc ) call usage()
            iarg = iarg + 1
            call GetArg ( iArg, SS )
            read(SS,*) miter
         elseif (index(argv,'-jiter') .gt. 0 ) then   ! can only be invoked w/ -ncf
            if ( iarg+1 .gt. argc ) call usage()
            iarg = iarg + 1
            call GetArg ( iArg, SS )
            read(SS,*) jiter
         elseif (index(argv,'-sens') .gt. 0 ) then
	    osens = .TRUE.
         elseif (index(argv,'-rmdupl') .gt. 0 ) then
	    rmdupl = .TRUE.
         elseif (index(argv,'-daily') .gt. 0 ) then
	    daily = .TRUE.
         elseif (index(argv,'-prepsigo') .gt. 0 ) then
	    adjsigo = .FALSE.
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

!     Consistency check
!     -----------------
      if (ncf) then
          call ods_obsdiags_setparam('ladjsigo',adjsigo)
      endif
      if ( miter>0 .or. jiter>0 ) then
        if (ncf) then
            if (osens) then
               print*, 'Extracting sensitivities and placing them in Xvec slot'
            else
               print*, 'Extracting sensitivities, calculating impact, and placing them in Xvec slot'
            endif
            call ods_obsdiags_init()
            call ods_obsdiags_setparam('lobsdiagsave',.true.)
            call ods_obsdiags_setparam('miter',miter)
            call ods_obsdiags_setparam('jiter',jiter)
            call ods_obsdiags_setparam('lobssens',osens)
        else
            print*, 'Inconsistent choice of command line options'
            print*, 'See usage. Aborting ...'
            call exit(1)
        endif
      endif

!     Echo the parameters
!     -------------------
      if (nkxs .gt. 0) then
          print *
          print *, 'Will select kx:'
          do i = 1, nkxs
              print *, '    ', kxs(i)
          end do
      end if
      if (nkts .gt. 0) then
          print *
          print *, 'Will select kt:'
          do i = 1, nkts
              print *, '    ', kts(i)
          end do
      end if
      if (lat2>=lat1) then
          print *
          print *, 'Will select lat:'
          print *, '     min lat = ', lat1
          print *, '     max lat = ', lat2
      end if
      if (lon2>=lon1) then
          print *
          print *, 'Will select lon:'
          print *, '     min lon = ', lon1
          print *, '     max lon = ', lon2
      end if
      if (lev2>=lev1) then
          print *
          print *, 'Will select lev:'
          print *, '     min lev = ', lev1
          print *, '     max lev = ', lev2
      end if
      if (t2>=t1) then
          print *
          print *, 'Will select within times:'
          print *, '     min time = ', t1
          print *, '     max time = ', t2
      end if
      if (nlev>0) then
          print *
          print *, 'Will select lev:'
          print *, '     min lev = ', levs(1)
          print *, '     max lev = ', levs(nlev)
      end if
      if (qcs .ne. 'ALL') then
          print *
          if (qcs .eq. 'GOOD') then
                           print *, 'Will select ''GOOD'' data only'
          elseif (qcs .eq. 'BAD') then
                           print *, 'Will select ''BAD'' data only'
          end if
      end if
      if (rmdupl) then
          print *
          print *, 'Will eliminate duplicates (kx,kt,time,lat,lon,lev,obs)'
      end if
      print *
      if (synhour .ge. 0)  print *, 'Will process synoptic hour ',
     .                             synhour, ' only'
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

!-------------------------------------------------------------------------

      subroutine usage()
      print *
      print *, 'Usage:'
      print *
      print *, 'odsselect [-kx KX] [-kt KT] [-lon LON1:LON2] ',
     .         '[-lat LAT1:LAT2] [-lev LEV1:LEV2] [-qc QC] ',
     .         '[-time t1:t2] [-ncf] [-rmdupl]',
     .         '[-synhour HH] [-daily] [-o ID] odsfile(s)'
      print *
      print *, 'where'
      print *
      print *,'-kx  KX           select data with kx=KX'
      print *,'-kx  KX1:KX2      select data with kx=KX1,..,KX2'
      print *,'              (default: all kx)'
      print *,'-kt  KT           select data with kt=KT'
      print *,'-kx  KT1:KT2      select data with kx=KT1,..,KT2'
      print *,'              (default: all kt)'
      print *,'-lon LON          select data with lon=LON'
      print *,'-lon LON:LON2     select data with LON1<=lon<=LON2'
      print *,'              (default: all lon)'
      print *,'-lat LAT          select data with lat=LAT'
      print *,'-lat LAT1:LAT2    select data with LAT1<=lat<=LAT2'
      print *,'              (default: all lat)'
      print *,'-lev LEVa LEVb .. select data with lev=LEV'
      print *,'-lev LEV1:LEV2    select data with LEV1<=lev<=LEV2'
      print *,'                 (surface data have lev=2000)'
      print *,'-time TIME        select data at time=TIME'
      print *,'-time T1:T2       select data with T1<=time<=T2'
      print *,'                 (-T1=T2=180)'
      print *,'              (default: all lev)'
      print *,'-qc QC            select data according to the following'
      print *,'                  quality control criteria:'
      print *,'                 -qc GOOD: only data that passed QC'
      print *,'                 -qc BAD:  only data that failed QC'
      print *,'                 -qc ALL:  all data regardless of QC'
      print *,'              (default: QC = GOOD)'
      print *,'-resetqcx FLAG resets all qc-excl flags FLAG to zero'
      print *,'              (default: leave flags alone)'
      print *,'-synhour HH       process synoptic hour HH only'
      print *,'              (default: all)'
      print *,'-daily            create daily output files'
      print *,'              (default: all output in single file)'
      print *,'-append           append to same synoptic time: '
      print *,'                  CAUTION: must specify -synhour'
      print *,'              (default: no append)'
      print *, '-ncf             specify when input files are non-ODS files'
      print *, '                (default: ignore)'
      print *, '-miter           specify number of iterations in GSI outer loop'
      print *, '                (default: miter=2, when -ncf opt is chosen)'
      print *, '-jiter           specify GSI inner iteration number'
      print *, '                (default: jiter=1, when -ncf opt is chosen)'
      print *, '-sens            extract observation sensitivity only (no impact)'
      print *, '                (default: ignore)'
      print *, '-rmdupl          remove duplicates (kx,kt,time,lat,lon,lev,obs)'
      print *, '                (default: do not)'
      print *, '-prepsigo       extract prep-sigo instead of adjusted sigo from diag files'
      print *, '                (default: adjusted sigo)'
      print *,'-o ID             use ID for naming output files'
      print *,'              (default: SELECT)'
      print *,' odsfile(s)       ODS file(s)'
      print *
      print *, 'You may include multiple -kx, -kt and -lev flags, as in'
      print *
      print *, ' odsselect -kt 1 -kt 4:6 test.ods'
      print *, ' and'
      print *, ' odsselect -lev 500 300 -lev 10:100 test.ods'
      print *
      print *, ' CAUTION: '
      print *, '   -jiter, -miter, osens opts can only be used in conjunction with -ncf'
      print *
      stop
      end
