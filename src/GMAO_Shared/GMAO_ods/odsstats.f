!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!-----------------------------------------------------------------------
!BOP
!
! !ROUTINE: odsstats:  Operates on specific attributes in ODS
!
! !INTERFACE:

      program odsstats

!
! !DESCRIPTION:
!           Reads ODS files, selects data according to criteria specified
!           via a resource file, operates on chosen attribute and summarizes
!           results by writing them to an ASC file.
!

      use m_odsmeta, only : KTMAX, KXMAX
      use m_ods
      use m_inpak90
      use m_chars,   only: lowercase
      use m_stdio,   only: stdout, stderr

      implicit NONE

!
! !REVISION HISTORY:
!
!     23Jan2008 Todling  Initial code
!     23Sep2008 Todling  Add stats summary and skip verbose
!     01Sep2009 Todling  Add obvar to handle attr other than xvec
!     13Nov2009 Todling  Accommulate global rms
!     13Apr2010 Gelaro   Accumulate numneg, increase NCMAX LVMAX for IASI 
!     15Aug2013 Todling/Daescu  Implement calc of imp from sensitvity
!     19Feb2014 Todling/Daescu  Implement calc of oma * xvec 
!     24Feb2014 Todling  Revisit Write_AccumStats (now to file and/or screen)
!
!EOP
!BOC

      character(len=*), parameter :: myname = 'odsstats'

      integer, parameter :: NFILES_MAX = 150 ! Max number of input files

      integer, parameter ::  NOMAX = 5       ! max no. of operations
      integer, parameter ::  NCMAX = 800     ! max no. of classes
      integer, parameter ::  NTMAX = 31 * 4  ! 4 syn times a day, for a month
      integer, parameter ::  LVMAX = 800     ! Max number of levels/channels

      logical, parameter :: debug = .false.
      logical            :: verb  = .false.

!     Local variables
!     ---------------
      integer nfiles, ifile, lf, isyn, ksyn, nymd, nhms, nkt, nkx, nlv, ncfound
      integer i, j, ic, nc, nt, ierr, nobs, nsel, synhour, nop, nops
      integer nymdb, nhmsb

      character(len=255) :: opers(nomax)          ! operations
      character(len=255) :: oclass(ncmax)         ! observations classes
      character(len=255) :: obvar                 ! get this attribute from ODS (usually, xvec)     
      character(len=255) :: RCfile
      integer            :: kttype(ktmax,ncmax)   !
      integer            :: kxtype(kxmax,ncmax)   !
      real               :: levlst(lvmax,ncmax)   ! levels/channels of selected obs
      real, pointer      :: ptr(:,:)              ! number of obs giving positive attribute
      real, target       :: obssum(ncmax,ntmax)   ! sum per class and syn hour
      real, target       :: negsum(ncmax,ntmax)   ! number of obs for which attr are negative
      real, target       :: obsgms(ncmax,ntmax)   ! global mean square of obs for each attr
      integer            :: obsnum(ncmax,ntmax)   ! number of obs per class and syn hour
      integer            :: nymda(ntmax)          ! all dates
      integer            :: nhmsa(ntmax)          ! all times
      integer            :: trange(2)             ! time interval of selected obs
      real               :: latrange(2)           ! latitudes  of selected obs
      real               :: lonrange(2)           ! longitudes of selected obs
      real               :: accum_obssum(ncmax)   ! accumulated number of obs per class
      real               :: accum_obsgms(ncmax)   ! accumulated global mean square of obs per class
      integer            :: accum_obsnum(ncmax)   ! accumulated impact per class
      integer            :: accum_obsneg(ncmax)   ! accumulated number of obs with neg attribute
 
      character*255 infile (NFILES_MAX) ! input filenames
      character*255 fileout, outfile, outodsfn   ! output filename
      character*80  ftype
      logical       anotherclass, allkxs, allevs, lstdv

!     storage for ODS:
!     ---------------
      type ( ods_vect ) ods
      type ( ods_vect ) odss

!     Option flags:
!     ------------
      synhour = -1       ! DEFAULT: process all synoptic times

      oclass = 'NONE'
      kxtype = 0
      kttype = 0
      levlst = -1.0
      lstdv = .false.

!     Parse in command line
!     ---------------------
      call init ( infile, nfiles_max, nfiles, RCfile, outfile, trange, 
     .            latrange, lonrange, verb )

!     Read in resource file
!     ---------------------
      call ObsOperSet_ ( RCfile, ierr )

!     Figure out range of synoptic time loop
!     --------------------------------------
      ksyn = 32767
      nt   = 0
      accum_obsnum = 0
      accum_obssum = 0.0
      accum_obsgms = 0.0
      accum_obsneg = 0
      obssum = 0.0
      negsum = 0.0
      obsgms = 0.0
      nymdb = -1
      nhmsb = -1

!     Loop over input files
!     ---------------------
      do ifile = 1, nfiles

!       Loop over all synoptic times on this file
!       -----------------------------------------
        do isyn = 1, ksyn

!         Read all data for this synoptic time
!         ------------------------------------

	  nymd = -1            ! get data for next synoptic time on file
	  nhms =  0
          call ODSNxTime ( trim(infile(ifile)), nymd, nhms )
          if ( nymd .eq. -1 ) then 
               if(verb) print *, 'End-Of-File'
               exit     ! skip to next file
          end if
	  
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

          if ( nymdb<0.and.nhmsb<0 ) then
               nymdb=nymd; nhmsb=nhms
          endif
          if (verb) then
              print *, 'number of obs read = ', nobs
              print *, 'date = ', nymd
              print *, 'time = ', nhms
          endif

!         Loop over observation classes
!         -----------------------------
          do nc = 1, ncfound
             if(verb) print *, trim(oclass(nc))

!            Select observations need for operating with present class
!            ---------------------------------------------------------
             call icount_this_ ( kttype(:,nc), nkt )
             if ( allkxs .and. allevs ) then
                  call ODS_Select ( ods, nobs, nsel, ierr,
     .                              odss=odss,
     .                              qcexcl=0, 
     .                                        kt_list=kttype(1:nkt,nc),
     .                                        time_range=trange,
     .                                         lat_range=latrange,
     .                                         lon_range=lonrange )
             else if ( allevs .and. (.not. allkxs ) ) then
                  call icount_this_ ( kxtype(:,nc), nkx )
                  call ODS_Select ( ods, nobs, nsel, ierr,
     .                              odss=odss,
     .                              qcexcl=0, kx_list=kxtype(1:nkx,nc),
     .                                        kt_list=kttype(1:nkt,nc),
     .                                        time_range=trange,
     .                                         lat_range=latrange,
     .                                         lon_range=lonrange )
             else 
                  call icount_this_ ( kxtype(:,nc), nkx )
                  call rcount_this_ ( levlst(:,nc), nlv )
                  call ODS_Select ( ods, nobs, nsel, ierr,
     .                              odss=odss,
     .                              qcexcl=0,  kx_list=kxtype(1:nkx,nc),
     .                                         kt_list=kttype(1:nkt,nc),
     .                                        lev_list=levlst(1:nlv,nc),
     .                                        time_range=trange,
     .                                         lat_range=latrange,
     .                                         lon_range=lonrange )
             endif

             obsnum(nc,nt) = nsel
             nymda(nt)     = nymd
             nhmsa(nt)     = nhms
             do nop = 1, nops
                if(opers(nop)=='sum')    call sumattr ( verb, odss, nsel, obvar, obssum(nc,nt) )
                if(opers(nop)=='numneg') call negattr ( verb, odss, nsel, obvar, negsum(nc,nt) )
                if(opers(nop)=='stdev')  then
                   call gmsattr ( verb, odss, nsel, obvar, obsgms(nc,nt) )
                   lstdv = .true.
                 endif
             enddo

!            Accumulate obs counts and impacts
!            ---------------------------------
             accum_obsnum(nc) = accum_obsnum(nc) + obsnum(nc,nt)
             accum_obssum(nc) = accum_obssum(nc) + obssum(nc,nt)
             accum_obsgms(nc) = accum_obsgms(nc) + obsgms(nc,nt)
             accum_obsneg(nc) = accum_obsneg(nc) + negsum(nc,nt)

!            When debugging, write out select ODS entries
!            ---------------------------------------------
             if ( debug ) then
                  write(outodsfn,'(2a,i8.8,a,i2.2,a)') trim(oclass(nc)), '.obs.', nymd, '_', nhms/10000, 'z.ods'
                  print *, 'nsel, nobs', nsel, odss%data%nobs
                  print *, 'Writing ods file for debug purposes: ', trim(outodsfn)
                  call ODS_Put ( trim(outodsfn), ftype, nymd, nhms, odss, ierr, append=.true. )
             endif

             call ODS_Clean ( odss, ierr ) 
          end do

!         Write to ASC file
!         -----------------
             if(verb) print *, 'calling Write_Stats'
          do nop = 1, nops
             if(opers(nop)=='sum') then
                ptr => obssum
             endif
             if(opers(nop)=='numneg') then
                ptr => negsum
             endif
             if(opers(nop)=='stdev') then
                ptr => obsgms
             endif
             fileout = trim(outfile) // '_' // trim(opers(nop)) // '.txt'
             call Write_Stats ( fileout, ptr, obsnum, oclass, ncfound, ncmax, nt, nymda, nhmsa, verb, ierr ) 
          enddo
             if(verb) print *, 'completed Write_Stats'

          call ODS_Clean ( ods, ierr )
             if ( ierr .ne. 0 ) then
                  print *, 'ODS_Clean error: ierr = ', ierr
             endif

        end do  ! loop over synoptic times

      end do  ! loop over files

!     Write out accumulated results
!     ------------------------------
      fileout = trim(outfile) // '_' // 'all.txt'
      call Write_AccumStats ( fileout, accum_obssum, accum_obsgms, accum_obsnum, accum_obsneg, 
     .                        oclass, ncfound, ncmax, nymdb, nhmsb, lstdv, verb, ierr ) 

      CONTAINS

!-------------------------------------------------------------------------
!         NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!-------------------------------------------------------------------------
!BOP
!
! !ROUTINE: ObsOpertSet --- General paramater setup for singular vector calc
! 
! !INTERFACE:
!
      subroutine ObsOperSet_ ( RCfile, stat )
 
! !USES: 
    
      Implicit None

! !INPUT PARAMETERS: 
!
      character(len=*), intent(in) :: RCfile  ! resource filename

! !OUTPUT PARAMETERS:

      integer,          intent(out) :: stat                  ! return error code


! !DESCRIPTION:  Initialize observation operations program.
!
! !REVISION HISTORY: 
!
!   23Jan2008  Todling    Initial code.
!   05Aug2008  Todling    Add table of levs/channels.
!
!EOP
!-------------------------------------------------------------------------

      character(len=*), parameter :: myname_ = myname//'ObsOperSet_'

      character(len=255) token, obsoperrc, tablename

      integer irow, j, n, kt, kx, jcnt, maxn, iret
      integer ii, lt, k1, k2, knext
      real    r1, r2, rnext

      stat = 0

!     Load resources file
!     -------------------
      obsoperrc = ' '
      call getenv('obs_opers.rc',OBSOPERRC)     ! Unix binding
      if(obsoperrc.eq.' ') obsoperrc=RCfile     ! default name
      call i90_loadf (trim(obsoperrc), iret)
      if( iret .ne. 0) then
          write(stderr,'(2a,i5)') myname_,': I90_loadf error, iret =',iret
          stat = 1
          return
      end if
      write(stdout,'(1x,  a  )') '---------------------------------------------------'
      write(stdout,'(1x, 2a  )') myname_, ': Reading resource file'
      write(stdout,'(1x,  a,/)') '---------------------------------------------------'

!     Define type of operations
!     -------------------------
      obvar = 'xvec'
      call I90_label('OBS*Variable:', iret)
      if (iret .eq. 0) then
          call I90_Gtoken ( token, iret )
          if ( iret==0 ) then
               obvar = lowercase(trim(token))
          endif
      end if
      write(stdout,'(2a)') ' Will extract information from this attribute: ', trim(obvar)

!     Define type of operations
!     -------------------------
      opers(1:nomax) = 'NONE'
      call I90_label('OBS*Operations:', iret)
      if (iret .eq. 0) then
          nops = 0
          do while ( iret==0 .and. jcnt<nomax )
             call I90_Gtoken ( token, iret )
             if ( iret==0 ) then
                  nops = nops + 1
                  opers(nops) = lowercase(trim(token))
                  write(stdout,'(2a)') ' Will perform the following operations: ', trim(opers(nops))
             endif
          enddo
      else
         nops      = 1
         opers(1)  = 'sum'   ! default: sum residuals
                             ! other possibilities are:
                             !  ave - average residuals
      end if

!     Read table with variable types
!     ------------------------------
      tablename = 'OBS*Class*Kts::'
      call I90_label(trim(tablename), iret)
      if (iret/=0) then
         write(stderr,'(2a,i5,2a)') myname_, ': I90_label error, iret=', iret,
     .                                       ': trying to read ', trim(tablename)
         stat = 2; return
      end if
      irow = 0
      ncfound = 0
      write(stdout,'(a)') ' Will process the following obs classes: '
      do while (iret==0)                   ! read table entries
         call I90_GLine ( iret )           ! iret=-1: end of file; +1: end of table
         if (iret==0.and.irow<ncmax) then  ! OK, we have next row of table
             irow = irow + 1

             call I90_GToken(token, iret ) ! obs class name
             if (iret/=0) then
                 write(stderr,'(2a,i5)') myname_, ': I90_GToken error, iret=', iret
             end if
             oclass(irow) = trim(token)
             write(stdout,'(1x,a)') trim(oclass(irow))

             jcnt=0
             ierr=0
             do  j = 1, ktmax
               call I90_GToken(token, ierr )
               if(ierr/=0) exit
               ii = index(token,':') ! token is single entry or range of entries
               lt = len_trim(token)
               if (ii==0) then       ! no colon, therefore single entry
                   read(token,*) knext
                   jcnt = jcnt + 1
                   kttype(jcnt,irow) = knext
               else                  ! colon, therefore kx1:kx2
                   read(token(1:ii-1),*) k1
                   read(token(ii+1:lt),*) k2
                   do knext = k1, k2
                      if (jcnt==KTMAX) then    ! check space
                        write(stderr,'(2a,i5)') myname,': increase KTMAX'
                        stat = 4; return
                      else if (ierr==0) then
                        jcnt = jcnt + 1
                        kttype(jcnt,irow) = knext
                      end if
                   end do
               end if
             enddo

             ncfound = max(ncfound,irow)
         end if
      end do
      nkt = jcnt
      print *, 'Will process the following Kt''s:'
      do j = 1, ncfound
         call icount_this_ ( kttype(:,j), n )
         print *, kttype(1:n,j)
      enddo

!     Read table with instruments types
!     ---------------------------------
      allkxs = .false.
      tablename = 'OBS*Class*Kxs::'
      call I90_label(trim(tablename), iret)
      if (iret/=0) then
         write(stderr,'(4a)') myname_, ': table ', trim(tablename),
     .                        'not found in RC file ... will take all KXs'
         allkxs = .true.
      end if
      if ( .not. allkxs ) then
        irow = 0
        do while (iret==0)                     ! read table entries
           call I90_GLine ( iret )             ! iret=-1: end of file; +1: end of table
           if (iret==0.and.irow<ncmax) then    ! OK, we have next row of table
               irow = irow + 1
  
               call I90_GToken ( token, iret ) ! obs class name
               if (iret/=0) then
                   write(stderr,'(2a,i5)') myname_, ': I90_GToken error, iret=', iret
               end if
               if(trim(token)/=trim(oclass(irow))) then     ! tables of kx and kt must be ordered 
                                                            ! in the same way, otherwise drop it
                  stat =4; return
               endif
  
               jcnt=0
               ierr=0
               do  j = 1, kxmax
                 call I90_GToken(token, ierr )
                 if(ierr/=0) exit
                 ii = index(token,':') ! token is single entry or range of entries
                 lt = len_trim(token)
                 if (ii==0) then       ! no colon, therefore single entry
                     read(token,*) knext
                     jcnt = jcnt + 1
                     kxtype(jcnt,irow) = knext
                 else                  ! colon, therefore kx1:kx2
                     read(token(1:ii-1),*) k1
                     read(token(ii+1:lt),*) k2
                     do knext = k1, k2
                        if (jcnt==KXMAX) then    ! check space
                          write(stderr,'(2a,i5)') myname,': increase KXMAX'
                          stat = 4; return
                        else if (ierr==0) then
                          jcnt = jcnt + 1
                          kxtype(jcnt,irow) = knext
                        end if
                     end do
                 end if
               enddo

           end if
        end do
        nkx = jcnt
        print *, 'Will process the following Kx''s:'
        do j = 1, ncfound
           call icount_this_ ( kxtype(:,j), n )
           print *, kxtype(1:n,j)
        enddo
      endif ! < all KXs >

!     Read table with levels/channels
!     -------------------------------
      allevs = .false.
      tablename = 'OBS*Class*Lev::'
      call I90_label(trim(tablename), iret)
      if (iret/=0) then
         write(stderr,'(4a)') myname_, ': table ', trim(tablename),
     .                        'not found in RC file ... will take all levels/channels'
         allevs = .true.
      end if
      if ( .not. allevs ) then
        irow = 0
        do while (iret==0)                     ! read table entries
           call I90_GLine ( iret )             ! iret=-1: end of file; +1: end of table
           if (iret==0.and.irow<ncmax) then    ! OK, we have next row of table
               irow = irow + 1
  
               call I90_GToken ( token, iret ) ! obs class name
               if (iret/=0) then
                   write(stderr,'(2a,i5)') myname_, ': I90_GToken error, iret=', iret
               end if
               if(trim(token)/=trim(oclass(irow))) then     ! tables of kx and kt must be ordered 
                                                            ! in the same way, otherwise drop it
                  stat =5; return
               endif
  
               jcnt=0
               ierr=0
               do  j = 1, lvmax
                 call I90_GToken(token, ierr )
                 if(ierr/=0) exit
                 ii = index(token,':') ! token is single entry or range of entries
                 lt = len_trim(token)
                 if (ii==0) then       ! no colon, therefore single entry
                     read(token,*) rnext
                     jcnt = jcnt + 1
                     levlst(jcnt,irow) = rnext
                 else                  ! colon, therefore kx1:kx2
                     read(token(1:ii-1),*) r1
                     read(token(ii+1:lt),*) r2
                     k1=nint(r1); k2=nint(r2)    ! levs/chanels converted to integer
                     do knext = k1, k2
                        if (jcnt==LVMAX) then    ! check space
                          write(stderr,'(2a,i5)') myname,': increase LVMAX'
                          stat = 4; return
                        else if (ierr==0) then
                          jcnt = jcnt + 1
                          levlst(jcnt,irow) = knext
                        end if
                     end do
                 end if
               enddo

           end if
        end do
        nlv = jcnt
        print *, 'Will process the following lev/chanel''s:'
        do j = 1, ncfound
           call rcount_this_ ( levlst(:,j), n )
           print *, levlst(1:n,j)
        enddo
      endif ! < all Levs >

!     release resource file:
!     ---------------------
      call I90_release()

      return
      end subroutine ObsOperSet_

      subroutine icount_this_ ( indx, n )
      implicit none
      integer, intent(in)  :: indx(:)
      integer, intent(out) :: n
      integer  i,m
      m=size(indx,1)
      n=0
      do i = 1, m
         if(indx(i)/=0) n=n+1 
      enddo
      end subroutine icount_this_

      subroutine rcount_this_ ( array, n )
      implicit none
      real,    intent(in)  :: array(:)
      integer, intent(out) :: n
      integer  i,m
      m=size(array,1)
      n=0
      do i = 1, m
         if(array(i)>0.0) n=n+1 
      enddo
      end subroutine rcount_this_


!EOC

      end program odsstats ! program odsstats

      subroutine sumattr ( verb, ods, nobs, attr, sum )
      use m_ods
      implicit none
      type(ods_vect) ods
      integer, intent(in) :: nobs
      character(len=*), intent(in) :: attr 
      real,    intent(inout) ::  sum
      logical, intent(in) :: verb
      integer  i

      sum = 0.0
      if (nobs==0) return ! nothing to do
      if ( trim(attr) == 'obs' ) then
         do i=1,nobs
            sum=sum+ods%data%obs(i)
         enddo
      endif
      if ( trim(attr) == 'omf' ) then
         do i=1,nobs
            sum=sum+ods%data%omf(i)
         enddo
      endif
      if ( trim(attr) == 'oma' ) then
         do i=1,nobs
            sum=sum+ods%data%oma(i)
         enddo
      endif
      if ( trim(attr) == 'xm' ) then
         do i=1,nobs
            sum=sum+ods%data%xm(i)
         enddo
      endif
      if ( trim(attr) == 'xvec' ) then
         do i=1,nobs
            sum=sum+ods%data%xvec(i)
         enddo
      endif
      if ( trim(attr) == 'imp_from_sens' .or.
     .     trim(attr) == 'xvecxomf'    ) then  ! this handles the case when xvec holds sensitivities
                                               ! rather than impacts themselves
         do i=1,nobs
            sum=sum+ods%data%xvec(i)*ods%data%omf(i) ! calculate impact
         enddo
      endif
      if ( trim(attr) == 'xvecxoma' ) then  ! this handles the case when xvec holds sensitivities
                                            ! rather than impacts themselves
         do i=1,nobs
            sum=sum+ods%data%xvec(i)*ods%data%oma(i) ! calculate impact
         enddo
      endif
      if ( verb ) then
          print *, 'Obs impact(sum):      ', sum
          print *, 'Obs impact(sum)/Nobs: ', sum/nobs
      endif

      end subroutine sumattr

      subroutine negattr ( verb, ods, nobs, attr, sum )
      use m_ods
      implicit none
      type(ods_vect) ods
      integer, intent(in) :: nobs
      character(len=*), intent(in) :: attr 
      real,    intent(inout) ::  sum
      logical, intent(in) :: verb
      integer  i

      sum = 0.0
      if (nobs==0) return ! nothing to do
      if ( trim(attr) == 'xvec' ) then
         do i=1,nobs
            if(ods%data%xvec(i)<0.0) sum=sum+1.0
         enddo
      endif
      if ( trim(attr) == 'imp_from_sens' .or.
     .     trim(attr) == 'xvecxomf'    ) then  ! this handles the case when xvec holds sensitivities
         do i=1,nobs
            if(ods%data%xvec(i)*ods%data%omf(i)<0.0) sum=sum+1.0
         enddo
      endif
      if ( verb ) then
          print *, '    Number of obs with positive impact:  ',  sum
          print *, 'Percentage of obs with positive impact:  ', (sum/nobs)*100.
      endif

      end subroutine negattr

      subroutine gmsattr ( verb, ods, nobs, attr, rsum )
      use m_ods
      implicit none
      type(ods_vect) ods
      integer, intent(in) :: nobs
      character(len=*), intent(in) :: attr 
      real,    intent(inout) ::  rsum
      logical, intent(in) :: verb
      integer  i

      rsum = 0.0 
      if (nobs==0) return ! nothing to do
      if ( trim(attr) == 'obs' ) then
         do i=1,nobs
            rsum=rsum+(ods%data%obs(i))**2
         enddo
      endif
      if ( trim(attr) == 'omf' ) then
         do i=1,nobs
            rsum=rsum+(ods%data%omf(i))**2
         enddo
      endif
      if ( trim(attr) == 'oma' ) then
         do i=1,nobs
            rsum=rsum+(ods%data%oma(i))**2
         enddo
      endif
      if ( trim(attr) == 'xm' ) then
         do i=1,nobs
            rsum=rsum+(ods%data%xm(i))**2
         enddo
      endif
      if ( trim(attr) == 'xvec' ) then
         do i=1,nobs
!           rsum=rsum+(ods%data%xvec(i))**2
            rsum=rsum+(max(1.e-15,abs(ods%data%xvec(i))))**2 ! hack: "-fp-model precise" flag
                                                             ! of intel compiler version 11 gives
                                                             ! underflow here in some cases - non-sense
         enddo
      endif
      if ( verb ) then
          print *, 'MS of impacts (', trim(attr), '): ', rsum
      endif

      end subroutine gmsattr

      subroutine Write_Stats ( outfile, obssum, obsnum, oclass, ncfound, ncmax, nt, nymd, nhms, verb, stat ) 
      use m_ioutil, only : luavail
      implicit none
      integer, intent(in)          :: ncfound, ncmax, nt
      integer, intent(in)          :: nymd(nt), nhms(nt)
      character(len=*), intent(in) :: outfile
      character(len=*), intent(in) :: oclass(ncfound)
      real,    intent(in)          :: obssum(ncmax,nt)
      integer, intent(in)          :: obsnum(ncmax,nt)
      logical, intent(in)          :: verb
      integer, intent(out)         :: stat

      integer i,j,ios,lu

      stat = 0
      lu=luavail()
      open(lu, file=outfile, form='formatted', iostat=ios)
      if ( ios .ne. 0 ) then
           print *, 'Cannot open file ', trim(outfile)
           stat = 1
           return
      else
           if(verb) print *, 'Writing ASC file ', trim(outfile)
      end if

      do j = 1, nt ! loop over times
         do i = 1, ncfound  ! loop ob classes
                            ! intentionally truncate oclass to 20 chars
            write(lu,'(i8.8,1x,i6.6,1x,a20,1x,i11,1x,1p,e11.4)') 
     .            nymd(j), nhms(j), oclass(i), obsnum(i,j), obssum(i,j)
         enddo
      enddo
      close(lu)

      end subroutine Write_Stats

      subroutine Write_AccumStats ( fname, obssum, obsgms, obsnum, obsneg, oclass, ncfound, ncmax, nymd, nhms, lstdev, 
     .                              verb, stat )
      use m_ioutil, only : luavail
      implicit none
      character(len=*), intent(in) :: fname
      integer, intent(in)          :: ncfound, ncmax
      character(len=*), intent(in) :: oclass(ncfound)
      real,    intent(in)          :: obssum(ncmax)
      real,    intent(in)          :: obsgms(ncmax)
      integer, intent(in)          :: obsnum(ncmax)
      integer, intent(in)          :: obsneg(ncmax)
      integer, intent(in)          :: nymd, nhms
      logical, intent(in)          :: lstdev
      logical, intent(in)          :: verb
      integer, intent(out)         :: stat

      integer,allocatable :: lu(:)
      integer i,j,ios,ll,n,nu
      real    xmean,rmean,tmp,xstdv

      stat = 0
      nu = 1
      if(verb) nu=2
      allocate(lu(nu))
      lu(1)=luavail()
      if(verb) lu(2)=6
      open(lu(1), file=fname, form='formatted', iostat=ios)

!     Write out accumulated results
!     -----------------------------
      if (lstdev) then
         do ll=1,nu
            write(lu(ll),'(44x,4a)') '#obs ', ' sum_impact ', '     numneg ', '    stdev'
         enddo
         do i = 1, ncfound  ! loop ob classes
            xstdv = 0.0
            n = obsnum(i)
            if (n>1) then
                xmean = obssum(i)/n
                rmean = obsgms(i)/n
                tmp   = abs(rmean-xmean**2)
                xstdv = sqrt(n*tmp/(n-1))
            endif
                            ! intentionally truncate oclass to 20 chars
            do ll=1,nu
            write(lu(ll),'(i8.8,1x,i6.6,1x,a20,1x,i11,1x,1p,e11.4,1x,i11,1x,1p,e11.4)') nymd, nhms, 
     .                                       oclass(i), obsnum(i), obssum(i), obsneg(i), xstdv
            enddo
         enddo
      else
         do ll=1,nu
            write(lu(ll),'(44x,3a)') '#obs ', ' sum_impact ', '     numneg '
            do i = 1, ncfound  ! loop ob classes
                               ! intentionally truncate oclass to 20 chars
               write(lu(ll),'(i8.8,1x,i6.6,1x,a20,1x,i11,1x,1p,e11.4,1x,i11)') 
     .               nymd, nhms, oclass(i), obsnum(i), obssum(i), obsneg(i)
         enddo
      enddo
      endif

      close(lu(1))
      deallocate(lu)

      end subroutine Write_AccumStats


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
      subroutine init ( infile, nfiles_max, nfiles, RCfile, outfile, 
     .                  trange, latrange, lonrange, verb )

! !USES:

! !INPUT PARAMETERS:
!
      implicit NONE
      integer,       intent(in)  :: nfiles_max

! !OUTPUT PARAMETERS:

      character(len=*), intent(out) :: RCfile
      character(len=*), intent(out) :: infile(nfiles_max)
      integer,          intent(out) :: nfiles
      character(len=*), intent(out) :: outfile
      integer,          intent(out) :: trange(2)
      real,             intent(out) :: latrange(2)
      real,             intent(out) :: lonrange(2)
      logical,          intent(out) :: verb
!
!
! !REVISION HISTORY:
!     22Jan2008 Todling - Initial code (stripped off odsselect)
!     19Feb2008 Todling - Add trange
!     16Apr2008 Todling - Add lat/lon ranges
!
!EOP
!BOC

      character*4, parameter :: myname_ = 'init'

      integer iret, i, ic, lt, lv, iarg, argc, iargc
      real swap
      character*255 argv
      character*10 SS

      RCfile  = 'odsstats.rc'
      outfile = 'odsstats'
      trange  = (/-9999,+9999/)   ! Default: includes obs in all time ranges
      latrange  = (/-90.,+90./)   ! Default: includes obs in all latitude ranges
      lonrange  = (/-180.,+180./) ! Default: includes obs in all longitude ranges 
      verb = .false.

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
         if (index(argv,'-o' ) .gt. 0 ) then
            if ( iarg+1 .gt. argc ) call usage()
            iarg = iarg + 1
            call GetArg ( iArg, outfile )
         else if (index(argv,'-verbose' ) .gt. 0 ) then
            verb = .true.
         else if (index(argv,'-rc' ) .gt. 0 ) then
            if ( iarg+1 .gt. argc ) call usage()
            iarg = iarg + 1
            call GetArg ( iArg, RCfile )
         elseif (index(argv,'-time') .gt. 0 ) then
            if ( iarg+1 .gt. argc ) call usage()
            iarg = iarg + 1
            call GetArg ( iArg, SS )
            ic = index(SS,':')     ! string is t1 or t1:t2
            lt = len_trim(SS)
            if (ic .eq. 0) then    ! no colon, therefore t1
                read(SS,*) trange(1)
                trange(2) = trange(1)
            else                   ! colon, therefore t1:t2
                read(SS(1:ic-1) ,*) trange(1)
                read(SS(ic+1:lt),*) trange(2)
            end if
            if (trange(2)<trange(1)) then    ! let's be nice to the user..
                swap = trange(2)
                trange(2) = trange(1)
                trange(1) = swap
            end if
         elseif (index(argv,'-lat') .gt. 0 ) then
            if ( iarg+1 .gt. argc ) call usage()
            iarg = iarg + 1
            call GetArg ( iArg, SS )
            ic = index(SS,':')     ! string is t1 or t1:t2
            lt = len_trim(SS)
            if (ic .eq. 0) then    ! no colon, therefore t1
                read(SS,*) latrange(1)
                latrange(2) = latrange(1)
            else                   ! colon, therefore t1:t2
                read(SS(1:ic-1) ,*) latrange(1)
                read(SS(ic+1:lt),*) latrange(2)
            end if
            if (latrange(2)<latrange(1)) then    ! let's be nice to the user..
                swap = latrange(2)
                latrange(2) = latrange(1)
                latrange(1) = swap
            end if
         elseif (index(argv,'-lon') .gt. 0 ) then
            if ( iarg+1 .gt. argc ) call usage()
            iarg = iarg + 1
            call GetArg ( iArg, SS )
            ic = index(SS,':')     ! string is t1 or t1:t2
            lt = len_trim(SS)
            if (ic .eq. 0) then    ! no colon, therefore t1
                read(SS,*) lonrange(1)
                lonrange(2) = lonrange(1)
            else                   ! colon, therefore t1:t2
                read(SS(1:ic-1) ,*) lonrange(1)
                read(SS(ic+1:lt),*) lonrange(2)
            end if
            if (lonrange(2)<lonrange(1)) then    ! let's be nice to the user..
                swap = lonrange(2)
                lonrange(2) = lonrange(1)
                lonrange(1) = swap
            end if
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

      print *
      print *, 'Input files: ', nfiles
      print *
      do i = 1, nfiles
         lv = len_trim(infile(i))
         print *, ' o ', infile(i)(1:lv)
      end do
      print *
      print *, 'Output filename: ', trim(outfile)

      return

      end ! subroutine init

!EOC

!-------------------------------------------------------------------------

      subroutine usage()
      print *
      print *, 'Usage:'
      print *
      print *, 'obsstats [-o ID] -rc RCfile odsfile(s)'
      print *
      print *, 'where'
      print *
      print *,'-o  ID         use ID for naming output files'
      print *,'                (default: obsstat.txt)'
      print *,'-verbose       sets verbose on (default: off)'
      print *,'-rc RCfile      resource file'
      print *,'                (default: obsstat.rc)'
      print *,'-time ti:tf    specify window of time to select obs from'
      print *,'                (default: -999:+999, all)'
      print *,'-lat  latA:latB specify latitude range'
      print *,'                (default: -90:+90, all)'
      print *,'-lon  lonA:lonB specify longitude range'
      print *,'                (default: -180:+180, all)'
      print *
      print *
      print *,' odsfile(s)    ODS file(s)'
      print *
      print *, 'NOTES: '
      print *, '-----  '
      print *
      print *, ' The following are entries allowed in the rc file'
      print *
      print *, '  1. Known OBS*Operations: (default: sum)'
      print *, '     sum     - sum all variables '
      print *, '     stddev  - stdandard deviation from mean of variables '
      print *, '     numneg  - sum all non-negative variables '
      print *
      print *, '  2. Known OBS*Variable: (default: xvec)'
      print *, '     obs,omf,oma,xm,xvec,xvecXomf,xvecXoma '
      print *, '     2a. xvecXomf,xvecXoma are only meaningful when xvec holds '
      print *, '         the sensitivities instead of the impacts.'
      print *
      stop
      end
