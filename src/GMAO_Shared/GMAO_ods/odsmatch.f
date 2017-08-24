
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!-----------------------------------------------------------------------
!BOP
!
! !ROUTINE: odsmatch:  Finds matching entries in two ods files.
!
! !INTERFACE:

      program odsmatch

!
! !USAGE: see the routine usage() below
!
! !DESCRIPTION:
!

! !USES:

      use m_MergeSorts
      use m_ods
      use m_odsmeta, only : KXMAX
      use m_odsmeta, only : KTMAX
      use m_odsmeta, only : X_ODSMATCH
      use m_inpak90, only : i90_loadf
      use m_inpak90, only : i90_label
      use m_inpak90, only : i90_Gline
      use m_inpak90, only : i90_Gtoken
      use m_inpak90, only : i90_release

      implicit NONE

!
! !REVISION HISTORY:
!
!     02May2005 Dee     - original code
!     22Jan2008 Todling - add obs impacts to xvec
!     11Aug2008 Redder  - Enabled program to handle data using height assigments
!     29Oct2008 Todling - Placed wired-in kx of Redder change into rc file
!     15Feb2009 Todling - Add opt to handle either sensitivy or impact
!                       - Add opt to generate sigmaO impacts
!     25Feb2009 Todling - Add ability to read range of kx's from rc file
!     31Jan2010 Todling - No longer uses ks to match 
!     23Feb2010 Todling - Add imp0hr (merge with das215)
!     21Apr2013 Todling - Add DFS (Lupu et al. 2011)
!     15Aug2013 Daescu/RT - Sigo impact had opposite sign
!
! !REMARKS:
!   1) I am temporarily adding the energy norm routine here, though I will
!   likely removed from here in the near future - this is left here for now not
!   to break my settings. RT
!
!EOP
!BOC

      character*7, parameter :: myname = 'odsmatch'


!     Local variables
!     ---------------
      integer, parameter :: MAXHOBS = 200  ! max no. of kx's using height for lev assignment 
      integer i, j, ii, jj, jt, k, nx, ny, m, nh
      integer ierr, isyn, nymd, nhms, synhour
      real    d
      logical addimp, obsimp, sigoimp, hassens, sclimp, imp0hr, dfs, enorm, esigo, so2xmsn2so
      logical regular

      character*255 RCfname   ! name of rc file
      character*255 odsfilex  ! master ods file
      character*255 odsfiley  ! matching ods file
      character*255 odsfilez  ! output ods file
      character*255 odsunmatch! output for unmatched data
      character*80  ftype 


!     storage for ODS:
!     ---------------
      type ( ods_vect ) odsx, odsy
      integer, allocatable :: is(:), js(:)
      integer :: lsthkx(MAXHOBS)
      integer :: lsthkt(MAXHOBS)

      real   , pointer :: latx(:), laty(:)
      real   , pointer :: lonx(:), lony(:)
      real   , pointer :: levx(:), levy(:)
      real   , pointer :: xlev(:), ylev(:)
      integer, pointer :: timx(:), timy(:)
      integer, pointer :: ktx (:), kty (:)
      integer, pointer :: kxx (:), kxy (:)
      integer, pointer :: qcxx(:), qcxy(:)
      integer, pointer :: qchx(:), qchy(:)
      real   , pointer :: xvcx(:), xvcy(:)
      real   , pointer :: xmx (:), xmy (:)
      real   , pointer :: oma (:)
      real   , pointer :: omf (:)
      real   , pointer :: omfx(:), omfy(:)
      real             :: ximp
      real(8) :: domf,doma

      integer, parameter :: im = 144
      integer, parameter :: jm = 91
      integer, parameter :: km = 72

      real     energy_scale
      external energy_scale

!     Parse command line and load resources
!     -------------------------------------
      call init ( addimp, obsimp, sigoimp, hassens, sclimp, imp0hr, dfs, enorm, esigo, so2xmsn2so,
     .            rcfname, odsfilex, odsfiley, odsfilez, odsunmatch )
      
!     Read RC file
!     ------------
      call getRC_ ( lsthkx, lsthkt )
      
!     Loop over all synoptic times on the master file
!     -----------------------------------------------
      do isyn = 1, 32767
      
!        Read all data for this synoptic time
!        ------------------------------------
         nymd = -1            ! get data for next synoptic time on file
         nhms =  0
         call ODSNxTime ( trim(odsfilex), nymd, nhms )
         if ( nymd .eq. -1 ) then 
              print *, 'End-Of-File'
              exit 
         end if	 
         call ODS_Get ( trim(odsfilex), nymd, nhms, ftype, odsx, ierr )
         if ( ierr .gt. 0 ) then
              print *, 'ODS_Get error: ierr = ', ierr, ' file = ', trim(odsfilex)
	      exit
         end if
	 
	 print *, '  nymd, nhms = ',  nymd, nhms
	 print *
	 
         call ODS_Get ( trim(odsfiley), nymd, nhms, ftype, odsy, ierr )
         if ( ierr .gt. 0 ) then
              print *, 'ODS_Get error: ierr = ', ierr, ' file = ', trim(odsfiley)
	      exit
         end if

         nx = odsx%data%nobs
         print *, 'Master   ODS: Number of obs = ', nx
         if ( nx .eq. 0 ) cycle

         latx => odsx%data%lat
         lonx => odsx%data%lon
         levx => odsx%data%lev
         timx => odsx%data%time
         ktx  => odsx%data%kt
         kxx  => odsx%data%kx
         xmx  => odsx%data%xm
         xvcx => odsx%data%xvec
	 
         qchx => odsx%data%qchist
         qcxx => odsx%data%qcexcl
         oma  => odsx%data%oma

         allocate ( xlev(nx), stat = ierr )
            if (ierr/=0) then
                print *, 'Allocate error: ierr = ', ierr, ' file = ', trim(odsfilex)
                call exit(7)
            endif
         xlev(1:nx) = levx(1:nx)
         do i = 1, nh
            where(lsthkt(i)==ktx.and.lsthkx(i)==kxx) xlev=xmx
         end do

	 
!        Set all oma in master ODS to missing	 
!        ------------------------------------	 
         regular=.not.(addimp.or.obsimp.or.sigoimp.or.imp0hr.or.dfs.or.so2xmsn2so)
	 if(regular          ) oma(1:nx) = obs_missing
	 if(obsimp.or.sigoimp) xvcx(1:nx) = 0.0
	 if(sclimp           ) xvcx(1:nx) = obs_missing
	 
         ny = odsy%data%nobs
         print *, 'Matching ODS: Number of obs = ', ny
         if ( ny .eq. 0 ) cycle

         laty => odsy%data%lat
         lony => odsy%data%lon
         levy => odsy%data%lev
         timy => odsy%data%time
         kty  => odsy%data%kt
         kxy  => odsy%data%kx
         qchy => odsy%data%qchist
         qcxy => odsy%data%qcexcl
         xmy  => odsy%data%xm
         xvcy => odsy%data%xvec
	 
         if (.not.(addimp.or.obsimp)) omf  => odsy%data%omf

         if (sclimp.or.hassens.or.imp0hr.or.dfs) then
             omfx => odsx%data%omf
             omfy => odsy%data%omf
         endif

         allocate ( ylev(ny), stat = ierr )
            if (ierr/=0) then
                print *, 'Allocate error: ierr = ', ierr, ' file = ', trim(odsfiley)
                call exit(7)
            endif
         ylev(1:ny) = levy(1:ny)
         do i = 1, nh
            where(lsthkt(i)==kty.and.lsthkx(i)==kxy) ylev=xmy
         end do
	 
!        Sort according to matching attributes
!        -------------------------------------
         allocate ( is(nx) )
         call IndexSet  ( nx, is )
         call IndexSort ( nx, is, latx(1:nx), descend=.false. )
         call IndexSort ( nx, is, lonx(1:nx), descend=.false. )
         call IndexSort ( nx, is, xlev(1:nx), descend=.false. )
         call IndexSort ( nx, is, timx(1:nx), descend=.false. )
         call IndexSort ( nx, is,  ktx(1:nx), descend=.false. )
         call IndexSort ( nx, is,  kxx(1:nx), descend=.false. )
 
         allocate ( js(ny) )
         call IndexSet  ( ny, js )
         call IndexSort ( ny, js, laty(1:ny), descend=.false. )
         call IndexSort ( ny, js, lony(1:ny), descend=.false. )
         call IndexSort ( ny, js, ylev(1:ny), descend=.false. )
         call IndexSort ( ny, js, timy(1:ny), descend=.false. )
         call IndexSort ( ny, js,  kty(1:ny), descend=.false. )
         call IndexSort ( ny, js,  kxy(1:ny), descend=.false. )
 
!        Find matching items
         i = 1
         j = 1
         k = 0
         do while (i .le. nx) ! for each obs in odsx:
            ii = is(i)
            do jt = j, ny     ! look for a match in odsy
               jj = js(jt)
               do
                  d = kxx (ii)-kxy (jj); if (d /= 0) exit
                  d = ktx (ii)-kty (jj); if (d /= 0) exit
                  d = timx(ii)-timy(jj); if (d /= 0) exit
                  d = xlev(ii)-ylev(jj); if (d /= 0) exit
                  d = lonx(ii)-lony(jj); if (d /= 0) exit
                  d = latx(ii)-laty(jj);             exit
               end do 
               if (d==0) then   ! all attributes match
                  if (addimp) then
                    if ( hassens ) then  ! if files have sensitivity ...
                                         ! ... calculate impacts and add
                         if( qcxx(ii)==0 .and. qcxy(jj)==0 ) then
                             xvcx(ii) = omfx(ii)*xvcx(ii) + omfy(jj)*xvcy(jj)
                         else if ( qcxx(ii)==0 .and. qcxy(jj)/=0 ) then
                            xvcx(ii) = omfx(ii)*xvcx(ii) 
                         else if ( qcxx(ii)/=0 .and. qcxy(jj)==0 ) then
                            xvcx(ii) = omfy(jj)*xvcy(jj) 
                            qcxx(ii) = 0
                         else
                            xvcx(ii) = 0.0
                         endif
                    else                 ! if files have impacts ...
                                         ! ... just add them
                         if( qcxx(ii)==0 .and. qcxy(jj)==0 ) then
                             xvcx(ii) = xvcx(ii) + xvcy(jj)
                         else if ( qcxx(ii)==0 .and. qcxy(jj)/=0 ) then
                            xvcx(ii) = xvcx(ii) 
                         else if ( qcxx(ii)/=0 .and. qcxy(jj)==0 ) then
                            xvcx(ii) = xvcy(jj) 
                            qcxx(ii) = 0
                         else
                            xvcx(ii) = 0.0
                         endif
                    endif
                  else if ( sclimp ) then
                    if( qcxx(ii)==0 .and. qcxy(jj)==0 ) then
                        xvcx(ii) = omfx(ii)*omfx(ii) - omfy(jj)*omfy(jj)
                        if ( esigo ) then
                           xvcx(ii) = xvcx(ii) / (xvcy(jj)*xvcy(jj)) ! scale by R^{-1}
                        endif
                        if ( enorm ) then
                           xvcx(ii) = energy_scale(esigo,im,jm,km,latx(ii),ktx(ii)) * xvcx(ii)
                        endif
                    else
                        xvcx(ii) = 0.0
                    endif
                  else if ( imp0hr ) then  ! second file irrelevant
                    if( qcxx(ii)==0 .and. xvcx(ii)<obs_missing) then
                       doma = oma (ii)
                       domf = omfx(ii)
                       if(abs(doma-obs_missing)<1.d-5.or.abs(domf-obs_missing)<1.d-5) then ! passive data is sometimes present
                          xvcx(ii) = 0.0
                       else
                          ximp     = doma*doma - domf*domf
                          xvcx(ii) = ximp / (xvcx(ii)*xvcx(ii)) ! scale by R^{-1}
                       endif
                    else
                       xvcx(ii)  = 0.0
                    endif
                  else if ( dfs ) then  ! second file irrelevant
                    if( qcxx(ii)==0 .and. xvcx(ii)<obs_missing) then
                       doma = oma (ii)
                       domf = omfx(ii)
                       if(abs(doma-obs_missing)<1.d-5.or.abs(domf-obs_missing)<1.d-5) then ! passive data is sometimes present
                          xvcx(ii) = 0.0
                       else
                          ximp     = (domf - doma)*doma ! = [h(xa)-h(xb)]*oma
                          xvcx(ii) = ximp / (xvcx(ii)*xvcx(ii)) ! scale by R^{-1}
                       endif
                    else
                       xvcx(ii)  = 0.0
                    endif
                  else if ( obsimp ) then
                    if ( qcxy(jj)==0 ) then
                         if (hassens) then
                             xvcx(ii) = omfy(jj)*xvcy(jj) ! calculate observation impact
                          else
                             xvcx(ii) = xvcy(jj)          ! 2nd file already contains desired obs impact
                         endif
                    else
                         xvcx(ii) = 0.0
                    endif
                  else if ( sigoimp ) then
                    if( qcxx(ii)==0 .and. qcxy(jj)==0 ) then
	                xvcx(ii) = -oma(ii)*xvcy(jj)
                    else
	                xvcx(ii) = 0.0
                    endif
                  else
                    if ( so2xmsn2so ) then
                       xmx(ii)  = xvcx(jj)
         	       xvcx(ii) = xvcy(jj)
                    else
                       if ( qcxx(ii)==0 .and. qcxy(jj)==0 ) then
             	          oma(ii)  = omf(jj)
                       else
             	          if(qcxx(ii)==0) qcxx(ii) = X_ODSMATCH ! if qcxx already non-zero, leave it alone
                       endif
                    endif
                  endif
	          j = min(jt + 1,ny)
	          i = i + 1
	          k = k + 1
	          exit         ! stop looking in odsy
 	       elseif (d<0 .OR. jt==ny) then  ! no match exists for this obs
	          j = jt
	          i = i + 1
	          exit         ! stop looking in odsy	       
	       end if
	    end do             ! loop on second file
	 end do                ! loop on first  file
	 
	 print *, '              Number of matched obs = ', k
	 
         if ( sclimp ) then
            where(odsx%data%xvec==obs_missing) odsx%data%qcexcl=X_ODSMATCH
         endif
      
         call ODS_Put ( odsfilez, ftype, nymd, nhms, odsx, ierr )

         call unmatched_()

         deallocate(is,js)
         deallocate(ylev)
         deallocate(xlev)

      end do  ! loop over synoptic times

      if ( associated(latx) ) nullify(latx,lonx,levx,timx,ktx,kxx)
      if ( associated(laty) ) nullify(laty,lony,levy,timy,kty,kxy)
      if ( associated(oma)  ) nullify(oma)
      if ( associated(omf)  ) nullify(omf)
      if ( associated(omfx) ) nullify(omfx)
      if ( associated(omfy) ) nullify(omfy)
      
      stop

      contains

      subroutine getRC_ ( lsthkx, lsthkt )

      integer :: lsthkx(:)
      integer :: lsthkt(:)

      character,parameter :: myname_ = myname//"::getRC_"
      integer i, i1, lt, k1, k2, iret, jret
      integer kxnext, ktnext, irow, icnt
      character*255 token,tablename
      logical allkxkts 

!     load resource file
!     ------------------
      call i90_loadf (trim(rcfname), iret)
          if (iret/=0) then
              print*, 'Error reading rc file: ', trim(rcfname)
              call exit(7)
          endif

!     Read table with instruments types
!     ---------------------------------
      nh = 0
      allkxkts = .false.
      tablename = 'list_height_based::'
      call I90_label(trim(tablename), iret)
      if (iret/=0) then
         write(6,'(4a)') myname_, ': table ', trim(tablename),
     .                        'not found in RC file ... will take all KXs'
         allkxkts = .true.
      end if
      if ( .not. allkxkts ) then
        irow = 0; icnt = 0
        do while (iret==0)                     ! read table entries
           call I90_GLine ( iret )             ! iret=-1: end of file; +1: end of table
           if (iret==0.and.irow<MAXHOBS) then    ! OK, we have next row of table
               irow = irow + 1
  
               call I90_GToken ( token, jret ) ! obs data type (kt)
               if (jret/=0) then
                   write(6,'(2a,i5)') myname_, ': I90_GToken error, jret=', jret
               end if
               read(token,*) ktnext
  
               jret=0
               do  j = 1, kxmax
                 call I90_GToken(token, jret )
                 if(jret/=0) exit
                 ii = index(token,':') ! token is single entry or range of entries
                 lt = len_trim(token)
                 if (ii==0) then       ! no colon, therefore single entry
                     read(token,*) kxnext
                     icnt = icnt + 1
                     if (icnt==MAXHOBS) then    ! check space
                         write(6,'(2a,i5)') myname,': increase MAXHOBS'
                         stop(999)
                     else if (jret==0) then
                         lsthkx(icnt) = kxnext
                         lsthkt(icnt) = ktnext ! same kt for all kx's in this row
                     end if
                 else                  ! colon, therefore k1:k2
                     read(token(1:ii-1),*) k1
                     read(token(ii+1:lt),*) k2
                     do kxnext = k1, k2
                        icnt = icnt + 1
                        if (icnt==MAXHOBS) then    ! check space
                            write(6,'(2a,i5)') myname,': increase MAXHOBS'
                            stop(999)
                        else if (jret==0) then
                            lsthkx(icnt) = kxnext
                            lsthkt(icnt) = ktnext ! same kt for all kx's in this row
                        end if
                     end do
                 end if
               enddo

           end if
        end do
        nh = icnt
        print *, 'Will process the following data as height-based obs:'
        print *, '          kt                 kx  '
        do j = 1, nh
           print *, lsthkt(j), lsthkx(j)
        enddo
      endif ! < all KXs >


!     release resource file
!     ---------------------
      call i90_release()

      end subroutine getRC_

      subroutine unmatched_

      integer ni,nj,it,ntot,iii,jjj
      type ( ods_vect ) odsu

      if ( trim(odsunmatch) == 'NONE' ) return

!     Take 1st file as master (odsx) ...
!     ----------------------------------
         i = 1
         j = 1
         k = 0
         nj= 0
         do while (i .le. nx) ! for each obs in odsx:
            ii = is(i)
            do jt = j, ny     ! look for a match in odsy
               jj = js(jt)
               do
                  d = kxx (ii)-kxy (jj); if (d /= 0) exit
                  d = ktx (ii)-kty (jj); if (d /= 0) exit
                  d = timx(ii)-timy(jj); if (d /= 0) exit
                  d = xlev(ii)-ylev(jj); if (d /= 0) exit
                  d = lonx(ii)-lony(jj); if (d /= 0) exit
                  d = latx(ii)-laty(jj);             exit
               end do 
               if (d==0) then   ! all attributes match
                  j = min(jt + 1,ny)
                  i = i + 1
                  k = k + 1
                  exit         ! stop looking in odsy
               elseif (d<0 .OR. jt==ny) then  ! no match exists for this obs
                  j = jt
                  i = i + 1
                  nj= nj+ 1
                  qchx(ii) = 28
                  exit         ! stop looking in odsy	       
               endif
            end do
         end do

!     Take 2nd file as master (odsx) ...
!     ----------------------------------
         i = 1
         j = 1
         k = 0
         ni= 0
         do while (j .le. ny) ! for each obs in odsy:
            jj = js(j)
            do it = i, nx     ! look for a match in odsy
               ii = is(it)
               do
                  d = kxy (jj)-kxx (ii); if (d /= 0) exit
                  d = kty (jj)-ktx (ii); if (d /= 0) exit
                  d = timy(jj)-timx(ii); if (d /= 0) exit
                  d = ylev(jj)-xlev(ii); if (d /= 0) exit
                  d = lony(jj)-lonx(ii); if (d /= 0) exit
                  d = laty(jj)-latx(ii);             exit
               end do 
               if (d==0) then   ! all attributes match
                  i = min(it + 1,nx)
                  j = j + 1
                  k = k + 1
                  exit         ! stop looking in odsx
               elseif (d<0 .OR. it==nx) then  ! no match exists for this obs
                  i = it
                  j = j + 1
                  ni= ni+ 1
                  qchy(jj) = 29
                  exit         ! stop looking in odsx	       
               endif
            end do
         end do
      print *, 'number of obs in 2nd file not in 1st: ',ni
      print *, 'number of obs in 1st file not in 2nd: ',nj

!     Initialized output ods structure
!     --------------------------------
      ntot=ni+nj
      call ODS_Init ( odsu, ntot, ierr, 
     .                odsx%meta%nkt, odsx%meta%nkx, odsx%meta%nqc, odsx%meta%ncr, nsyn=odsx%meta%nsyn )

!     copy metadata from input ODS
!     ----------------------------
      odsu%meta%kt_names  = odsx%meta%kt_names
      odsu%meta%kt_units  = odsx%meta%kt_units
      odsu%meta%kx_names  = odsx%meta%kx_names
      odsu%meta%kx_meta   = odsx%meta%kx_meta
      odsu%meta%qcx_names = odsx%meta%qcx_names

      ii=0
      do i=1,nx
         iii = is(i)
         if(qchx(iii)==28) then
           ii=ii+1
           odsu%data%kid(ii)=ii
           odsu%data%kx(ii)=odsx%data%kx(iii)
           odsu%data%kt(ii)=odsx%data%kt(iii)
           odsu%data%ks(ii)=odsx%data%ks(iii)
           odsu%data%xm(ii)=odsx%data%xm(iii)
           odsu%data%lat(ii)=odsx%data%lat(iii)
           odsu%data%lon(ii)=odsx%data%lon(iii)
           odsu%data%lev(ii)=odsx%data%lev(iii)
           odsu%data%obs(ii)=odsx%data%obs(iii)
           odsu%data%omf(ii)=odsx%data%omf(iii)
           odsu%data%oma(ii)=odsx%data%oma(iii)
           odsu%data%xvec(ii)=odsx%data%xvec(iii)
           odsu%data%time(ii)=odsx%data%time(iii)
           odsu%data%qcexcl(ii)=odsx%data%qcexcl(iii)
           odsu%data%qchist(ii)=odsx%data%qchist(iii)
         endif
      enddo
      jj=ii
      do j=1,ny
         jjj = js(j)
         if(qchy(jjj)==29) then
           jj=jj+1
           odsu%data%kid(jj)=jj
           odsu%data%kx(jj)=odsy%data%kx(jjj)
           odsu%data%kt(jj)=odsy%data%kt(jjj)
           odsu%data%ks(jj)=odsy%data%ks(jjj)
           odsu%data%xm(jj)=odsy%data%xm(jjj)
           odsu%data%lat(jj)=odsy%data%lat(jjj)
           odsu%data%lon(jj)=odsy%data%lon(jjj)
           odsu%data%lev(jj)=odsy%data%lev(jjj)
           odsu%data%obs(jj)=odsy%data%obs(jjj)
           odsu%data%omf(jj)=odsy%data%omf(jjj)
           odsu%data%oma(jj)=odsy%data%oma(jjj)
           odsu%data%xvec(jj)=odsy%data%xvec(jjj)
           odsu%data%time(jj)=odsy%data%time(jjj)
           odsu%data%qcexcl(jj)=odsy%data%qcexcl(jjj)
           odsu%data%qchist(jj)=odsy%data%qchist(jjj)
         endif
      enddo
      if(jj/=odsu%data%nobs) then
         print*, 'Something is amiss, aborting ...'
         stop
      endif
      print*, 'odd data pieces ',jj-count(odsu%data%qchist==28)-count(odsu%data%qchist==29)
      if ( odsu%data%nobs>0 ) then
         print*, 'total number of unmatched observations: ',odsu%data%nobs
         print*, 'qch = 28: data found in 2nd file but absent in 1st file '
         print*, 'qch = 29: data found in 1st file but absent in 2nd file '
         call ODS_Put ( trim(odsunmatch), ftype, nymd, nhms, odsu, ierr )
      else
         print*, 'all data matched, no unmatched output to write'
      endif

      end subroutine unmatched_

      end ! program odsmatch

!EOC

      real function energy_scale ( esigo,im,jm,km,lat,kt )
      implicit none

      logical, intent(in) :: esigo
      integer, intent(in) :: im
      integer, intent(in) :: jm
      integer, intent(in) :: km
      real,    intent(in) :: lat
      integer, intent(in) :: kt

!     These have to be wired in to prevent meaningless
!     dependency of ODS on hermes
!     ------------------------------------------
      real, parameter ::  cpm    = 1004.64            ! J / (kg * K)
      real, parameter ::  alhl   = 2.499e6            ! water latent hear of vaporization (J/kg)
      real, parameter ::  rgas   = 287.04
      real, parameter ::  pstd   = 1000.0
      real, parameter ::  tstd   = 280.0              ! K
      real, parameter ::  cpO3   = 817.               ! O3 heat capacity at cnst pressure J/(kg*K)
      real, parameter ::  o3mw   = 3 * 15.9994 / 1000.! molecular weight of O3 (kg/mole)
      real, parameter ::  o3Lv   = 316.3e3            ! latent heat of vaporization (J/kg)
      real, parameter ::  O3Mass = 3.e15              ! Apprx Mass O3 in Atmosphere (g)
      real, parameter ::  o3ppmv2= 100.               ! peak stratospheric mix ratio squared

      logical, parameter :: grid = .true.

      real, save ::  scale_nil, scale_pressure, scale_temperature
      real, save ::  scale_wvapor, scale_o3mixratio, scale_wind
      real, save ::  scale_speed, scale_pcp
      real, save ::  pi, dx, dy, dz
      real hgrid,ggrid,cosmid,radlat,radlatm,radlatp,cosp,cosm,o3moles,o3Lv2,alhl2

      logical, save :: first = .true.

      if ( first ) then ! Coeffs are missing some sort of cosine weighting
	 print *, '              Scaling by total energy coefficients'
         pi = 4.0*atan(1.0)
!        dx=360./im; dy=90./(jm-1.);dz=1./km
         dx=2.*pi/im; dy=pi/(jm-1.); dz=1./km
         o3moles = O3Mass / o3mw
         o3Lv2   = o3Lv * o3Lv / o3ppmv2                 ! note: OmF for O3 is ppmv
         alhl2   = alhl * alhl * 1.e-6                   ! note: OmF for Q in g/kg transformed to kg/kg
! E norm defined with 1/2 ...
         scale_nil        = 0.0
         if(esigo) then ! when R^{-1} is used to scale omf, the following gives them units of J/kg
            scale_pressure   = 0.5 * rgas * tstd * pstd**2; !print*, "pressure = ", scale_pressure
            scale_temperature= 0.5 * cpm * tstd;            !print*, "temperature = ", scale_temperature
            scale_wvapor     = 0.5 * alhl2 * (cpm * tstd);  !print*, "wvapor = ", scale_wvapor
            scale_speed      = 0.5;                         !print*, "speed = ", scale_speed
            scale_wind       = 0.5;                         !print*, "winds = ", scale_wind
            scale_o3mixratio = 0.5 * o3Lv2*(cpo3 * tstd)  ;  print*, "o3mxrat = ", scale_o3mixratio
            scale_pcp        = 0.5;                          print*, "pcp     = ", scale_pcp
         else           ! these are the actual energy-like scaling coefficients
            scale_pressure   = 0.5 * rgas * tstd / pstd**2; !print*, "pressure = ", scale_pressure
            scale_temperature= 0.5 * cpm / tstd;            !print*, "temperature = ", scale_temperature
            scale_wvapor     = 0.5 * alhl2 / (cpm * tstd);  !print*, "wvapor = ", scale_wvapor
            scale_speed      = 0.5;                         !print*, "speed = ", scale_speed
            scale_wind       = 0.5;                         !print*, "winds = ", scale_wind
            scale_o3mixratio = 0.5 * o3Lv2/(cpo3 * tstd)  ;  print*, "o3mxrat = ", scale_o3mixratio
            scale_pcp        = 0.5;                          print*, "pcp     = ", scale_pcp
         endif

         first = .false.
      endif

! need to do poles right ...
      if ( grid ) then
           radlat  = lat*pi/180.
           radlatm = radlat-dy
           radlatp = radlat+dy
!          cosmid  = 1.0*(cos(radlatp)-cos(radlatm)) / dy
           cosp    = (sin(radlatp)-sin(radlat )) / dy
           cosm    = (sin(radlat )-sin(radlatm)) / dy
           cosmid  = 0.5 * (cosm + cosp)
           hgrid = cosp*dx*dy
           ggrid = cosp*hgrid*dz
      else
           hgrid = 1.0
           ggrid = 1.0
      endif

      select case (kt)
      case (4)                 ! u
         energy_scale = ggrid*scale_wind
      case (5)                 ! v
         energy_scale = ggrid*scale_wind
      case (11)                ! water vapor
         energy_scale = ggrid*scale_wvapor
      case (12)                ! speed
         energy_scale = ggrid*scale_speed
      case (17)                ! pcp
         energy_scale = hgrid*scale_pcp
      case (21)                ! ozone
         energy_scale = ggrid*scale_o3mixratio
      case (22)                ! ozone
         energy_scale = ggrid*scale_o3mixratio
      case (33)                ! surf pressure
         energy_scale = hgrid*scale_pressure
      case (40)                ! temperature and tb
         energy_scale = ggrid*scale_temperature
      case (44)                ! temperature and tb
         energy_scale = ggrid*scale_temperature
      case default
         energy_scale = scale_nil
      end select

      end function energy_scale

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!-----------------------------------------------------------------------
!BOP
! !ROUTINE: init: initialize odsmatch
!
! !DESCRIPTION:
!
! !INTERFACE:
!
      subroutine init ( addimp, obsimp, sigoimp, hassens, sclimp, imp0hr, dfs, enorm, esigo, 
     .                  so2xmsn2so, rcfname,
     .                  odsfilex, odsfiley, odsfilez, odsunmatch )

      implicit NONE
!
! !OUTPUT PARAMETERS:

      logical,       intent(out) :: addimp
      logical,       intent(out) :: obsimp
      logical,       intent(out) :: sigoimp
      logical,       intent(out) :: hassens
      logical,       intent(out) :: sclimp
      logical,       intent(out) :: imp0hr
      logical,       intent(out) :: dfs
      logical,       intent(out) :: enorm
      logical,       intent(out) :: esigo
      logical,       intent(out) :: so2xmsn2so
      character*255, intent(out) :: rcfname
      character*255, intent(out) :: odsfilex
      character*255, intent(out) :: odsfiley
      character*255, intent(out) :: odsfilez
      character*255, intent(out) :: odsunmatch
!
!
! !REVISION HISTORY:
!
!     02May2005 D.Dee   - Initial code
!     23Feb2010 Todling - Add imp0hr (merge with das215)
!     21Apr2013 Todling - Add DFS (degrees of freedom for signal)
!     14Aig2013 Todling - Add so2xmsn2so: move sigo to xm slot; 
!                         place sensitivity from 2nd file in sigo slot;
!                         this is so odsstats may work without change
!
!EOP
!BOC

      character(len=*), parameter :: myname = 'init'
      character(len=*), parameter :: RCfile = 'odsmatch.rc'

      integer, parameter :: nfiles_max=3
      integer i, iarg, argc, iargc, nfiles
      character*255 argv
      character*255 infile(nfiles_max)

      addimp  = .false.
      obsimp  = .false.
      sigoimp = .false.
      hassens = .false.
      sclimp  = .false.
      imp0hr  = .false.
      dfs     = .false.
      enorm   = .false.
      esigo   = .false.
      so2xmsn2so = .false.
      rcfname = RCfile
      odsunmatch = 'NONE'
      
!     Parse command line
!     ------------------
      argc =  iargc()
      if ( argc .lt. nfiles_max ) call usage()
      nfiles = 0
      iarg = 0
      do i = 1, 32767
         iarg = iarg + 1
         if ( iarg .gt. argc ) go to 111
         call GetArg ( iArg, argv )
         if (index(argv,'-addimp' ) .gt. 0 ) then
             addimp = .true.
         else if (index(argv,'-obsimp' ) .gt. 0 ) then
             obsimp = .true.
         else if (index(argv,'-sigoimp' ) .gt. 0 ) then
             sigoimp = .true.
         else if (index(argv,'-hassens' ) .gt. 0 ) then
             hassens = .true.
         else if (index(argv,'-sclimp' ) .gt. 0 ) then
             sclimp = .true.
         else if (index(argv,'-imp0hr' ) .gt. 0 ) then
             imp0hr = .true.
         else if (index(argv,'-dfs' ) .gt. 0 ) then
             dfs = .true.
         else if (index(argv,'-enorm' ) .gt. 0 ) then
             enorm = .true.
         else if (index(argv,'-esigo' ) .gt. 0 ) then
             esigo = .true.
         else if (index(argv,'-so2xmsn2so' ) .gt. 0 ) then
             so2xmsn2so = .true.
         else if (index(argv,'-unmatched' ) .gt. 0 ) then
             iarg = iarg + 1
             call GetArg ( iArg, odsunmatch )
         else if (index(argv,'-rc' ) .gt. 0 ) then
             iarg = iarg + 1
             call GetArg ( iArg, rcfname )
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
      if ( nfiles .lt. nfiles_max ) call usage()

      odsfilex = infile(1)
      odsfiley = infile(2)
      odsfilez = infile(3)
 
!     Consistency check
!     -----------------
      if ( addimp .and. obsimp ) call usage()
      if ( addimp .and. sclimp ) call usage()
      if ( addimp .and. sigoimp) call usage()
      if ( obsimp .and. sclimp ) call usage()
      if ( obsimp .and. sigoimp) call usage()
      if ( obsimp .and. imp0hr ) call usage()
      if ( sclimp .and. imp0hr ) call usage()
      if ( dfs    .and. imp0hr ) call usage()
      if ( obsimp .and. dfs    ) call usage()
      if ( sclimp .and. dfs    ) call usage()

!     Echo the parameters
!     -------------------
      print *
      print *, 'Master   file: ', trim(odsfilex)
      print *, 'Matching file: ', trim(odsfiley)
      print *, 'Output   file: ', trim(odsfilez)
      print *

      return

      end ! subroutine init

!EOC

      subroutine usage()
      print *
      print *, 'odsmatch - Find matching observations in two ODS files.'
      print *
      print *, 'Usage:  odsmatch [options] odsfilex odsfiley odsfilez'
      print *
      print *, 'where'
      print *
      print *, ' -hassens      indicates 2nd input ods file has sensitivities in Xvec slot'
      print *, ' -addimp       adds oma slots of master and matching file'
      print *, ' -obsimp       places obs impacts in the xvec slot of output file'
      print *, ' -sigoimp      places sigmaO impacts in the xvec slot of output file'
      print *, ' -sclimp       calculate impact from omf and scale per energy weights'
      print *, '                 imp = omf(t|ta)Eomf(t|ta)-omf(t|tb)Eomf(t|tb) '
      print *, '                 beware of sequence of files in command line'
      print *, ' -imp0hr      calculate 0-hr impact: omaR^{-1}oma-omfR^{-1}omf'
      print *, ' -dfs         calculate degrees of freedom for signal (DFS)'
      print *, ' -esigo       calculate impact from omf and scale with R^{-1}'
      print *, '                 imp = omf(t|ta)R^{-1}omf(t|ta)-omf(t|tb)R^{-1}omf(t|tb) '
      print *, '                 beware of sequence of files in command line'
      print *, '                    >>> see note 6'
      print *, ' -so2xmsn2so  move sigo from master to xm slot and sensitivity from 2nd file '
      print *, '              to sigo slot of master; leave oma slot untouched; '
      print *, '              CAUTION:  this only applies in special cases'
      print *, ' -rc RCfile   name of resource file (default: odsmatch.rc)'
      print *, ' -unmatched   filename to write out unmatched data to (default: NONE)'
      print *
      print *, ' odsfilex      master ODS file, oma to be redefined'
      print *, ' odsfiley      matching ODS file, to be searched for matching entries'
      print *, ' odsfilez      name of output ODS file'
      print *
      print *, 'NOTES'
      print *, '1) Both input files must be true ODS post-anal files (not GSI diag files)'
      print *, '2) A match is defined if all of the following attributes are equal:'
      print *, '        obs, kx, kt, time, lev, lon, lat.'
      print *, '   If an observation in the matching file is found to match an observation'
      print *, '   in the master file, then its omf attribute is copied to the oma attribute'
      print *, '   in the output file.'
      print *, '3) When -obsimp opt is used, the xvec slot of the output file contains the'
      print *, '   obs impacts from matching file; all other entries remain as in master file'
      print *, '4) When -addimp opt is used, the xvec slot of the output file contains the'
      print *, '   the sum of the xvec slots from master and matching file for matched obs'
      print *, '5) Ops -addimp and -obsimp cannot be used together'
      print *, '6) When opts esigo and sclimp are used together the impact is given by:'
      print *, '      imp = omf(t|ta)DR^{-1}Domf(t|ta)-omf(t|tb)DR^{-1}Domf(t|tb) '
      print *, '   whenre D=sqrt(E) with E being diagonal.'
      print *, '7) Opts addimp and sens do not go together, for obvious reasons'
      print *, '8) DFS cannot be mixed up with diagnostic type (such as imp0hr, obsimp, etc)'
      stop
      end
