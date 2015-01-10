!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!-----------------------------------------------------------------------
!BOP
!
! !ROUTINE: ods_rmdupl:  Remove duplicates from an ODS structure.
!
! !INTERFACE:

      subroutine ods_rmdupl ( ods )

      use m_ods,          only : ods_vect
      use m_die,          only : die
      use m_SortingTools, only : indexSet,indexSort

      implicit NONE
      
! !INPUT/OUTPUT PARAMETERS:

      type(ods_vect),   intent(inout) :: ods
      
! !DESCRIPTION:
!
!     Removes duplicates from an ods structure. Duplicates are
!     defined based on equality of kx,kt,time,lat,lon,lev,obs 
!
! !REVISION HISTORY:
!
!     09Dec2004 (Dee) - original code
!
!EOP
!BOC

      character*10, parameter :: myname = 'ods_rmdupl'

!     Pointers to ods attributes:
!     --------------------------
      integer, pointer :: kt(:)       ! data type index
      integer, pointer :: kx(:)       ! data source index
      integer, pointer :: ks(:)       ! sounding index
      real   , pointer :: lon(:)      ! longitude of obs (degrees)
      real   , pointer :: lat(:)      ! latitude of obs (degrees)
      real   , pointer :: lev(:)      ! pressure level of obs (hPa)
      integer, pointer :: time(:)     ! time
      real   , pointer :: obs(:)      ! observation
      real   , pointer :: OmF(:)      ! obs-minus-fcst (O-F)
      real   , pointer :: OmA(:)      ! obs-minus-ana  (O-A)
      real   , pointer :: Xvec(:)     ! metadata
      real   , pointer :: xm(:)       ! metadata
      integer, pointer :: qcx(:)      ! exclusion mark
      integer, pointer :: qch(:)      ! history mark
      
      integer, allocatable :: indx(:)
      integer :: nobs, ier, i, is

!     Assign pointers to ODS attributes
!     ---------------------------------
      kt   => ods%data%kt
      kx   => ods%data%kx
      ks   => ods%data%ks
      lon  => ods%data%lon
      lat  => ods%data%lat
      lev  => ods%data%lev
      time => ods%data%time
      obs  => ods%data%obs
      OmF  => ods%data%OmF
      OmA  => ods%data%OmA
      Xvec => ods%data%Xvec
      xm   => ods%data%xm
      qcx  => ods%data%qcexcl
      qch  => ods%data%qchist
      
      nobs = ods%data%nobs
      
      allocate ( indx(nobs), stat=ier )
      if ( ier/=0 ) call die ( myname, 'allocate()', ier )

      call IndexSet  ( nobs, indx )
      call IndexSort ( nobs, indx,  obs(1:nobs), descend=.false. )
      call IndexSort ( nobs, indx,  lat(1:nobs), descend=.false. )
      call IndexSort ( nobs, indx,  lon(1:nobs), descend=.false. )
      call IndexSort ( nobs, indx,  lev(1:nobs), descend=.false. )
      call IndexSort ( nobs, indx, time(1:nobs), descend=.false. )
      call IndexSort ( nobs, indx,   kt(1:nobs), descend=.false. )
      call IndexSort ( nobs, indx,   kx(1:nobs), descend=.false. )
  
      kt  (1:nobs) = kt  ( (/ (indx(i), i=1,nobs) /) )
      kx  (1:nobs) = kx  ( (/ (indx(i), i=1,nobs) /) )
      ks  (1:nobs) = ks  ( (/ (indx(i), i=1,nobs) /) )
      lon (1:nobs) = lon ( (/ (indx(i), i=1,nobs) /) )
      lat (1:nobs) = lat ( (/ (indx(i), i=1,nobs) /) )
      lev (1:nobs) = lev ( (/ (indx(i), i=1,nobs) /) )
      time(1:nobs) = time( (/ (indx(i), i=1,nobs) /) )
      obs (1:nobs) = obs ( (/ (indx(i), i=1,nobs) /) )
      OmF (1:nobs) = OmF ( (/ (indx(i), i=1,nobs) /) )
      OmA (1:nobs) = OmA ( (/ (indx(i), i=1,nobs) /) )
      xm  (1:nobs) = xm  ( (/ (indx(i), i=1,nobs) /) )
      qcx (1:nobs) = qcx ( (/ (indx(i), i=1,nobs) /) )
      qch (1:nobs) = qch ( (/ (indx(i), i=1,nobs) /) )
      Xvec(1:nobs) = Xvec( (/ (indx(i), i=1,nobs) /) )

      is = 1
      indx(is) = 1
      do i = 2, nobs
         if ( kx  (i)==kx  (i-1) .AND. kt (i)==kt (i-1) .AND. 
     .        time(i)==time(i-1) .AND. lat(i)==lat(i-1) .AND.
     .        lon (i)==lon (i-1) .AND. lev(i)==lev(i-1) .AND.
     .        obs (i)==obs (i-1) ) then
!             print *, '***Duplicate: kx, kt, time, lat, lon, lev, obs = ', 
!     .                    kx(i), kt(i), time(i), lat(i), lon(i), lev(i), obs(i)
             cycle
	 end if
         is = is + 1
         indx(is) = i
      end do

      nobs = is

      do is = 1, nobs
         i  = indx(is)
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
         omf(is)  = omf(i)
         oma(is)  = oma(i)
         xvec(is) = xvec(i)
      end do
      
      ods%data%nobs = nobs
      
      deallocate ( indx )
           
!     Nullify pointers
!     -----------------
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
            
      end subroutine ods_rmdupl
