!-------------------------------------------------------------------------
!         NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!-------------------------------------------------------------------------
!BOP
!
! !MODULE:  m_odsdata --- Implements observation data stream data class.
!
! !INTERFACE:
!

   module  m_MPodsdata

! !USES:

   use m_MPodsmeta
   use m_mergesorts

   implicit none

! !PUBLIC TYPES:
!
   PRIVATE
   PUBLIC  obs_vect        ! atomic attributes (lat, lon, etc.)
!
! !PUBLIC MEMBER FUNCTIONS:
!
   PUBLIC  Init,    odsd_init
   PUBLIC  Clean,   odsd_clean
   PUBLIC  Resize,  odsd_resize
   PUBLIC  MaskOut, odsd_MaskOut, ods_MaskOut   ! keep old names for backward 
   PUBLIC  Reorder, odsd_Reorder, obs_Reorder   ! compatibility 
   PUBLIC  MoveUp,  odsd_MoveUp,  obs_MoveUp
   PUBLIC  Permute, odsd_Permute, obs_Permute
   PUBLIC  Sort,    odsd_Sort,    obs_Sort
!
! !PUBLIC DATA MEMBERS:
!
   real, parameter :: obs_missing = 1.0E+15 ! obs.  missing value

   PUBLIC  obs_missing      

!
! !DESCRIPTION:
!
!  This module defines the observation data stream data class.
!  It relies on the ODS library and HDF.
!
! !REVISION HISTORY:
!
!  10Jun2002 Zaslavsky  Created from m_ods module.
!
!  18Jun2002 Zaslavsky nsyn moved from m_odsmeta module to this module.
!
!EOP
!-------------------------------------------------------------------------

!  Internal dimensions
!  -------------------
 
   integer, parameter :: NOBS_MAX = 450 * 1000
   integer, parameter :: NSYN_DEF =   4

!  Observation atomic atttributes (lat, lon, lev, etc.)
!  ----------------------------------------------------
   type obs_vect

      integer          :: nobs      ! actual number of observations
      integer          :: nhalo     ! # obs in halo region (owned by other PEs)
      integer          :: nvct      ! vector size allocated (nobs .le. nvct)
      integer          :: nymd      ! date (year-month-date)
      integer          :: nhms      ! time (hour-mimute-second)
      integer          :: nsyn      ! actual number of synoptic times per day

      integer, pointer :: gid(:)    ! Obs identification index (global)
      integer, pointer :: kid(:)    ! Obs identification index

      real,    pointer :: lat(:)    ! latitute     of obs (degrees)
      real,    pointer :: lon(:)    ! longitude    of obs (degrees)
      real,    pointer :: lev(:)    ! level        of obs (hPa)
      integer, pointer :: kx(:)     ! data source index
      integer, pointer :: kt(:)     ! data type   index
      integer, pointer :: ks(:)     ! sounding    index
      real,    pointer :: xm(:)     ! atomic metadata (depends on kx)
      integer, pointer :: time(:)   ! time (relative to the input/output
                                    !   synoptic date and time)
      real,    pointer :: obs(:)    ! observation value (units depend on kt)
      real,    pointer :: OmF(:)    ! obs minus forecast (O-F) innovations
      real,    pointer :: OmA(:)    ! obs minus analysis (O-A) residuals
      real,    pointer :: Xvec(:)   ! PSAS CG solution vector
      integer, pointer :: qcexcl(:) ! On-line QC exclusion flag
      integer, pointer :: qchist(:) ! On-line QC history flag

      end type obs_vect

   Interface Init                ! Renamed for backward compatibility
      module procedure ODSD_Init
   end Interface

   Interface Clean                ! Renamed for backward compatibility
      module procedure ODSD_Clean
   end Interface

   Interface Resize
      module procedure ODSD_Resize
   end Interface

   Interface ODS_MaskOut                ! Renamed for backward compatibility
      module procedure ODSD_Maskout
   end Interface

   Interface Obs_Reorder                ! Renamed for backward compatibility
      module procedure ODSD_Reorder
   end Interface

   Interface Obs_MoveUp                 ! Renamed for backward compatibility
      module procedure ODSD_MoveUp
   end Interface

   Interface Obs_Permute
      module procedure ODSD_Permute
   end Interface

   Interface MaskOut           
      module procedure ODSD_Maskout
   end Interface

   Interface Reorder                
      module procedure ODSD_Reorder
   end Interface

   Interface MoveUp                
      module procedure ODSD_MoveUp
   end Interface

   Interface Permute
      module procedure ODSD_Permute
   end Interface

   Interface Sort                
      module procedure ODSD_Sort
   end Interface

   Interface Obs_Sort                   ! Renamed for backward compatibility
      module procedure ODSD_Sort
   end Interface


CONTAINS

!-------------------------------------------------------------------------
!         NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE:  ODSD_Init --- Allocate memory for ODS data.
!
! !INTERFACE:
!
  subroutine ODSD_Init ( y, nobs, rc,    &
                         nsyn )          ! optional

! !USES:

  implicit none

! !INPUT PARAMETERS:
  integer, intent(in)           :: nobs    ! number of observations; if
                                           ! nobs<0  nobs will be reset to
                                           ! a large internal value.
                                           ! does not include halo obs!!

  integer, intent(in), optional :: nsyn    ! number of synoptic times per
                                           ! day                            

! !INPUT/OUTPUT PARAMETERS:

  type(obs_vect), intent(inout) ::  y      ! ODS vector to be allocated

! !OUTPUT PARAMETERS:

  integer, intent(out)          ::  rc     ! Error return code:
                                           !  0 - all is well
                                           !  2 - error allocating attributes

! !DESCRIPTION:
!
!  Allocates memory for an ODS vector.
!
! !REVISION HISTORY:
! 
!  10 Jun2002  Zaslavsky  Created from ODS_Init
!  19Jun2002   Zaslavsky  Added checks to avoid memory leaks.
!
!EOP
!-------------------------------------------------------------------------

  character(len=*), parameter :: myname = 'odsd_init'
  integer ios, nvct, i, nnsyn, nobs_tot


! Avoid memory leaks
! ------------------
  if ( associated ( y%gid    ) .or.  &
       associated ( y%kid    ) .or.  &
       associated ( y%lat    ) .or.  &
       associated ( y%lon    ) .or.  &
       associated ( y%lev    ) .or.  &
       associated ( y%kx     ) .or.  &
       associated ( y%kt     ) .or.  &
       associated ( y%ks     ) .or.  &
       associated ( y%xm     ) .or.  &
       associated ( y%time   ) .or.  &
       associated ( y%obs    ) .or.  &
       associated ( y%OmF    ) .or.  &
       associated ( y%OmA    ) .or.  &
       associated ( y%Xvec   ) .or.  &
       associated ( y%qcexcl ) .or.  &
       associated ( y%qchist ) )     &
  then
       rc = 1 
       return 
  end if

! Check for halo region
! ---------------------

  nobs_tot = nobs
  y%nhalo = 0

! Reset nobs if user specifies less than zero value
! -------------------------------------------------
  if ( nobs_tot .lt. 0 ) nobs_tot = NOBS_MAX


! Allocated memory for ODS atomic attributes
! ------------------------------------------
  if ( nobs_tot .gt. 0 ) then
       nvct = nobs_tot
  else
       nvct = 1
  endif

  if ( present(nsyn) ) then
       nnsyn = nsyn
  else
       nnsyn = NSYN_DEF
  end if

  y%nobs = nobs
  y%nvct = nvct
  y%nymd = -1
  y%nhms = -1
  y%nsyn = nnsyn

  allocate    ( y%gid(nvct),      &
                y%kid(nvct),      &
                y%lat(nvct),      &
                y%lon(nvct),      &
                y%lev(nvct),      &
                y%kx(nvct),       &
                y%kt(nvct),       &
                y%ks(nvct),       &
                y%xm(nvct),       &
                y%time(nvct),     &
                y%obs(nvct),      &
                y%OmF(nvct),      &
                y%OmA(nvct),      &
                y%Xvec(nvct),     &
                y%qcexcl(nvct),   &
                y%qchist(nvct),   &
                stat = ios )
     if ( ios .ne. 0 ) then
        rc = 2
        return
     end if

! Set observation index
! ---------------------
  y%gid = (/ (i, i=1,nvct) /)  ! global ID may be overwritten later
  y%kid = (/ (i, i=1,nvct) /)

! All done
! --------
  rc = 0
  return

  end subroutine odsd_init



!-------------------------------------------------------------------------
!         NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE:  ODSD_Clean - Deallocates ODS data
!
! !INTERFACE:
!
      subroutine ODSD_Clean ( y, rc )

! !USES:

      implicit none

! !INPUT/OUTPUT PARAMETERS:

  type(obs_vect), intent(inout) ::  y      ! ODS data to be allocated

! !OUTPUT PARAMETERS:

  integer, intent(out)          ::  rc     ! Error return code:
                                           !  0 - all is well
                                           !  2 - error deallocating attributes

! !DESCRIPTION:
!
!  Frees memory used by an ODS data.
!
! !REVISION HISTORY:
!
!  10Jun2002   Zaslavsky  Created from ODS_Clean.
!  19Jun2002   Zaslavsky  Added checks to avoid memory leaks.
!
!EOP
!-------------------------------------------------------------------------

  character(len=*), parameter :: myname = 'odsd_clean'
  integer ios1, ios2, ios3, ios4, ios5, ios6, ios7, ios8, ios9, ios10, &
          ios11, ios12, ios13, ios14, ios15, ios0 

  rc = 0

! Deallocated memory for ODS atomic attributes
! --------------------------------------------
  y%nobs = 0
  y%nvct = 0
  y%nymd = -1
  y%nhms = -1
  y%nsyn = 0
  
  if ( associated ( y%gid ) ) deallocate ( y%gid,  stat = ios0 )
  if ( associated ( y%kid ) ) deallocate ( y%kid,  stat = ios1 )     
  if ( associated ( y%lat ) ) deallocate ( y%lat, stat = ios2 )       
  if ( associated ( y%lon ) ) deallocate ( y%lon,  stat = ios3 )     
  if ( associated ( y%lev ) ) deallocate ( y%lev, stat = ios4 )      
  if ( associated ( y%kx ) ) deallocate ( y%kx, stat = ios5 )       
  if ( associated ( y%kt ) ) deallocate ( y%kt,  stat = ios6 )      
  if ( associated ( y%ks ) ) deallocate ( y%ks, stat = ios7 )       
  if ( associated ( y%xm ) ) deallocate ( y%xm, stat = ios8 )       
  if ( associated ( y%time ) ) deallocate ( y%time, stat = ios9 )     
  if ( associated ( y%obs ) ) deallocate ( y%obs, stat = ios10 )      
  if ( associated ( y%OmF ) ) deallocate ( y%OmF, stat = ios11 )      
  if ( associated ( y%OmA ) ) deallocate ( y%OmA, stat = ios12 )      
  if ( associated ( y%Xvec ) ) deallocate ( y%Xvec, stat = ios13 )     
  if ( associated ( y%qcexcl ) ) deallocate ( y%qcexcl, stat = ios14 )   
  if ( associated ( y%qchist ) ) deallocate ( y%qchist, stat = ios15 )   

  if ( ( ios0  .ne. 0 ) .or. &
       ( ios1  .ne. 0 ) .or. &
       ( ios2  .ne. 0 ) .or. &
       ( ios3  .ne. 0 ) .or. &
       ( ios4  .ne. 0 ) .or. &
       ( ios5  .ne. 0 ) .or. &
       ( ios6  .ne. 0 ) .or. &
       ( ios7  .ne. 0 ) .or. &
       ( ios8  .ne. 0 ) .or. &
       ( ios9  .ne. 0 ) .or. &
       ( ios10 .ne. 0 ) .or. &
       ( ios11 .ne. 0 ) .or. &
       ( ios12 .ne. 0 ) .or. &
       ( ios13 .ne. 0 ) .or. &
       ( ios14 .ne. 0 ) .or. &
       ( ios15 .ne. 0 ) )    &
  then
         rc = rc + 1
  end if


  return

  end subroutine  odsd_clean



!-------------------------------------------------------------------------
!         NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE:  ODSD_Resize --- Revise the size of the ODS data set
!
! !INTERFACE:
!
  subroutine ODSD_Resize( data, nobs_new, is_halo, rc )

! !USES:
  use parutilitiesmodule, only : gid
  implicit none

! !INPUT PARAMETERS:

  integer, intent(in)           ::  nobs_new ! new size
  logical, intent(in)           ::  is_halo  ! if halo region

! !INPUT/OUTPUT PARAMETERS:
  type(obs_vect), intent(inout) ::  data   ! Previous ODS data

! !OUTPUT PARAMETERS:
  integer, intent(out)          ::  rc     ! Error return code:
                                           !  0 - all is well
                                           !  1 - inconsistent halo size

! !DESCRIPTION:
! \label{MODS:Enlarge}
!  Revises the size of an ODS vector data segment.  If the new
!  size is larger and it is a halo region, nobs stays the same
!  and nhalo makes up for the additional size, but if it is not
!  a halo, nobs still stays the same (the user is responsible for
!  changing this) and nhalo is 0.  If the new size is smaller and
!  and a halo is specified, then nobs should still be equal to
!  the new size (if not the size is inconsistent rc=1) and the
!  halo is set back to zero.  If the new size is smaller and
!  halo is not specified, the last obs are thrown away.
!
! !BUGS:
!  This routine is extremely dangerous to use directly by the
!  user, size it can make the ODS decomposition invalid.  The
!  use of ODS_Resize, ODS_Ghost and ODS_Unghost is strongly
!  advised.
!
! !REVISION HISTORY:
!
!  25Oct2002   Sawyer     Creation from ODS_Ghost
!
!EOP
!-------------------------------------------------------------------------

  character(len=*), parameter :: myname = 'odsd_resize'
  integer nobs, i
  integer, pointer :: tmp_int(:)
  real, pointer    :: tmp_real(:)

  rc = 0
  nobs = data%nobs
  if ( nobs_new >= nobs) then
     if ( is_halo ) data%nhalo = nobs_new - nobs
!        if not halo, changing data%nobs is up to the user!!!
  else  !  nobs_new < nobs
     if ( is_halo ) then
        data%nhalo = 0
        rc = 1   ! this case is inconsistent, but not wrong
     endif
     nobs = nobs_new         ! observations are being thrown away
     data%nobs = nobs
  end if

  data%nvct = nobs_new

    allocate( tmp_int(nobs_new) )
    tmp_int(1:nobs) = data%gid(1:nobs)
    deallocate( data%gid )
    data%gid => tmp_int
    nullify( tmp_int )

    allocate( tmp_int(nobs_new) )
    tmp_int(1:nobs) = data%kid(1:nobs)
    deallocate( data%kid )
    data%kid => tmp_int
    nullify( tmp_int )

    allocate( tmp_real(nobs_new) )
    tmp_real(1:nobs) = data%lat(1:nobs)
    deallocate( data%lat )
    data%lat => tmp_real
    nullify( tmp_real )

    allocate( tmp_real(nobs_new) )
    tmp_real(1:nobs) = data%lon(1:nobs)
    deallocate( data%lon )
    data%lon => tmp_real
    nullify( tmp_real )

    allocate( tmp_real(nobs_new) )
    tmp_real(1:nobs) = data%lev(1:nobs)
    deallocate( data%lev )
    data%lev => tmp_real
    nullify( tmp_real )

    allocate( tmp_int(nobs_new) )
    tmp_int(1:nobs) = data%kx(1:nobs)
    deallocate( data%kx )
    data%kx => tmp_int
    nullify( tmp_int )

    allocate( tmp_int(nobs_new) )
    tmp_int(1:nobs) = data%kt(1:nobs)
    deallocate( data%kt )
    data%kt => tmp_int
    nullify( tmp_int )

    allocate( tmp_int(nobs_new) )
    tmp_int(1:nobs) = data%ks(1:nobs)
    deallocate( data%ks )
    data%ks => tmp_int
    nullify( tmp_int )

    allocate( tmp_real(nobs_new) )
    tmp_real(1:nobs) = data%xm(1:nobs)
    deallocate( data%xm )
    data%xm => tmp_real
    nullify( tmp_real )

    allocate( tmp_int(nobs_new) )
    tmp_int(1:nobs) = data%time(1:nobs)
    deallocate( data%time )
    data%time => tmp_int
    nullify( tmp_int )

    allocate( tmp_real(nobs_new) )
    tmp_real(1:nobs) = data%obs(1:nobs)
    deallocate( data%obs )
    data%obs => tmp_real
    nullify( tmp_real )

    allocate( tmp_real(nobs_new) )
    tmp_real(1:nobs) = data%OmF(1:nobs)
    deallocate( data%OmF )
    data%OmF => tmp_real
    nullify( tmp_real )

    allocate( tmp_real(nobs_new) )
    tmp_real(1:nobs) = data%OmA(1:nobs)
    deallocate( data%OmA )
    data%OmA => tmp_real
    nullify( tmp_real )

    allocate( tmp_real(nobs_new) )
    tmp_real(1:nobs) = data%Xvec(1:nobs)
    deallocate( data%Xvec )
    data%Xvec => tmp_real
    nullify( tmp_real )

    allocate( tmp_int(nobs_new) )
    tmp_int(1:nobs) = data%qcexcl(1:nobs)
    deallocate( data%qcexcl )
    data%qcexcl => tmp_int
    nullify( tmp_int )

    allocate( tmp_int(nobs_new) )
    tmp_int(1:nobs) = data%qchist(1:nobs)
    deallocate( data%qchist )
    data%qchist => tmp_int
    nullify( tmp_int )

  return

  end subroutine odsd_resize


!-------------------------------------------------------------------------
!         NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE:  ODSD_MaskOut - Mask out observations given
!
! !INTERFACE:
!
      subroutine ODSD_MaskOut ( y, nobs, boxes, nboxes, X_FLAG, qcx, nset, &
                               outside )   ! optional

! !USES:

      implicit none

! !INPUT PARAMETERS:

      type(obs_vect), intent(inout)  :: y       ! obs vector
      integer,     intent(in)        :: nobs    ! no. of obs
      integer,     intent(in)        :: nboxes  ! no. of boxes

                                                ! Hyper cubes:
      real,        intent(in)        :: boxes(2,6,*)
                                                ! First dimension:
                                                !   1 - beginning
                                                !   2 - end
                                                ! Second dimension:
                                                !   1 - kt range
                                                !   2 - kx range
                                                !   3 - latitude  range
                                                !   4 - longitude range
                                                !   5 - level range
                                                !   6 - time stamp range



      integer,     intent(in)        :: X_FLAG  ! QC flag value to set

      logical, intent(in), OPTIONAL  :: outside ! If true, observations
                                                !  outside boxes will be
                                                !  masked out

! !OUTPUT PARAMETERS:

  integer,        intent(out)  ::  qcx(nobs)     ! Typically, qc exclusion
                                                 !  flag

  integer,        intent(out)  ::  nset          ! no. of obs set with X_FLAG

! !DESCRIPTION:
!
!  Maskout observations given attribute ranges. By default, the output
!  vector is set with {\tt X\_FLAG} if the observation lies INSIDE the
!  boxes. When this routine is called with {\tt outside = .true.}, the
!  output vector is set with {\tt X\_FLAG} if the observation lies
!  OUTSIDE the boxes.
!
!
! !REVISION HISTORY:
!
!  03dec1999   da Silva   Initial code
!  13dec1999   da Silva   Removed timestamp check for now.
!  30jan2002 Sienkiewicz  m_ods time problem fixed - uncomment time check
!  10jun2002 Zaslavsky    Renamed from ODS_MaskOut to ODSD_MaskOut. Became 
!                         a part of m_odsdata module.
!
!EOP
!-------------------------------------------------------------------------

  character(len=*), parameter :: myname = 'odsd_maskout'

  logical inunion, inside
  integer n, k

  real rlondn,rlonup
  real rlatdn,rlatup
  real rlevdn,rlevup
  integer kxdn,kxup
  integer ktdn,ktup
  integer tmdn,tmup

  kxdn(k)   = nint(boxes(1,1,k))	! Data source beg
  kxup(k)   = nint(boxes(2,1,k))	! Data source end
  ktdn(k)   = nint(boxes(1,2,k))	! Data type beg
  ktup(k)   = nint(boxes(2,2,k))	! Data type end
  rlatdn(k) = boxes(1,3,k)		! South
  rlatup(k) = boxes(2,3,k)		! North
  rlondn(k) = boxes(1,4,k)		! West
  rlonup(k) = boxes(2,4,k)		! East
  rlevdn(k) = boxes(2,5,k)		! Low (high pressure)
  rlevup(k) = boxes(1,5,k)		! High (low pressure)
  tmdn(k)   = nint(boxes(1,6,k))	! earliest time stamp
  tmup(k)   = nint(boxes(2,6,k))	! latest time stamp



! Maskout observations inside boxes, unless caller wants otherwise
! ----------------------------------------------------------------
  inside = .true.
  if ( present(outside) ) then
       inside = (.not. outside)
  end if

  nset = 0
  if ( nboxes .eq. 0 .or. nobs .eq. 0 ) return

  do n = 1, nobs

     if ( y%qcexcl(n) .eq. 0 ) then ! only examine obs which are CLEAR

        k=0
        inunion=.false.
        do while(.not.inunion.and.k.lt.nboxes)
           k=k+1
           inunion=                                                   &
     		rlondn(k).le. y%lon(n).and. y%lon(n).le.rlonup(k) .and.  &
     		rlatdn(k).le. y%lat(n).and. y%lat(n).le.rlatup(k) .and.  &
                kxdn(k).le.  y%kx(n).and.  y%kx(n).le.  kxup(k) .and.  &
                ktdn(k).le.  y%kt(n).and.  y%kt(n).le.  ktup(k)  .and.  &
                tmdn(k).le.y%time(n).and.y%time(n).le.  tmup(k)

           !             Do not check levels for sfc obs
           !             -------------------------------
           if(inunion .and. y%kt(n).ne.ktSLP .and. &
                y%kt(n).ne.ktUS .and. &
                y%kt(n).ne.ktVS )     &
                inunion =                  &
    		rlevup(k).le. y%lev(n).and. y%lev(n).le.rlevdn(k)

        end do

!       Mark qcx whether obs is inside or outside of boxes
!       --------------------------------------------------
        if ( inside ) then
           if ( inunion ) then
                qcx(n) = X_FLAG
                nset = nset + 1
           end if
        else
           if ( .not. inunion ) then
                qcx(n) = X_FLAG
                nset = nset + 1
           end if
        end if


     endif

   end do

end subroutine ODSD_MaskOut




!-------------------------------------------------------------------------
!         NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE:  ODSD_MoveUp() --- Move selected observations to front.
!
! !INTERFACE:
!
      subroutine ODSD_MoveUp ( y, nobs, nmoved, rc,  &
                                 attrname, attrvalues, &
                                 idx )       ! optional

!
! !INPUT PARAMETERS:
!
      integer , intent(in)      ::  nobs       ! total number of obs
      character*(*), intent(in) ::  attrname   ! name of attribute
      integer , intent(in)      ::  attrvalues(:)  ! selection values
!
! !INPUT/OUTPUT PARAMETERS:
!
      type ( obs_vect ) , intent (inout) :: y
!
! !OUTPUT PARAMETERS:
!
      integer , intent(out)     ::  nmoved     ! number of data moved up
      integer , intent(out)     ::  rc         ! rc=0: ok
                                               ! rc=1: allocate problem
                                               ! rc=2: invalid attribute
                                               ! rc=3: not enough space for idx
      integer , intent(out), optional ::  idx(:)  ! reordering index
!
! !DESCRIPTION: Moves observations with selected attribute values to the
!     front of the list. Only works for integer attributes.
!
!
! !REVISION HISTORY:
!
!  10May2002   Dee     Initial code based on GEOS-2 ODS_MoveUp
!  14May2002   Dee     Added optional output parameter idx
!  10jun2002 Zaslavsky Renamed from OBS_MoveUp to ODSD_MoveUp
!
!EOP
!-------------------------------------------------------------------------

!     Local work space
!     ----------------
      integer, allocatable ::  indx(:)
      logical, allocatable ::  selected(:)

      character(len=*), parameter :: myname='OBS_MoveUp'
      integer ios, i, nvalues

      rc = 0
      nmoved = 0

!     Nothing to do
!     -------------
      if ( nobs == 0 ) return

!     Allocate memory
!     ---------------
      allocate ( indx(nobs), selected(nobs), stat=ios )
      if ( ios/=0 ) then
          rc = 1
          return
      end if

      selected = .false.
      nvalues = size(attrvalues)

      select case (attrname)
      case ('kid')
           do i = 1, nvalues
              where (y%kid    ==attrvalues(i)) selected = .true.
           end do
      case ('kt')
           do i = 1, nvalues
              where (y%kt     ==attrvalues(i)) selected = .true.
           end do
      case ('kx')
           do i = 1, nvalues
              where (y%kx     ==attrvalues(i)) selected = .true.
           end do
      case ('ks')
           do i = 1, nvalues
              where (y%ks     ==attrvalues(i)) selected = .true.
           end do
      case ('time')
           do i = 1, nvalues
              where (y%time   ==attrvalues(i)) selected = .true.
           end do
      case ('qcexcl')
           do i = 1, nvalues
              where (y%qcexcl ==attrvalues(i)) selected = .true.
           end do
      case ('qchist')
           do i = 1, nvalues
              where (y%qchist ==attrvalues(i)) selected = .true.
           end do
      case default
           rc = 2   ! invalid attribute name
           return
      end select

      nmoved = count(selected)

      indx = (/ (i,i=1,nobs) /)
      if ( nmoved > 0 ) then
           indx = (/ pack(indx,mask=selected), pack(indx,mask=.not. selected) /)
           call ODSD_Reorder ( nobs, y, indx )
      end if

!     Return indx if requested
!     ------------------------
      if ( present ( idx ) ) then
           if ( size(idx) >= nobs ) then
                idx(1:nobs) = indx
           else
                rc = 3   ! not enough space
           end if
      end if

!     Deallocate memory
!     ---------------
      deallocate ( selected, indx )

      end subroutine ODSD_MoveUp

!-------------------------------------------------------------------------
!         NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE:  ODSD_Reorder() --- Reorder observation vector.
!
! !INTERFACE:
!
      subroutine ODSD_Reorder ( ntot, y, indx, is_ghosted )

!
! !INPUT PARAMETERS:
!
      integer , intent(in)               :: ntot          ! number of obs
      logical, intent(in), optional      :: is_ghosted    ! ignore ghost region
!
! !INPUT/OUTPUT PARAMETERS:
!
      type ( obs_vect ) , intent (inout) :: y             ! data
      integer , intent(inout)            :: indx(ntot)    ! sorting index
!
! !DESCRIPTION:  Reorders the first ntot elements of each array
!                attribute of the observation vector y according
!                to the index vector indx.  If the region is ghosted
!                {\tt is\_ghosted} the result is split up into a
!                ghosted and non-ghosted part and the permutation
!                {\tt indx} is changed accordingly.  Only use this
!                option if operating in SPMD mode, and if you really
!                understand the background behind ghosted observations.
!
! !REVISION HISTORY:
!
!  10May2002 (Dee)   Initial code based on GEOS-2 ODS_Reorder
!  10Jun2002 Zaslavsky  Renamed from OBS_Reorder to ODSD_Reorder
!  04Nov2002 Sawyer  Added is_ghosted to ignore ghost region in reorder
!
!EOP
!-------------------------------------------------------------------------

      character(len=*), parameter :: myname='OBS_Reorder'

      integer i, count, counthalo
      integer, allocatable :: rev_tmp(:)

!     Nothing to do
!     -------------
      if ( ntot==0 ) return

!     Apply sorting permutation to data arrays
!     ----------------------------------------
      if ( present(is_ghosted) ) then
!
! The index array is split into ghost and inherent sections
!
         allocate( rev_tmp(ntot) )
         do i=1,ntot
            rev_tmp(indx(i)) = i   ! Rev_tmp is the reverse permutation
         enddo
         count = 0             ! point to the start of the observations
         counthalo = y%nobs    ! point to the start of halo region
         do i=1, ntot
            if ( rev_tmp(i) > y%nobs ) then   ! part of halo region
               counthalo = counthalo+1
               indx(rev_tmp(i)) = counthalo
            else                              ! part of inherent region
               count = count+1
               indx(rev_tmp(i)) = count
            endif
         enddo
         print *, "counthalo", counthalo, " should be ", y%nobs+y%nhalo
         print *, "count", count, " should be ", y%nobs
         deallocate( rev_tmp )
      endif
      y%kid    (1:ntot) = y%kid    ( (/ (indx(i), i=1,ntot) /) )
      y%kt     (1:ntot) = y%kt     ( (/ (indx(i), i=1,ntot) /) )
      y%kx     (1:ntot) = y%kx     ( (/ (indx(i), i=1,ntot) /) )
      y%ks     (1:ntot) = y%ks     ( (/ (indx(i), i=1,ntot) /) )
      y%lon    (1:ntot) = y%lon    ( (/ (indx(i), i=1,ntot) /) )
      y%lat    (1:ntot) = y%lat    ( (/ (indx(i), i=1,ntot) /) )
      y%lev    (1:ntot) = y%lev    ( (/ (indx(i), i=1,ntot) /) )
      y%time   (1:ntot) = y%time   ( (/ (indx(i), i=1,ntot) /) )
      y%obs    (1:ntot) = y%obs    ( (/ (indx(i), i=1,ntot) /) )
      y%OmF    (1:ntot) = y%OmF    ( (/ (indx(i), i=1,ntot) /) )
      y%OmA    (1:ntot) = y%OmA    ( (/ (indx(i), i=1,ntot) /) )
      y%Xvec   (1:ntot) = y%Xvec   ( (/ (indx(i), i=1,ntot) /) )
      y%xm     (1:ntot) = y%xm     ( (/ (indx(i), i=1,ntot) /) )
      y%qcexcl (1:ntot) = y%qcexcl ( (/ (indx(i), i=1,ntot) /) )
      y%qchist (1:ntot) = y%qchist ( (/ (indx(i), i=1,ntot) /) )

      return

      End subroutine ODSD_Reorder


!-------------------------------------------------------------------------
!         NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE:  ODSD_Permute() --- Permute index vector using attribute
!
! !INTERFACE:
!
      subroutine ODSD_Permute( y, nobs, nmoved, rc,                        &
                              attrname, attrvalues, perm ) 

!
! !INPUT PARAMETERS:
!
      integer , intent(in)      ::  nobs       ! active obs
      type ( obs_vect ) , intent (in) :: y     ! observations
      character*(*), intent(in) ::  attrname   ! name of attribute
      integer , intent(in)      ::  attrvalues(:)  ! selection values
!
! !INPUT/OUTPUT PARAMETERS:
!
      integer , intent(inout)   ::  perm(:)    ! permutation vector
!
! !OUTPUT PARAMETERS:
!
      integer , intent(out)     ::  nmoved     ! # permutations <= nobs
      integer , intent(out)     ::  rc         ! rc=0: ok
                                               ! rc=1: allocate problem
                                               ! rc=2: invalid attribute
                                               ! rc=3: not enough space for idx
!
! !DESCRIPTION:   This routine is closely related to OBS_MoveUp, however
!    \begin{enumerate} 
!       \item the observation vector is never reshuffled,
!       \item the permutation vector provides indirect addressing
!             into the obs vector.
!       \item the permutation vector may be longer than nobs,
!             and may already be permuted,
!       \item nmoved <=  nobs <= size(perm)
!    \end{enumerate}
!
! !REVISION HISTORY:
!
!  10Sep2002 Sawyer   Initial code based on ODS_MoveUp
!
!EOP
!-------------------------------------------------------------------------

!     Local work space
!     ----------------
      integer, allocatable ::  indx(:)
      logical, allocatable ::  selected(:)
      logical :: flag

      character(len=*), parameter :: myname='ODSD_Permute'
      integer ios, i, j, nvalues, ntotal

      rc = 0
      nmoved = 0
      ntotal = size(perm)

!     Nothing to do
!     -------------
      if ( nobs == 0 ) return

!     Allocate memory
!     ---------------
      allocate ( indx(nobs), selected(nobs), stat=ios )
      if ( ios/=0 ) then
          rc = 1
          return
      end if

      selected = .false.
      nvalues = size(attrvalues)

      select case (attrname)
      case ('kid')
        do j = 1, nobs
           do i = 1, nvalues
              if (y%kid(perm(j)) == attrvalues(i)) selected(j) = .true.
           end do
        end do
      case ('kt')
        do i = 1, nvalues
           do j = 1, nobs
              if (y%kt(perm(j)) == attrvalues(i)) selected(j) = .true.
           end do
        end do
      case ('kx')
        do j = 1, nobs
           do i = 1, nvalues
              if (y%kx(perm(j)) == attrvalues(i)) selected(j) = .true.
           end do
        end do
      case ('ks')
        do j = 1, nobs
           do i = 1, nvalues
              if (y%ks(perm(j)) == attrvalues(i)) selected(j) = .true.
           end do
        end do
      case ('time')
        do j = 1, nobs
           do i = 1, nvalues
              if (y%time(perm(j)) == attrvalues(i)) selected(j) = .true.
           end do
        end do
      case ('qcexcl')
        do j = 1, nobs
           do i = 1, nvalues
              if (y%qcexcl(perm(j)) == attrvalues(i)) selected(j) = .true.
           end do
        end do
      case ('qchist')
        do j = 1, nobs
           do i = 1, nvalues
              if (y%qchist(perm(j)) == attrvalues(i)) selected(j) = .true.
           end do
        end do
      case default
           rc = 2   ! invalid attribute name
           return
      end select

      nmoved = count(selected)

      indx = (/ (i,i=1,nobs) /)
      if ( nmoved > 0 ) then
        indx = (/ pack(indx,mask=selected), pack(indx,mask=.not. selected) /)
      end if

      perm(1:nobs) = perm( (/ (indx(i), i=1,nobs) /) )

!     Deallocate memory
!     ---------------
      deallocate ( selected, indx )

      end subroutine ODSD_Permute

!-------------------------------------------------------------------------
!         NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE:  ODSD_Sort - Sort observations by attribute values
!
! !INTERFACE:
!
   subroutine ODSD_Sort ( y     , nobs  , rc    , attr1 ,         &
                         attr2 , attr3 , attr4 , attr5 , attr6 , & ! optional
                         attr7 , attr8 , attr9 , attr10, attr11, & ! optional
                         attr12, attr13, attr14, attr15,         & ! optional
                         idx, descend )                            ! optional

! !USES:

  implicit none

! !INPUT PARAMETERS:

  type(obs_vect), intent(inout) :: y       ! obs vector
  integer, intent(in)           :: nobs    ! number of obs
  character(len=*), intent(in)  :: attr1   ! first attribute to sort on

  character(len=*), intent(in), optional :: attr2  ! next attribute to sort on
  character(len=*), intent(in), optional :: attr3  ! next attribute to sort on
  character(len=*), intent(in), optional :: attr4  ! next attribute to sort on
  character(len=*), intent(in), optional :: attr5  ! next attribute to sort on
  character(len=*), intent(in), optional :: attr6  ! next attribute to sort on
  character(len=*), intent(in), optional :: attr7  ! next attribute to sort on
  character(len=*), intent(in), optional :: attr8  ! next attribute to sort on
  character(len=*), intent(in), optional :: attr9  ! next attribute to sort on
  character(len=*), intent(in), optional :: attr10 ! next attribute to sort on
  character(len=*), intent(in), optional :: attr11 ! next attribute to sort on
  character(len=*), intent(in), optional :: attr12 ! next attribute to sort on
  character(len=*), intent(in), optional :: attr13 ! next attribute to sort on
  character(len=*), intent(in), optional :: attr14 ! next attribute to sort on
  character(len=*), intent(in), optional :: attr15 ! next attribute to sort on

  logical, intent(in), optional :: descend  ! sort in descending order if true

! !OUTPUT PARAMETERS:

  integer, intent(out), optional :: idx(:)     ! sort index vector
  integer, intent(out)           :: rc         ! return code:
                                               !  0: all ok
                                               !  1: allocate error
                                               !  2: invalid attribute
                                               !  3: not enough space for idx

! !DESCRIPTION: Sorts the input observation vector y according based on
!   one or more attribute values. Sorts on multiple attribute values are
!   performed in the order in which they appear in the list of calling
!   arguments. Optionally returns the sort index vector.
!
! !REVISION HISTORY:
!
!  20May2002   Dee        Initial code
!  30Jul2002   Zaslavsky  Moved to m_odsdata module and renamed.
!  30Jul2002   Zaslavsky  Included "gid" case.
!
!EOP
!-------------------------------------------------------------------------

   character(len=*), parameter :: myname = 'ODSD_Sort'
   integer,          parameter :: nattr  = 15

   integer, allocatable :: indx(:)
   character*6 attrname
   integer i,n

   allocate ( indx(nobs), stat=rc )
   if (rc/=0) then
       rc = 1
       return
   end if

   call IndexSet ( nobs, indx )

!  Sort on each attribute
!  ----------------------
   do n = 1, nattr

      if ( n== 1                       ) attrname = attr1
      if ( n== 2 .AND. present(attr2 ) ) attrname = attr2
      if ( n== 3 .AND. present(attr3 ) ) attrname = attr3
      if ( n== 4 .AND. present(attr4 ) ) attrname = attr4
      if ( n== 5 .AND. present(attr5 ) ) attrname = attr5
      if ( n== 6 .AND. present(attr6 ) ) attrname = attr6
      if ( n== 7 .AND. present(attr7 ) ) attrname = attr7
      if ( n== 8 .AND. present(attr8 ) ) attrname = attr8
      if ( n== 9 .AND. present(attr9 ) ) attrname = attr9
      if ( n==10 .AND. present(attr10) ) attrname = attr10
      if ( n==11 .AND. present(attr11) ) attrname = attr11
      if ( n==12 .AND. present(attr12) ) attrname = attr12
      if ( n==13 .AND. present(attr13) ) attrname = attr13
      if ( n==14 .AND. present(attr14) ) attrname = attr14
      if ( n==15 .AND. present(attr15) ) attrname = attr15

      select case (trim(attrname))
      case ('gid'    )
         call IndexSort ( nobs, indx, y%gid   (1:nobs), descend )
      case ('kid'    )
         call IndexSort ( nobs, indx, y%kid   (1:nobs), descend )
      case ('kt'     )
         call IndexSort ( nobs, indx, y%kt    (1:nobs), descend )
      case ('kx'     )
         call IndexSort ( nobs, indx, y%kx    (1:nobs), descend )
      case ('ks'     )
         call IndexSort ( nobs, indx, y%ks    (1:nobs), descend )
      case ('lon'    )
         call IndexSort ( nobs, indx, y%lon   (1:nobs), descend )
      case ('lat'    )
         call IndexSort ( nobs, indx, y%lat   (1:nobs), descend )
      case ('lev'    )
         call IndexSort ( nobs, indx, y%lev   (1:nobs), descend )
      case ('time'   )
         call IndexSort ( nobs, indx, y%time  (1:nobs), descend )
      case ('obs'    )
         call IndexSort ( nobs, indx, y%obs   (1:nobs), descend )
      case ('omf'    )
         call IndexSort ( nobs, indx, y%omf   (1:nobs), descend )
      case ('oma'    )
         call IndexSort ( nobs, indx, y%oma   (1:nobs), descend )
      case ('xvec'   )
         call IndexSort ( nobs, indx, y%xvec  (1:nobs), descend )
      case ('xm'     )
         call IndexSort ( nobs, indx, y%xm    (1:nobs), descend )
      case ('qcexcl' )
         call IndexSort ( nobs, indx, y%qcexcl(1:nobs), descend )
      case ('qchist' )
         call IndexSort ( nobs, indx, y%qchist(1:nobs), descend )
      case ('default')
         rc = 2
      end select

   end do

   call ODSD_Reorder ( nobs, y, indx )

!  Return index vector if requested
!  --------------------------------
   if ( present ( idx ) ) then
        if ( size(idx) >= nobs ) then
             idx(1:nobs) = indx
        else
             rc = 3   ! not enough space
        end if
   end if

   deallocate ( indx )

   end subroutine ODSD_Sort

end module m_MPodsdata
