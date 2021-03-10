!-------------------------------------------------------------------------
!         NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!-------------------------------------------------------------------------
!BOP
!
! !MODULE:  m_mpods --- Implements Parallel version of ODS using PILGRIM
! 
! !INTERFACE:
!

    module  m_mpods

! !USES:
   use decompmodule, only  : decomptype
   use ghostmodule, only   : ghosttype
   use m_zeit, only        : zeit_ci, zeit_co
   use m_MPodsmeta, only     : ods_meta, NCR_MAX
   use m_MPodsdata, only     : obs_vect, maskout, ods_maskout, obs_missing
   use m_odstransfer, only : transfertype
   implicit none

   PRIVATE
   PUBLIC  ODS_Meta           ! ODS global metadata (kx names, etc.)
   PUBLIC  Obs_Vect           ! atomic attributes (lat, lon, etc.)
   PUBLIC  ODS_Vect           ! ODS vector consisting of ODS_meta and obs_vect

! !PUBLIC MEMBER FUNCTIONS:
!
   PUBLIC  Init,    ODS_Init
   PUBLIC  Clean,   ODS_Clean
   PUBLIC  Ghost,   ODS_Ghost
   PUBLIC  Unghost, ODS_Unghost
   PUBLIC  Resize,  ODS_Resize
   PUBLIC  Refresh, ODS_Refresh
   PUBLIC  Put,     ODS_put
   PUBLIC  Get,     ODS_Get
   PUBLIC  Merge,   ODS_Merge
   PUBLIC  MaskOut, ODS_MaskOut   ! Deprecated!  use m_odsdata, only: ODSD_MaskOut
   PUBLIC  Tally,   ODS_Tally
   PUBLIC  Select,  ODS_Select
   PUBLIC  ODS_Transfer
   PUBLIC  ODS_Get1pe
   PUBLIC  ODS_Put1pe

!
! !PUBLIC DATA MEMBERS:
!
   PUBLIC  obs_missing      ! observation missing value
                            ! Deprecated!  use m_odsdata, only: obs_missing

!
! !DESCRIPTION:
! \label{MODS:Mod}
!  This module defines the observation data stream vector class.
!  It relies on the ODS library and HDF.
!  {\bf SPMD Notes:}  By definition {\it all} PUBLIC routines are
!  data parallel.  That is, they all work on or return distributed
!  ODS vectors (although it is possible to define decompositions
!  which happen to place all the data on one PE).  As a consequence
!  all PEs must call these routines in unison, even if they have
!  ``no'' work (e.g. nobs=0):  there may be collective communication
!  in which all PEs must participate.
!
!
! !REVISION HISTORY: 
!
!  02Sep1999 da Silva  Created from Ricardo's m_ods_structure, with several
!                      changes to remove circular dependencies on future 
!                      analysis modules, making it to reflect the data
!                      structures of an actual ODS file.
!  15Nov1999 da Silva  Added ktRH.
!  03Dec1999 da Silva  Added ODS_Merge()
!  05Jan2000 da Silva  Removed definition of kt???
!  13Dec2000 C. Redder Add the component, nsyn, to ods_meta, and the
!                      parameter constant, NSYN_DEF. 
!  22Jan2001 da Silva  Added Dick's ods_select()
!  16Feb2001 Todling   Added Xvec as an attribute of ODS.
!  21Mar2002 Sawyer    Added decomposition to ods_vect
!  10May2002 Sawyer    Added ODS_GetSPMD, ODS_PutSPMD wrappers
!  10May2002 Dee       Added Obs_Reorder and Obs_MoveUp (based on GEOS-2
!                      ODS_Reorder, ODS_MoveUp)
!  13May2002 Sawyer    Added transfertype, ODS_Transfer_init/Destroy.
!  07Jun2002 Sawyer    Removed SPMD distinctions, ODS_Merge SPMD complaint
!  10Jun2002 Zaslavsky Three modules m_odsmeta, m_odsdata and m_ods were
!                      created from the "old" m_ods.
!  18Jun2002 Zaslavsky Made changes related to relocation of nsyn from 
!                      metadata to data.
!  05Jul2002 Sawyer    Initial additions to support SPMD QC
!  26Aug2002 Zaslavsky Changed intent from OUT to INOUT for several objects. 
!  17Sep2002 Sawyer    Upgrade of documentation, passes protex, LaTeX
!  25Oct2002 Sawyer    Moved transfer facilities to m_odstransfer
!  04Dec2002 Sawyer    Added halo_trans to ods_vect; bug fixes
!  05Dec2002 Sawyer    Separated ODS_Ghost completely from ODS_Get
!  05Feb2003 Sawyer    Fixed bug in ODS_Merge, ODS_Get: lat_max bidirect.
!  13Feb2003 Sawyer    New: ODS_Partition_, Rewritten: ODS_GetM_ 
!  29May2003 Sawyer    ods_get1pe, ods_put1pe made public
!
!EOP
!-------------------------------------------------------------------------


!  Overloaded Interfaces
!  ---------------------
   Interface ODS_Get
       module procedure ODS_Get1_
       module procedure ODS_GetM_
   end Interface

   Interface Init
       module procedure ODS_Init
   end Interface

   Interface Clean
       module procedure ODS_Clean
   end Interface

   Interface Ghost
       module procedure ODS_Ghost
   end Interface

   Interface Unghost
       module procedure ODS_Unghost
   end Interface

   Interface Resize
       module procedure ODS_Resize
   end Interface

   Interface Refresh
       module procedure ODS_Refresh
   end Interface

   Interface Put
       module procedure ODS_Put
   end Interface

   Interface Get
       module procedure ODS_Get1_
       module procedure ODS_GetM_
   end Interface

   Interface Tally
       module procedure ODS_Tally
   end Interface

   Interface ODS_Merge
       module procedure Merge ! NOTE: this odd interface is to cope with
                              !       IRIX64 compiler glitch
   end Interface

   Interface Select
       module procedure ODS_Select
   end Interface

!  ODS vector
!  ----------
   type  ods_vect
      real     :: lat_max              ! Maximum latitude of background grid
      real     :: lat_min              ! Minmum latitude of background grid
      type(decomptype) decomp          ! Fundamental decomposition (excluding halo)
      type(ghosttype) halo_decomp      ! Used if halo is defined
      type(transfertype) halo_trans    ! Transfer pattern; used if halo is defined
      type(ods_meta) meta              ! ODS file metadata (kx_names, etc.)
      type(obs_vect) data              ! Data vectors (lat, lon, obs, omf, etc.)
   end type ods_vect

   integer, external ::  ods_julian
   integer, external ::  ods_caldat
   integer, external ::  ods_parmi


CONTAINS

!-------------------------------------------------------------------------
!         NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE:  ODS_Init --- Allocate memory for ODS vector.
!
! !INTERFACE:
!
  subroutine ODS_Init ( ods, nobs, rc,       &
                        nkt, nkx,  nqc, ncr, nsyn, dist, lmin, lmax ) ! Optional

! !USES:
  use m_MPodsmeta, only : odsm_init
  use m_MPodsdata, only : odsd_init
  use decompmodule, only :  decompcreate

  implicit none

! !INPUT/OUTPUT PARAMETERS:

  type(ods_vect), intent(out)   ::  ods    ! ODS vector to be allocated

  integer, intent(inout)        :: nobs    ! number of observations; if
                                           ! nobs<0  nobs will be reset to
                                           !         a large internal value.
                                           ! does not include halo obs !!!

                                           ! Size of tables and atributes:
  integer, intent(inout), optional :: nkt  !  number of KT's 
  integer, intent(inout), optional :: nkx  !  number of KX's 
  integer, intent(inout), optional :: nqc  !  number of QC flags 
  integer, intent(inout), optional :: ncr  !  number of chars in ODS tables
  integer, intent(inout), optional :: nsyn !  number of synoptic times per
                                           !    day
  integer, intent(in), optional :: dist(:) ! distribution on all PEs
  real,    intent(in), optional :: lmin    ! latitude minimum
  real,    intent(in), optional :: lmax    ! latitude minimum

! !OUTPUT PARAMETERS:

  integer, intent(out)          ::  rc     ! Error return code:
                                           !  0 - all is well
                                           !  1 - error allocating metadata
                                           !  2 - error allocating attributes

! !DESCRIPTION:
! \label{MODS:Init}
!  Allocates memory for an ODS vector.
!
! !REVISION HISTORY:
!
!  16Nov1998   Todling    Initial code
!  02Sep1999   da Silva   Revised names, removed debris.
!  13Sep2000   da Silva   Fixed bug in optional parameters
!  13Dec2000   Redder     Added optional input argument, nsyn
!  16Feb2001   Todling    Changed to account for Xvec addition
!  10Jun2002   Zaslavsky  Split into ODSM_Init, ODSD_Init and ODS_Init
!  18Jun2002   Zaslavsky  Made changes related to relocation of nsyn from 
!                         metadata to data
!  20Dec2002   Sawyer     Moved decompcreate here, use 'dist' variable
!
!EOP
!-------------------------------------------------------------------------

  character(len=*), parameter :: myname = 'ods_init'
  integer ios, nvct, npes

  call ODSM_Init (ods%meta, rc, nkt, nkx, nqc, ncr )
  if ( rc .eq. 0 ) then 
     call ODSD_Init (ods%data, nobs, rc, nsyn )
  endif

  if ( present(dist) ) then
     npes = size(dist)
     call decompcreate(npes,dist,ods%decomp)
  endif

  if ( present(lmin) ) then
     ods%lat_min = lmin
  else
     ods%lat_min = -90.0
  endif

  if ( present(lmax) ) then
     ods%lat_max = lmax
  else
     ods%lat_max = 90.0
  endif

  return

  end subroutine ods_init



!-------------------------------------------------------------------------
!         NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE:  ODS_Clean - Deallocates ODS vector
!
! !INTERFACE:
!
      subroutine ODS_Clean ( ods, rc )

! !USES:
      use m_MPodsmeta, only          : odsm_clean
      use m_MPodsdata, only          : odsd_clean
      use m_odstransfer, only      : transfer_clean
      use decompmodule, only       : decompdefined, decompfree
      use ghostmodule, only        : ghostdefined, ghostfree
      implicit none

! !INPUT/OUTPUT PARAMETERS:

  type(ods_vect), intent(inout) ::  ods    ! ODS vector to be allocated

! !OUTPUT PARAMETERS:

  integer, intent(out)          ::  rc     ! Error return code:
                                           !  0 - all is well
                                           !  1 - error deallocating metadata
                                           !  2 - error deallocating attributes

! !DESCRIPTION:
! 
!  Frees memory used by an ODS vector.
!
! !REVISION HISTORY:
!
!  07sep1999   da Silva   Initial code
!  16Feb2001   Todling    Changed to account for Xvec addition
!  28Mar2001   Todling    Made it go through whole clean even in error.
!  10jun2002   Zaslavsky  Split into ODSM_Clean, ODSD_Clean and ODS_Clean
!  29Oct2002   Sawyer     Added halo decomposition and comm. pattern
!
!EOP
!-------------------------------------------------------------------------

  character(len=*), parameter :: myname = 'ods_clean'
  integer ios

  rc = 0

  call ODSM_Clean (ods%meta, rc )
  if ( rc /= 0 ) return
  call ODSD_Clean (ods%data, rc )
  if ( rc /= 0 ) return
  if ( decompdefined( ods%decomp ) ) call decompfree( ods%decomp )

  if ( ghostdefined( ods%halo_decomp ) ) then
    call ghostfree( ods%halo_decomp )
    call transfer_clean( ods%halo_trans, rc )
  endif

  return

  end subroutine  ods_clean



!-------------------------------------------------------------------------
!         NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE:  ODS_Ghost --- Extend ODS with ghost region
!
! !INTERFACE:
!
  subroutine ODS_Ghost( ods, lat_south, lat_north, rc, no_transfer )

! !USES:
  use parutilitiesmodule, only : gid, gsize, SUMOP, commglobal,            &
                                 parcollective, parexchangevector
  use decompmodule, only : decompglobaltolocal
  use ghostmodule, only :  ghostcreate, ghostdefined
  use m_odstransfer, only : Transfer_Init, Transfer_Perform
  use m_MPodsdata, only : odsd_resize, odsd_sort
  implicit none

! !INPUT PARAMETERS:
  real, intent(in)              :: lat_south ! latitude extent to south
  real, intent(in)              :: lat_north ! latitude extent to north
  logical, optional, intent(in) :: no_transfer ! if present, do not transfer
  
! !INPUT/OUTPUT PARAMETERS:
  type(ods_vect), intent(inout) ::  ods     ! Previous ODS vector

! !OUTPUT PARAMETERS:
  integer, intent(out)          ::  rc     ! Error return code:
                                           !  0 - all is well
                                           !  1 - error in enlargement

! !DESCRIPTION:
! \label{MODS:Ghost}
!  Extends an ODS vector with a ghost region.  The number of
!  inherent observations and the fundamental decomposition 
!  remains the same. By default it updates the ghost 
!  region.  This can be turned off with {\tt no\_transfer}.
!
!  If the ghost zone is to be redefined, it is necessary to
!  call {\tt ODS\_Unghost} to clear the current zone.
!
! !REVISION HISTORY:
!
!  20Nov2002   Sawyer        Creation (prototype for discussion)
!  05Dec2002   Sawyer        Bug fixes; testing; becomes baseline.
!
!EOP
!-------------------------------------------------------------------------

  real, allocatable :: lat_bounds(:,:)
  integer           :: ipe, nobs_inh, nhalo, ntot, i, cnt, oldcnt
  integer, allocatable :: ind_in(:), ind_out(:), tags(:)
  integer, allocatable :: ind_in_cnt(:), ind_out_cnt(:)
  character(len=*), parameter :: myname = 'ods_ghost'
  integer           :: loc, pe, nsize

  nobs_inh = ods%data%nobs

  rc = 0

!
! Check to see if the ghost descriptor is already defined
! -------------------------------------------------------
  if ( .not.ghostdefined(ods%halo_decomp) ) then
     nsize = 2*int((lat_north+lat_south)/(ods%lat_max-ods%lat_min)+2.0)*nobs_inh
     allocate( lat_bounds(gsize,2) )
     allocate( ind_in_cnt(gsize), ind_out_cnt(gsize) )
     allocate( ind_in(nsize), ind_out(nsize) )

     lat_bounds = 0.0
     lat_bounds(gid+1,1) = ods%lat_max
     lat_bounds(gid+1,2) = ods%lat_min
     call parcollective( commglobal, SUMOP, gsize, 2, lat_bounds )

!
! Find the observations needed for the PE to the north and south
! --------------------------------------------------------------
     cnt=0
     oldcnt = 0
     do ipe = 1, gsize
        if ( ipe /= gid+1 ) then
          do i=1, nobs_inh
            if ( (lat_bounds(ipe,1)+lat_north >= ods%data%lat(i)) .and. &
              (ods%data%lat(i) >= lat_bounds(ipe,1)) ) then  ! North zone
              cnt = cnt+1
              ind_in(cnt) = ods%data%gid(i)
            endif
          enddo
          do i=1, nobs_inh
            if ( (lat_bounds(ipe,2)-lat_south <= ods%data%lat(i)) .and. &
              (ods%data%lat(i) <= lat_bounds(ipe,2)) ) then  ! South zone
              cnt = cnt+1
              ind_in(cnt) = ods%data%gid(i)
            endif
          enddo
        endif
        ind_in_cnt(ipe) = cnt-oldcnt
        oldcnt = cnt
     enddo

     print *, "gid", gid, "cnt=", cnt, "should be <=", nsize

!
! Send the indices (or "tags") to the requesting PEs
! --------------------------------------------------
     call parexchangevector( commglobal, ind_in_cnt, ind_in,             &
                                         ind_out_cnt, ind_out ) 

     nhalo = sum( ind_out_cnt )
     ntot = nobs_inh + nhalo

     allocate(tags(ntot))
     do i=1, nobs_inh
       tags(i) = ods%data%gid(i)
     enddo
     do i=1, nhalo
       tags(nobs_inh+i) = ind_out(i)
     enddo

!
! Create a ghost decomposition from the overlapping points, added
! additional space to the data regions, form the transfer pattern
! ---------------------------------------------------------------
     call ghostcreate( ods%decomp, gid, ntot, tags, ods%halo_decomp )
     deallocate(tags)

     call odsd_resize( ods%data, ntot, .true., rc )
     if ( rc /= 0 ) return

     call transfer_init(ods%halo_decomp, ods%halo_trans, rc )
     if ( rc /= 0 ) return

     deallocate( ind_in, ind_out )
     deallocate( ind_in_cnt, ind_out_cnt )
     deallocate( lat_bounds )
  endif

! Fill in ghost region unless no_transfer is present

  if ( .not. present( no_transfer ) ) then
    call transfer_perform( ods%halo_trans, ods%data, ods%data, rc )
  endif
  return

  end subroutine ODS_Ghost

!-------------------------------------------------------------------------
!         NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE:  ODS_Unghost --- Remove the ghost region from ods vector
!
! !INTERFACE:
!
  subroutine ODS_Unghost( ods, rc )

! !USES:
  use parutilitiesmodule, only : gid
  use ghostmodule, only :  ghostfree, ghostdefined
  use m_odstransfer, only : Transfer_Clean
  use m_MPodsdata, only : odsd_resize
  implicit none

! !INPUT/OUTPUT PARAMETERS:
  type(ods_vect), intent(inout) ::  ods     ! Previous ODS vector

! !OUTPUT PARAMETERS:
  integer, intent(out)          ::  rc     ! Error return code:
                                           !  0 - all is well
                                           !  1 - error in data reduction
                                           !  2 - halo was not defined

! !DESCRIPTION:
! \label{MODS:Unghost}
!  Removes the ghost region from the ODS vector.  Slims down
!  vector to size "nobs".  Observations outside 1..nobs are lost.
!  But the decomposition remains the same.
!
! !REVISION HISTORY:
!
!  17Jul2002   Sawyer     Initial code, from ODS_Init
!  26Oct2002   Sawyer     Now uses ODSD_Resize
!
!EOP
!-------------------------------------------------------------------------

  character(len=*), parameter :: myname = 'ods_unghost'

  call odsd_resize( ods%data, ods%data%nobs, .true., rc )
  if ( rc /= 0 ) return

  rc = 0
  if ( ghostdefined( ods%halo_decomp ) ) then
!
! Destroy the ghost decomposition and the communication pattern
! -------------------------------------------------------------
    call ghostfree( ods%halo_decomp )
    call transfer_clean( ods%halo_trans, rc )
  else
    rc = 2
  endif

  return

  end subroutine ods_unghost

!-------------------------------------------------------------------------
!         NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE:  ODS_Resize --- Change the number of observations
!
! !INTERFACE:
!
  subroutine ODS_Resize( ods, nobs_new, rc )

! !USES:
  use parutilitiesmodule, only : gid, gsize, commglobal, parcollective, SUMOP
  use decompmodule, only :  decompcreate, decompfree
  use m_MPodsdata, only : odsd_resize
  implicit none

! !INPUT PARAMETERS:

  integer, intent(in)           ::  nobs_new ! new number of obs

! !INPUT/OUTPUT PARAMETERS:
  type(ods_vect), intent(inout) ::  ods     ! Previous ODS vector

! !OUTPUT PARAMETERS:
  integer, intent(out)          ::  rc     ! Error return code:
                                           !  0 - all is well
                                           !  1 - error in enlargement

! !DESCRIPTION:
! \label{MODS:Resize}
!  Extends an ODS vector with extra space which is {\em not}
!  for a ghost region.  This calls ODSD\_Resize, then updates
!  the decomposition.  This functionality may be useful in the
!  future for super-obbing.
!
! !REVISION HISTORY:
!
!  31Oct2002   Sawyer     Initial code, from ODS_Ghost
!
!EOP
!-------------------------------------------------------------------------

  character(len=*), parameter :: myname = 'ods_resize'
  integer, allocatable :: num_in_lat(:)

  call decompfree( ods%decomp )
  allocate( num_in_lat( gsize ) )
  num_in_lat = 0
  num_in_lat( gid+1 ) = nobs_new
  call parcollective( commglobal, SUMOP, gsize, num_in_lat )
  call decompcreate( gsize, num_in_lat, ods%decomp )
  deallocate( num_in_lat )
  call odsd_resize( ods%data, nobs_new, .false., rc )  ! no halo region

  return

  end subroutine ods_resize


!-------------------------------------------------------------------------
!         NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE:  ODS_Refresh --- Force the ODS to be decomposition consistent
!
! !INTERFACE:
!
  subroutine ODS_Refresh( ods, rc )

! !USES:
  use parutilitiesmodule, only : gid, gsize, commglobal, parcollective, SUMOP
  use decompmodule, only :  decompcreate, decompfree
  implicit none

! !INPUT/OUTPUT PARAMETERS:
  type(ods_vect), intent(inout) ::  ods     ! Previous ODS vector

! !OUTPUT PARAMETERS:
  integer, intent(out)          ::  rc     ! Error return code:
                                           !  0 - all is well

! !DESCRIPTION:
! \label{MODS:Refresh}
!   This routine refreshes the decomposition of an ODS vector with
!   the current number of observations.  It does not test if the
!   the decomposition is already consistent.  It is advised not to
!   use this routine to ``clean'' things up -- if the decomposition
!   is inconsistent it is an indication that there is a bug in the
!   code.
!
! !REVISION HISTORY:
!
!  31Oct2002   Sawyer     Initial code, from ODS_Ghost
!
!EOP
!-------------------------------------------------------------------------

  character(len=*), parameter :: myname = 'ods_refresh'
  integer, allocatable :: num_in_lat(:)

  call decompfree( ods%decomp )
  allocate( num_in_lat( gsize ) )
  num_in_lat = 0
  num_in_lat( gid+1 ) = ods%data%nobs
  call parcollective( commglobal, SUMOP, gsize, num_in_lat )
  call decompcreate( gsize, num_in_lat, ods%decomp )
  deallocate( num_in_lat )
  rc = 0
  return

  end subroutine ods_refresh

!-------------------------------------------------------------------------
!         NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE:  ODS_GetM_ --- Reads data from an ODS file
! 
! !INTERFACE:

  subroutine ODS_GetM_ ( nfiles, fnames, nymd, nhms, ftypes, ods, rc,      &
                         iymd, ihms, lat_max, contents )   ! Optional

!
! !USES
!
    use parutilitiesmodule, only : gid, gsize, commglobal, BCSTOP,         &
                                   parcollective
    Implicit NONE

! !INPUT PARAMETERS: 
!
    integer, intent(in)            :: nfiles  ! number of files
    character(len=*), intent(in)   :: fnames(nfiles) ! ODS file names
    integer, intent(in)            :: nymd    ! year-month-day, e.g., 19990701
    integer, intent(in)            :: nhms    ! hour-min-sec,   e.g., 120000
    integer, intent(in), optional  :: iymd(nfiles)   ! year-month-day for 
                                                     ! every file
    integer, intent(in), optional  :: ihms(nfiles)   ! hour-min-sec for 
                                                     ! every file
    real, intent(inout), optional  :: lat_max(:)     ! bidirectional
                                              ! max lat for each PE
                                              ! if = 0, define and return
                                              ! "optimal" values (see below)
    character(len=*), intent(in), optional :: contents
                                              ! type of file contents
                                              ! 'ods' (default if not present)
                                              ! 'oms' (skin temperature only)
                                              ! 'oms++' (the works, not impl.)
                                              ! 'buffer' (not yet implemented)
! !OUTPUT PARAMETERS:
!
    character(len=*), intent(out)   :: ftypes(nfiles) 
                                              ! ODS file type: 'pre_anal'
                                              !  or 'post_anal'. If 'unknown'
                                              !  the could not be read

    type(ods_vect), intent(inout)   :: ods    ! ODS vector

    integer, intent(out)           ::  rc     ! Error return code:
                                              ! = 0    - all is well
                                              ! <1000  - could not read file # rc
                                              ! > 999  - rc-999 is rc from ODS_Merge 

! !DESCRIPTION:
! \label{MODS:GetM}
!  This routine reads data from 1 or more ODS files, allocates the
!  necessary memory for the ODS vector, and loads the data for the
!  synoptic time (nymd,nhms). Check the contents of {\tt ftypes} to
!  check whether individual files could be read sucessfully. 
!
!  {\bf SPMD Notes:}  The decomposition of the ODS
!    vector can be specified by an array {\tt lat\_max} containing
!    the maximum latitude an observation may have on a given
!    PE.  If {\tt lat\_max} is absent, the observations are sorted
!    by latitude and an approximately equal number handed out to
!    each PE.  Note this implicitly creates a latitude strip
!    decomposition.
!
! !REVISION HISTORY: 
!
!  03Dec1999  da Silva  Initial code.
!  09Mar2000  da Silva  Now is OK if some of the files could not be read,
!                       the only requiment is that at least one file can
!                       read. 
!  12Dec2999  Redder    Revised prologue to reflect the added flexibility
!                       of arbitrarily setting the ODS file parameter,
!                       nsyn.
!  18Jun2002  Zaslavsky Allowed multiple synoptic date and times for 
!                       multiple files, iymd(:), ihms(:)
!  26Oct2002  Sawyer    Upgraded for SPMD
!  05Dec2002  Sawyer    Separated out ODS_Ghost (removed dl argument)
!  30Jan2003  Sawyer    Added ROMS capability (skin temperature only)
!  12Feb2003  Sawyer    Rewritten to read/merge on 1 PE, then scatter
!
!EOP
!-------------------------------------------------------------------------
! !LOCAL VARIABLES
  type (ods_vect)      :: ods_g
  integer              :: i, ier, n
  type(ods_vect)       :: tods(nfiles)   ! to hold data for each file

  character(len=*), parameter :: myname = 'ods_getm_'

! Start of code

  if ( gid == 0 ) then       ! On 1 PE
     n = 1
     do i=1, nfiles
        if ( .not. present(contents) ) then
           if ( present( iymd ) .and. present( ihms ) ) then
              call ODS_Get1pe ( fnames(i), iymd(i), ihms(i), ftypes(i), tods(n), rc )
           else
              call ODS_Get1pe ( fnames(i), nymd, nhms, ftypes(i), tods(n), rc )
           endif
        else
           select case (contents)
           case ('ods')
              if ( present( iymd ) .and. present( ihms ) ) then
                 call ODS_Get1pe ( fnames(i), iymd(i), ihms(i), ftypes(i), tods(n), rc )
              else
                 call ODS_Get1pe ( fnames(i), nymd, nhms, ftypes(i), tods(n), rc )
              endif
           case ('oms') 
              if ( present( iymd ) .and. present( ihms ) ) then
                 call ROMS_Get1pe_ ( fnames(i), iymd(i), ihms(i), ftypes(i), tods(n), rc, &
                                     only_ts = .true. )
              else
                 call ROMS_Get1pe_ ( fnames(i), nymd, nhms, ftypes(i),   &
                                     tods(n), rc, only_ts = .true. )
              endif
           case ('oms++')          ! not yet supported
              if ( present( iymd ) .and. present( ihms ) ) then
                 call ROMS_Get1pe_ ( fnames(i), iymd(i), ihms(i), ftypes(i), &
                                     tods(n), rc )
              else
                 call ROMS_Get1pe_ ( fnames(i), nymd, nhms, ftypes(i),       &
                                     tods(n), rc )
              endif
           case ('buffer')         ! not yet supported
              rc = 6
           case default            ! other values not supported
              rc = 5
           end select
        endif
        if ( rc .eq. 0 ) then
           n = n + 1
        else
           ftypes(i) = 'unknown'
        end if
     enddo
     n = n - 1   ! number of files successfully read


! Merge ODS files
! ---------------
     call ODS_Merge ( tods, n, ods_g, ier, single=.true. )
     if ( ier .ne. 0 ) then
        rc = 999 + ier
     end if


! Free memory used to hold data for individual files
! --------------------------------------------------
     do i = 1, n
        call ODS_Clean ( tods(i), ier )  ! ignore ier
     end do
  endif

  do i = 1,nfiles
     call parcollective( commglobal, BCSTOP, ftypes(i) )
  enddo

! Partition the global ODS vector over all PEs
! --------------------------------------------
  if ( present( lat_max ) ) then
     call ods_partition_ (ods_g, ods, ier, lat_max )
  else
     call ods_partition_ (ods_g, ods, ier )
  endif

  if ( ier .ne. 0 ) then
       rc = 999 + ier
       return
  else
       rc = 0
  end if


! All done
! --------
  return

end subroutine ODS_GetM_


!-------------------------------------------------------------------------
!         NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE:  ODS_Get1_ --- Reads data from a single ODS file
! 
! !INTERFACE:

  subroutine ODS_Get1_ ( fname, nymd, nhms, ftype, ods, rc,             &
                         ncid, lat_max, contents )        ! Optional

!
! !USES
!
    use parutilitiesmodule, only : gid, gsize, commglobal, BCSTOP,      &
                                   parcollective
    Implicit NONE

! !INPUT PARAMETERS: 
!
    character(len=*), intent(in)   :: fname   ! ODS file name
    integer, intent(in)            :: nymd    ! year-month-day, e.g., 19990701
    integer, intent(in)            :: nhms    ! hour-min-sec,   e.g., 120000

    integer, intent(in), optional  :: ncid    ! ODS file id as returned by
                                              ! ODS_Open(); in this case
                                              ! fname is not used
    real, intent(inout), optional  :: lat_max(:) ! max lat for each PE
                                              ! if = 0, define and return
                                              ! "optimal" values (see below)
    character(len=*), intent(in), optional :: contents
                                              ! type of file contents
                                              ! 'ods' (default if not present)
                                              ! 'oms' (skin temperature only)
                                              ! 'oms++' (the works, not impl)
                                              ! 'buffer'

! !OUTPUT PARAMETERS:
!
    character(len=*), intent(out)  :: ftype   ! ODS file type: 'pre_anal'
                                              !  or 'post_anal'

    type(ods_vect), intent(inout)  :: ods     ! ODS vector

    integer, intent(out)           :: rc      ! Error return code:
                                              !  0 - all is well
                                              !  1 - could not open file
                                              !  2 - could not get table sizes
                                              !  3 - could not get no. obs
                                              !  4 - could not allocate memory
                                              !  5 - could not read tables
                                              !  6 - could not read attributes (pre)
                                              !  7 - could not read attributes (post)
                                              !  8 - could not obtain first julian day on file

! !DESCRIPTION:
! \label{MODS:Get1_}
!    Get an ODS vector from a file.  The decomposition of the ODS
!    vector can be specified by an array {\tt lat\_max} containing
!    the maximum latitude an observation may have on a given
!    PE.  If {\tt lat\_max} is absent, the observations are sorted
!    by latitude and an approximately equal number handed out to
!    each PE.  Note this implicitly creates a latitude strip
!    decomposition.
!
! !REVISION HISTORY: 
!
!  10May2002  Sawyer     Initial code, SPMD wrapper for ODS_Get
!  06Jun2002  Sawyer     Made into standard ODS_Get1_
!  17Sep2002  Sawyer     Added logic for ghost region
!  31Oct2002  Sawyer     Extensive changes to ghosting mechanism
!  05Dec2002  Sawyer     Separated out ODS_Ghost (removed dl argument)
!  20Dec2002  Sawyer     Moved decompcreate out, into ods_init
!  30Jan2002  Sawyer     Added ROMS capability (skin temperature only!)
!  12Feb2003  Sawyer     Moved partition/scatter to ODS_Partition_
!
!EOP
!-------------------------------------------------------------------------

! !LOCAL VARIABLES
  type (ods_vect)      :: ods_g

!
! PE 0 reads the ODS file
! -----------------------     
     if ( gid == 0 ) then
        if ( .not. present(contents) ) then
           if ( present( ncid ) ) then
              call Ods_get1pe ( fname, nymd, nhms, ftype, ods_g, rc, ncid )
           else
              call Ods_get1pe ( fname, nymd, nhms, ftype, ods_g, rc )
           endif
        else
           select case (contents)
           case ('ods')
              if ( present( ncid ) ) then
                 call Ods_get1pe ( fname, nymd, nhms, ftype, ods_g, rc, ncid )
              else
                 call Ods_get1pe ( fname, nymd, nhms, ftype, ods_g, rc )
              endif
           case ('oms') 
              call ROMS_Get1pe_ ( fname, nymd, nhms, ftype, ods_g, rc, &
                                  only_ts = .true. )
           case ('oms++')        ! not yet supported
              call ROMS_Get1pe_ ( fname, nymd, nhms, ftype, ods_g, rc )
           case ('buffer')       ! not yet supported
              rc = 6
              return
           case default          ! others not allowed
              rc = 5
              return
           end select
        endif
     end if 

     if ( present( lat_max ) ) then
        call ods_partition_( ods_g, ods, rc, lat_max )
     else 
        call ods_partition_( ods_g, ods, rc )
     end if
     call parcollective( commglobal, BCSTOP, ftype )

  return
  end subroutine ODS_Get1_

!-------------------------------------------------------------------------
!         NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE:  ODS_Put --- Writes to an ODS file
! 
! !INTERFACE:

  subroutine ODS_Put ( fname, ftype,  nymd, nhms, ods, rc, &
                           ncid,  append, new )         ! Optional
!
! !USES:
!
    use m_MPodsdata, only : odsd_sort
    use parutilitiesmodule, only : gid, gsize, commglobal, SUMOP, &
                                   parcollective
    use m_odstransfer, only : TransferType, Transfer_Init, Transfer_Clean
    Implicit NONE

! !INPUT PARAMETERS: 
!
    character(len=*), intent(in)   :: fname   ! ODS file name
    character(len=*), intent(in)   :: ftype   ! ODS file type: 'pre_anal'
                                              !  or 'post_anal'
    integer, intent(in)            :: nymd    ! year-month-day, e.g., 19990701
    integer, intent(in)            :: nhms    ! hour-min-sec,   e.g., 120000

    integer, intent(in), optional  :: ncid    ! ODS file id as returned by
                                              ! ODS_Open(); in this case
                                              ! fname is not used

    type(ods_vect), intent(inout)  :: ods     ! ODS vector

    logical, intent(in), optional  :: append  ! If true, the data is appended
                                              ! to the ODS file at this time.
                                              ! Default: append = .false.

    logical, intent(in), optional  :: new     ! If true, file will be created
                                              ! whether it exists or not
                                              ! Default: new    = .false.

!
! !OUTPUT PARAMETERS:
!
    integer, intent(out)           :: rc      ! Error return code:
                                              !  0 - all is well
                                              !  1 - could not open file
                                              !  2 - could not create file
                                              !  3 - could not append file
                                              !  4 - could not write a
                                              !     'pre_anal' observation
                                              !      attribute
                                              !  5 - could not write a
                                              !     'post_anal' observation
                                              !      attribute
                                              !  6 - error in allocating
                                              !      scratch space.
                                              !  7 - could not obtain the
                                              !      number of observations
                                              !      already on file
                                              !  8 - could not obtain the
                                              !      the first Julian day
                                              !      on file.
!
! !DESCRIPTION:
! \label{MODS:Put}
!   This routine is an SPMD wrapper for {\tt ODS\_Put1pe\_}.  If the 
!   optional argument {\tt lat\_max} is is left off, it assumes the 
!   {\tt ods\_vect} is defined in a global manner on PE 0.  If 
!   {\tt lat\_max} exists, it gathers the ods vector on PE zero 
!   with a transfer.
!
! !REVISION HISTORY: 
!
!  10May2002  Sawyer  Initial code.
!  18Jun2002 Zaslavsky Made changes related to relocation of nsyn from 
!                      metadata to data
!
!EOP
!-------------------------------------------------------------------------

   integer nobs, id
   integer, allocatable :: num_in_lat(:)
   type (ods_vect) :: ods_g
   type (transfertype) :: gather

      nobs = ods%data%nobs
      call parcollective( commglobal, SUMOP, nobs )
      allocate( num_in_lat( gsize ) )
      num_in_lat = 0
      num_in_lat(1) = nobs
      call ods_init  ( ods_g, nobs, rc,                                 &
                       nkt=ods%meta%nkt, nkx=ods%meta%nkx,              &
                       nqc=ods%meta%nqc, ncr=ods%meta%ncr,              &
                       nsyn=ods%data%nsyn, dist=num_in_lat )
      deallocate( num_in_lat )
      if ( rc .ne. 0 ) then
         print *, 'Error initializing ods', 'rc = ', rc
      end if


!
! Gather the ODS vector
! ---------------------
      call Transfer_Init( ods%decomp, ods_g%decomp, gather, rc )
      call ODS_Transfer( gather, ods, ods_g, rc )
      call Transfer_clean( gather, rc )

      if ( gid == 0 ) then

         call odsd_sort( ods_g%data, nobs, rc, 'gid' )

         if ( present( ncid ) ) then
            if ( present( append ) ) then
               if (present( new ) ) then 
                  call ods_put1pe( fname, ftype,  nymd, nhms, ods_g, rc, &
                                    ncid=ncid, append=append, new=new )
               else
                  call ods_put1pe( fname, ftype,  nymd, nhms, ods_g, rc, &
                                    ncid=ncid, append=append )
               endif
            else
               if (present( new ) ) then 
                  call ods_put1pe( fname, ftype,  nymd, nhms, ods_g, rc, &
                                    ncid=ncid, new=new )
               else
                  call ods_put1pe( fname, ftype,  nymd, nhms, ods_g, rc, ncid=ncid )
               endif
            endif
         else 
            if ( present( append ) ) then
               if (present( new ) ) then 
                  call ods_put1pe( fname, ftype,  nymd, nhms, ods_g, rc, &
                                    append=append, new=new )
               else
                  call ods_put1pe( fname, ftype,  nymd, nhms, ods_g, rc, &
                                    append=append )
               endif
            else
               if (present( new ) ) then 
                  call ods_put1pe( fname, ftype,  nymd, nhms, ods_g, rc, &
                                    new=new )
               else
                  call ods_put1pe( fname, ftype,  nymd, nhms, ods_g, rc )
               endif
            endif
         endif
      endif
      call ods_clean( ods_g, rc )

   return
   end subroutine ODS_Put

!-------------------------------------------------------------------------
!         NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE:  ODS_Merge --- Merges 2 or more ODS vectors 
!
! !INTERFACE:
!
      subroutine Merge ( ods_in, nods, ods_out, rc, nymd, nhms, single )

! !USES:
      use parutilitiesmodule, only : commglobal, gid, gsize, MAXOP,       &
                                     parcollective
      implicit none

! !INPUT PARAMETERS:

  integer,        intent(in)   ::  nods          ! number of input ODS vectors

  type(ods_vect), intent(in)   ::  ods_in(nods)  ! Many ODS vectors
  integer, optional, intent(in)::  nymd, nhms    ! synoptic date and time for
                                                 ! the merged file
  logical, optional, intent(in)::  single        ! if single task

! !OUTPUT PARAMETERS:

  type(ods_vect), intent(inout) ::  ods_out       ! Merge ODS vector
  integer,        intent(out)   ::  rc            ! Error return code:
                                                  !  0 - all is well
                                                  !  1 - could not allocate mem

! !DESCRIPTION:
! \label{MODS:Merge}
!    Merges 2 or more ODS vectors.  This supports both single PE and SPMD
!    modes, however, in the latter, it is assumed that all observations
!    on a given PE are in the same latitude band, i.e., that 
!    \verb$ods%lat_max$ and \verb$ods%lat_min$ are the same on all PEs.
!
! !BUGS:
!
!     LIMITATIONS:  The merged vector has the same metadata 
!     as the first vector, that is the metadata is not really merged. 
!     Furthermore: it is assumed that all ODS vectors have the same
!     latitude minima and maxima (ods%lat_max, ods%lat_min).  
!
!     The sounding index ks of the resulting ODS vector is unique
!     if the input ODS vectors are unique.  Finally, the "single"
!     argument for 1 PE mode (needed in ODS_GetM_) is a hack.
!
! !REVISION HISTORY:
!
!  03dec1999   da Silva   Initial code
!  23dec1999   da Silva   Fixed metadata size bug.
!  26jul2000   da  Silva  Removed kid from merging since kid is set
!                         to 1:nvct during ods_init().
!  16Feb2001   Todling    Changed to account for Xvec addition
!  04Oct2001   G Lou      Added nsyn for RUC
!  11Jun2002   Sawyer     SPMD version
!  18Jun2002   Zaslavsky  Made changes related to relocation of nsyn from 
!                         metadata to data.
!  21Jun2002   Zaslavsky  Included time adjustment.
!  20Dec2002   Sawyer     Initial components for SPMD implementation
!  05Feb2003   Sawyer     Added initialization of data%gid (bug fix)
!  12Feb2003   Sawyer     Added "single" pe mode for ods_getm_
!
!EOP
!-------------------------------------------------------------------------

  character(len=*), parameter :: myname = 'ods_merge'
  
  integer i, ier, n, n1, n2, ks_offset, tmpmaxval
  integer nkt, nkx, nqc, ncr, nobs, nnsyn, nnymd, nnhms, mymd, mhms, k
  integer, allocatable :: num_in_lat(:)
  real :: lmin, lmax, lminmax, lmaxmin

  rc = 0

! Get total number of observations
! --------------------------------
  nobs = 0
  do i = 1, nods
     nobs = nobs + ods_in(i)%data%nobs
  end do

  allocate( num_in_lat(gsize) )
  num_in_lat = 0
  num_in_lat(gid+1) = nobs
  if ( .not. single ) then
    call parcollective( commglobal, MAXOP, gsize, num_in_lat )
  endif

! Metadata sizes
! --------------
  nkt=maxval(ods_in(1:nods)%meta%nkt) 
  nkx=maxval(ods_in(1:nods)%meta%nkx)
  nqc=maxval(ods_in(1:nods)%meta%nqc)
  ncr=maxval(ods_in(1:nods)%meta%ncr)
  lmax=maxval(ods_in(1:nods)%lat_max)
  lmaxmin=minval(ods_in(1:nods)%lat_max)
  lminmax=maxval(ods_in(1:nods)%lat_min)
  lmin=minval(ods_in(1:nods)%lat_min)

! Allocate memory for merged ODS file
! -----------------------------------
  nnsyn = ods_in(1)%data%nsyn
  call ODS_Clean ( ods_out, ier )      ! ignore ier
  call ODS_Init  ( ods_out, nobs, ier, nkt, nkx, nqc, ncr, nnsyn, &
                   dist=num_in_lat, lmin=lmin, lmax=lmax )

! Initialize data%gid  these tags MUST BE UNIQUE
! ----------------------------------------------
  do i = 2,gsize
    num_in_lat(i) = num_in_lat(i) + num_in_lat(i-1)  ! offsets
  enddo
  if (gid > 0 ) then
    ods_out%data%gid = (/ (i, i=num_in_lat(gid)+1,num_in_lat(gid+1)) /)
  else
    ods_out%data%gid = (/ (i, i=1,num_in_lat(1)) /)
  endif
  deallocate( num_in_lat )

  if ( ier .ne. 0 ) then
       rc = 1
       return
  end if

! Set merged metadata from first ODS 
! ----------------------------------
  ods_out%meta%kt_names  = ods_in(1)%meta%kt_names
  ods_out%meta%kt_units  = ods_in(1)%meta%kt_units
  ods_out%meta%kx_names  = ods_in(1)%meta%kx_names
  ods_out%meta%kx_meta   = ods_in(1)%meta%kx_meta
  ods_out%meta%qcx_names = ods_in(1)%meta%qcx_names

  if ( nobs .le. 0 ) return

  if( present( nymd ) .and. present ( nhms ) ) then
     nnymd = nymd
     nnhms = nhms
  else
     nnymd = ods_in(1)%data%nymd
     nnhms = ods_in(1)%data%nhms
  endif

! Set merged obs vector
! ---------------------
  n2 = 0
  ks_offset = 0
  do i = 1, nods
   
     n  = ods_in(i)%data%nobs
     if ( n .eq. 0 ) cycle
     n1 = n2 + 1
     n2 = n1 + n - 1

     ods_out%data%lat(n1:n2)    = ods_in(i)%data%lat(1:n)
     ods_out%data%lon(n1:n2)    = ods_in(i)%data%lon(1:n)
     ods_out%data%lev(n1:n2)    = ods_in(i)%data%lev(1:n)
     ods_out%data%kx(n1:n2)     = ods_in(i)%data%kx(1:n)
     ods_out%data%kt(n1:n2)     = ods_in(i)%data%kt(1:n)
     ods_out%data%xm(n1:n2)     = ods_in(i)%data%xm(1:n)
     ods_out%data%time(n1:n2)   = ods_in(i)%data%time(1:n)

!    Adjust time
     mymd = ods_in(i)%data%nymd
     mhms = ods_in(i)%data%nhms
     if ( ( nnymd .ne. mymd ) .or. (nnhms .ne. mhms ) ) then  
         do k = n1, n2
              call ODS_NewRefTime ( mymd, mhms,            &
                                    nnymd, nnhms,          &
                                    ods_out%data%time(k:k) &  
                                  )
         end do
     end if

     ods_out%data%obs(n1:n2)    = ods_in(i)%data%obs(1:n)
     ods_out%data%OmF(n1:n2)    = ods_in(i)%data%OmF(1:n)
     ods_out%data%OmA(n1:n2)    = ods_in(i)%data%OmA(1:n)
     ods_out%data%Xvec(n1:n2)   = ods_in(i)%data%Xvec(1:n)
     ods_out%data%qcexcl(n1:n2) = ods_in(i)%data%qcexcl(1:n)
     ods_out%data%qchist(n1:n2) = ods_in(i)%data%qchist(1:n)

!    Make sure sounding index is unique  (local to PE only)
!    ------------------------------------------------------
     ods_out%data%ks(n1:n2)     = ks_offset + ods_in(i)%data%ks(1:n)
     tmpmaxval = maxval ( ods_in(i)%data%ks(1:n) )   ! local maximum
     if ( .not. present( single ) ) then             ! global maximum
        call parcollective( commglobal, MAXOP, tmpmaxval )
     endif
     ks_offset = ks_offset + tmpmaxval
               
  end do

!
! The inconsistent case where the lat_max/lat_min do not correspond
! is tolerated but produces a non-zero error code:
! ------------------------------------------------

  if ( lmax > lmaxmin .or. lmin < lminmax ) then
     rc = 3
  endif


! All done
! --------
  return

end subroutine Merge



!-------------------------------------------------------------------------
!         NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!-------------------------------------------------------------------------
! BOP
!
! !ROUTINE:      ODS_NewRefTime - Adjust sampling time for new reference date and time
!
! !INTERFACE:
      subroutine ODS_NewRefTime ( nymd1, nhms1,        &
                                  nymd2, nhms2, Time )
!
! !INPUT PARAMETERS:
      implicit NONE
      integer, intent (in) :: nymd1, nhms1 ! Old and ...
      integer, intent (in) :: nymd2, nhms2 ! ... new reference date and time
!
! !OUTPUT PARAMETERS:
      integer, intent (inout), dimension (:) :: &
                              Time         ! Time (in minutes) since the
                                           !   reference date and time
!
! !DESCRIPTION: 
! \label{MODS:NewRefTime}
!     Adjusts the time to account for the change in the reference
!     synoptic date and time.
!
! !REVISION HISTORY: 
!     12Dec2000   Redder   Origional version
!
! EOP
!-------------------------------------------------------------------------

      integer  iVal, NVal, AdjTime, dSec

      AdjTime = ( ODS_Julian ( nymd1 )                    &
              -   ODS_Julian ( nymd2 ) ) * 60 * 24        &
              + ( nhms1 / 10000 - nhms2 / 10000 ) * 60    &
              +   mod ( nhms1, 10000 ) / 100              &
              -   mod ( nhms2, 10000 ) / 100

      dSec    =   mod ( nhms1, 100 ) - mod ( nhms2, 100 )
      if      ( dSec .lt. -30 ) then
              AdjTime = AdjTime - 1
      else if ( dSec .ge.  30 ) then
              AdjTime = AdjTime + 1
      end if

      NVal    = size ( Time )
      do iVal = 1, NVal
         Time ( iVal ) = Time ( iVal ) + AdjTime

      end do

      return
      end subroutine ODS_NewRefTime

!..............................................................


!-------------------------------------------------------------------------
!         NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE:  ODS_Tally - Prints out ODS summary
!
! !INTERFACE:
!

   subroutine ODS_Tally ( lu, ods, nobs, rc, &
                          kt_only )  ! optional

! !USES:

      implicit none

! !INPUT PARAMETERS:

      integer,     intent(in)        :: lu      ! unit number for ASCII
                                                !  OUTPUT
      type(ods_vect), intent(in)     :: ods     ! ODS vector
      integer,     intent(in)        :: nobs    ! no. of obs

      logical, intent(in), OPTIONAL  :: kt_only ! only summary by kt


! !OUTPUT PARAMETERS:
!

      integer,     intent(out)       :: rc      ! error code:
                                                !  0   - all is well
                                                ! ...  - abnormal exit

! !DESCRIPTION:
! \label{MODS:Tally}
!  Prints out a summary of an ODS vector.
!
! !REVISION HISTORY:
!
!  13dec1999   da Silva   Initial code based on OBSMRY() from the PSAS
!                         library.
!
!EOP
!-------------------------------------------------------------------------

   character(len=*), parameter :: myname = 'ods_tally'


!  Local storage for counters
!  --------------------------
   integer, allocatable ::  ksums(:,:,:)
   integer, allocatable ::  ktots(:,:)

   integer ktmax, kxmax, ios, kgrand(2), i, j, k, n, ksum
   
   logical :: kx_summary 



   rc = 0
   if ( nobs .eq. 0 ) return


!  Determine type of summary
!  -------------------------
   kx_summary = .true.
   if ( present(kt_only) ) then
        if ( kt_only ) kx_summary = .false.
   end if

!  Allocate local workspace
!  ------------------------
   kxmax = maxval ( ods%data%kx(1:nobs) )
   ktmax = maxval ( ods%data%kt(1:nobs) )
   if ( kxmax .gt. ods%meta%nkx .or. &
        ktmax .gt. ods%meta%nkt        ) then
      rc = 1
      return
   end if
   allocate ( ksums(kxmax,ktmax,2), ktots(ktmax,2), stat=ios )
   if ( ios .ne. 0 ) then
      rc = 2
      return
   end if

  
!  Initialize counters
!  -------------------
   ktots = 0
   ksums = 0

!  Count up the data by source and type
!  ------------------------------------
   do n = 1, nobs
      ksums(ods%data%kx(n),ods%data%kt(n),2) = ksums(ods%data%kx(n),ods%data%kt(n),2) + 1
      ktots(ods%data%kt(n),2)                = ktots(ods%data%kt(n),2)       + 1
      if ( ods%data%qcexcl(n) .eq. 0 ) then
        ksums(ods%data%kx(n),ods%data%kt(n),1) = ksums(ods%data%kx(n),ods%data%kt(n),1) + 1
        ktots(ods%data%kt(n),1)                = ktots(ods%data%kt(n),1)       + 1
      end if
   end do


!  Write the observation summary sub-tables for each data source
!  -------------------------------------------------------------
   if ( kx_summary ) then
   do i = 1, kxmax

      ksum = 0
      do j = 1, ktmax
         ksum = ksum + ksums(i,j,2)
      end do

      if( ksum.gt.0 ) then

         write(lu,'(a)') myname // ':'
         write(lu,900) myname, i, trim(ods%meta%kx_names(i))
900      format(a,': Kx = ',I3,' ',a)
         
         do j = 1, ktmax
            
            if( ksums(i,j,2).ne.0 ) then
               write(lu,910) myname, j, trim(ods%meta%kt_names(j)), ksums(i,j,1:2)
910            format(a,':     kt = ',I3,' ',a32,2i10)
            endif
            
         end do
         
      endif
      
   end do
   end if 


!  Write the observation summary table for all data sources
!  --------------------------------------------------------
   kgrand = 0
   write(lu,'(a)') myname // ':'
   write(lu,900) myname, 9999, 'ALL KXs '

   do j = 1, ktmax

      kgrand = kgrand + ktots(j,1:2)
      if ( ktots(j,2) .ne. 0 ) then
         write(lu,910) myname, j, trim(ods%meta%kt_names(j)), ktots(j,1:2)
      end if

   end do
   write(lu,'(a)') myname // ':'
   write(lu,910) myname, 9999, 'ALL KTs ', kgrand 
   write(lu,'(a)') myname // ':'

   deallocate ( ksums, ktots )
   
   rc = 0
   return

 end subroutine ODS_Tally


!-------------------------------------------------------------------------
!         NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE:  ODS_Select - Select observations by attribute values
!
! !INTERFACE:
!
   subroutine ODS_Select ( ods, nobs, nsel, rc,                    &
                           odss,                                   & ! optional
                           complement,                             & ! optional
                              kid      ,   kx      ,     kt      , & ! optional
                               ks      , time      , qcexcl      , & ! optional
                           qchist      ,  lat      ,    lon      , & ! optional
                              lev      ,   xm      ,    obs      , & ! optional
                              OmF      ,  OmA      ,   Xvec      , & ! optional
                              kid_list ,   kx_list ,     kt_list , & ! optional
                               ks_list , time_list , qcexcl_list , & ! optional
                           qchist_list ,  lat_list ,    lon_list , & ! optional
                              lev_list ,   xm_list ,    obs_list , & ! optional
                              OmF_list ,  OmA_list ,   Xvec_list , & ! optional
                              kid_range,   kx_range,     kt_range, & ! optional
                               ks_range, time_range, qcexcl_range, & ! optional
                           qchist_range,  lat_range,    lon_range, & ! optional
                              lev_range,   xm_range,    obs_range, & ! optional
                              OmF_range,  OmA_range,   Xvec_range  ) ! optional

! !USES:
  use parutilitiesmodule, only : commglobal, gid, gsize, MAXOP,       &
                                 parcollective
  implicit none

! !INPUT PARAMETERS:

  type(ods_vect), intent(inout) :: ods        ! obs vector
  integer, intent(in)           :: nobs       ! number of obs
  logical, intent(in), optional :: complement ! apply complement of selection criteria

                                              ! selection criteria:
                                              ! ------------------

  integer, intent(in), optional :: kid          !  identification index
  real,    intent(in), optional :: lat          !  latitute
  real,    intent(in), optional :: lon          !  longitude
  real,    intent(in), optional :: lev          !  level
  integer, intent(in), optional :: kx           !  data source index
  integer, intent(in), optional :: kt           !  data type index
  integer, intent(in), optional :: ks           !  sounding index
  real,    intent(in), optional :: xm           !  metadata
  integer, intent(in), optional :: time         !  time
  real,    intent(in), optional :: obs          !  obs value
  real,    intent(in), optional :: OmF          !  obs minus forecast
  real,    intent(in), optional :: OmA          !  obs minus analysis
  real,    intent(in), optional :: Xvec         !  PSAS CG solution vector
  integer, intent(in), optional :: qcexcl       !  QC exclusion flag
  integer, intent(in), optional :: qchist       !  QC history flag

  integer, intent(in), optional :: kid_list(:)    !  list of identification indices
  real,    intent(in), optional :: lat_list(:)    !  list of latitutes
  real,    intent(in), optional :: lon_list(:)    !  list of longitudes
  real,    intent(in), optional :: lev_list(:)    !  list of levels
  integer, intent(in), optional :: kx_list(:)     !  list of data source indices
  integer, intent(in), optional :: kt_list(:)     !  list of data type indices
  integer, intent(in), optional :: ks_list(:)     !  list of sounding indices
  real,    intent(in), optional :: xm_list(:)     !  list of metadata
  integer, intent(in), optional :: time_list(:)   !  list of times
  real,    intent(in), optional :: obs_list(:)    !  list of obs values
  real,    intent(in), optional :: OmF_list(:)    !  list of obs minus forecasts
  real,    intent(in), optional :: OmA_list(:)    !  list of obs minus analyses
  real,    intent(in), optional :: Xvec_list(:)   !  list of PSAS CG sol vec
  integer, intent(in), optional :: qcexcl_list(:) !  list of QC exclusion flags
  integer, intent(in), optional :: qchist_list(:) !  list of QC history flags

  integer, intent(in), optional :: kid_range(2)     !  range of identification indices
  real,    intent(in), optional :: lat_range(2)     !  range of latitutes
  real,    intent(in), optional :: lon_range(2)     !  range of longitudes
  real,    intent(in), optional :: lev_range(2)     !  range of levels
  integer, intent(in), optional :: kx_range(2)      !  range of data source indices
  integer, intent(in), optional :: kt_range(2)      !  range of data type indices
  integer, intent(in), optional :: ks_range(2)      !  range of sounding indices
  real,    intent(in), optional :: xm_range(2)      !  range of metadata
  integer, intent(in), optional :: time_range(2)    !  range of times
  real,    intent(in), optional :: obs_range(2)     !  range of obs values
  real,    intent(in), optional :: OmF_range(2)     !  range of obs minus forecasts
  real,    intent(in), optional :: OmA_range(2)     !  range of obs minus analyses
  real,    intent(in), optional :: Xvec_range(2)    !  range of PSAS CG sol vecs
  integer, intent(in), optional :: qcexcl_range(2)  !  range of QC exclusion flags
  integer, intent(in), optional :: qchist_range(2)  !  range of QC history flags

! !OUTPUT PARAMETERS:

  type(ods_vect), intent(inout), optional :: odss     ! selected obs vector
  integer, intent(out)                    :: nsel     ! number of selected obs
  integer, intent(out)                    :: rc       ! return code:
                                                    !  0: all ok
                                                    !  1: allocate error

! !DESCRIPTION: Selects from the input observation vector ods according to
!   a given set of conditions, specified in terms of values and/or enumerated
!   lists and/or ranges of any of the atomic attributes. For example,
!\begin{verbatim}
!     call ODS_Select ( ods, nobs, nsel, rc,       &
!                       qcexcl=0, kx_list=(/14,16/), lev_range=(/150.,250./) )
!\end{verbatim}
!   selects observations that simultaneously satisfy all specified criteria.
!   Alternatively,
!\begin{verbatim}
!     call ODS_Select ( ods, nobs, nsel, rc,       &
!                       complement=.true.          &
!                       qcexcl=0, kx_list=(/14,16/), lev_range=(/150.,250./) )
!\end{verbatim}
!   selects observations that do not satisfy any of the specified conditions.
!
!   This procedure only sees the first nobs observations in the input observation
!   vector ods; nobs may be less than the actual number of observations in ods.
!
!   On return, nsel is the number of selected observations, and the first nobs
!   observations in ods are reordered such that the selected observations are
!   located at the front. Thus this procedure can be applied successively in
!   order to combine selection criteria, for example as in
!\begin{verbatim}
!     call ODS_Select ( ods, nobs, n   , rc, kt=6 )
!     call ODS_Select ( ods, n   , nsel, rc, complement=.true., qcexcl=0 )
!\end{verbatim}
!   which will select all height observations that did not pass quality
!   control.
!
!   An output obs vector can be specified for the selected observations:
!\begin{verbatim}
!     call ODS_Select ( ods, nobs, nsel, rc,       &
!                       odss,                      &
!                       qcexcl=0, kx_list=(/14,16/), lev_range=(/150.,250./) )
!\end{verbatim}
!   in which case the original ods will remain unchanged.
!
!     Notes:
!\begin{itemize}
!\item This routine must be called with keyword arguments for the
!      selection criteria
!\item Each range includes the endpoints, BUT
!\item It is up to the user to deal with possible effects of floating
!      point arithmetic on equality tests for real attributes. For
!      example, use
!\begin{verbatim}
!      lev_range=(/500.-epsilon(1.),500.+epsilon(1.)/)
!\end{verbatim}
!\item There is no check for inconsistent conditions; the result will
!      be \verb|nsel=0, rc=0|
!\item Memory for output obs vector is allocated by this routine
!\end{itemize}
!
! !REVISION HISTORY:
!
!  01Nov2000   Dee        Initial code
!  16Feb2001   Todling    Changed to account for Xvec addition
!  23Mar2001   Redder     Correct latex bugs in prologue
!  04Oct2001   G Lou      Added nsyn for RUC   
!  18Jun2002   Zaslavsky  Made changes related to relocation of nsyn from 
!                         metadata to data
!
!EOP
!-------------------------------------------------------------------------

   character(len=*), parameter :: myname = 'ODS_Select'

!  Local storage for selection:
!  ---------------------------
   logical, allocatable ::  mask(:), selected(:)

!  Etc.
!  ---
   integer ier, i
   integer nnkt, nnkx, nnqc, nncr
   integer nnsyn
   integer, allocatable :: num_in_lat(:)

   if ( nobs .eq. 0 ) then
        nsel = 0
        rc   = 0
        return
   end if

   allocate ( mask(nobs), selected(nobs), stat=ier )
   if ( ier .ne. 0 ) then
        rc = 1
        return
   end if

!  All data are selected to begin with
!  -----------------------------------
   selected = .true.

!  Select on single attribute values
!  ---------------------------------
   if (present(kid   )) selected = selected .and. ods%data%kid   (1:nobs) .eq. kid
   if (present(lat   )) selected = selected .and. ods%data%lat   (1:nobs) .eq. lat
   if (present(lon   )) selected = selected .and. ods%data%lon   (1:nobs) .eq. lon
   if (present(lev   )) selected = selected .and. ods%data%lev   (1:nobs) .eq. lev
   if (present(kx    )) selected = selected .and. ods%data%kx    (1:nobs) .eq. kx
   if (present(kt    )) selected = selected .and. ods%data%kt    (1:nobs) .eq. kt
   if (present(ks    )) selected = selected .and. ods%data%ks    (1:nobs) .eq. ks
   if (present(xm    )) selected = selected .and. ods%data%xm    (1:nobs) .eq. xm
   if (present(time  )) selected = selected .and. ods%data%time  (1:nobs) .eq. time
   if (present(obs   )) selected = selected .and. ods%data%obs   (1:nobs) .eq. obs
   if (present(OmF   )) selected = selected .and. ods%data%OmF   (1:nobs) .eq. OmF
   if (present(OmA   )) selected = selected .and. ods%data%OmA   (1:nobs) .eq. OmA
   if (present(Xvec  )) selected = selected .and. ods%data%Xvec  (1:nobs) .eq. Xvec
   if (present(qcexcl)) selected = selected .and. ods%data%qcexcl(1:nobs) .eq. qcexcl
   if (present(qchist)) selected = selected .and. ods%data%qchist(1:nobs) .eq. qchist

!  Select on ranges of attribute values
!  ------------------------------------
   if (present(kid_range   )) selected = selected .and. &
       kid_range(1)    .le. ods%data%kid   (1:nobs) .and. &
                            ods%data%kid   (1:nobs) .le. kid_range(2)
   if (present(lat_range   )) selected = selected .and. &
       lat_range(1)    .le. ods%data%lat   (1:nobs) .and. &
                            ods%data%lat   (1:nobs) .le. lat_range(2)
   if (present(lon_range   )) selected = selected .and. &
       lon_range(1)    .le. ods%data%lon   (1:nobs) .and. &
                            ods%data%lon   (1:nobs) .le. lon_range(2)
   if (present(lev_range   )) selected = selected .and. &
       lev_range(1)    .le. ods%data%lev   (1:nobs) .and. &
                            ods%data%lev   (1:nobs) .le. lev_range(2)
   if (present(kx_range    )) selected = selected .and. &
       kx_range(1)     .le. ods%data%kx    (1:nobs) .and. &
                            ods%data%kx    (1:nobs) .le. kx_range(2)
   if (present(kt_range    )) selected = selected .and. &
       kt_range(1)     .le. ods%data%kt    (1:nobs) .and. &
                            ods%data%kt    (1:nobs) .le. kt_range(2)
   if (present(ks_range    )) selected = selected .and. &
       ks_range(1)     .le. ods%data%ks    (1:nobs) .and. &
                            ods%data%ks    (1:nobs) .le. ks_range(2)
   if (present(xm_range    )) selected = selected .and. &
       xm_range(1)     .le. ods%data%xm    (1:nobs) .and. &
                            ods%data%xm    (1:nobs) .le. xm_range(2)
   if (present(time_range  )) selected = selected .and. &
       time_range(1)   .le. ods%data%time  (1:nobs) .and. &
                            ods%data%time  (1:nobs) .le. time_range(2)
   if (present(obs_range   )) selected = selected .and. &
       obs_range(1)    .le. ods%data%obs   (1:nobs) .and. &
                            ods%data%obs   (1:nobs) .le. obs_range(2)
   if (present(OmF_range   )) selected = selected .and. &
       OmF_range(1)    .le. ods%data%OmF   (1:nobs) .and. &
                            ods%data%OmF   (1:nobs) .le. OmF_range(2)
   if (present(OmA_range   )) selected = selected .and. &
       OmA_range(1)    .le. ods%data%OmA   (1:nobs) .and. &
                            ods%data%OmA   (1:nobs) .le. OmA_range(2)
   if (present(Xvec_range  )) selected = selected .and. &
       Xvec_range(1)   .le. ods%data%Xvec  (1:nobs) .and. &
                            ods%data%Xvec  (1:nobs) .le. Xvec_range(2)
   if (present(qcexcl_range)) selected = selected .and. &
       qcexcl_range(1) .le. ods%data%qcexcl(1:nobs) .and. &
                            ods%data%qcexcl(1:nobs) .le. qcexcl_range(2)
   if (present(qchist_range)) selected = selected .and. &
       qchist_range(1) .le. ods%data%qchist(1:nobs) .and. &
                            ods%data%qchist(1:nobs) .le. qchist_range(2)

!  Select on lists of attribute values
!  -----------------------------------
   if (present(kid_list   )) &
       call slist_int_ (ods%data%kid   , nobs, kid_list   , size(kid_list   ))
   if (present(lat_list   )) &
       call slist_real_(ods%data%lat   , nobs, lat_list   , size(lat_list   ))
   if (present(lon_list   )) &
       call slist_real_(ods%data%lon   , nobs, lon_list   , size(lon_list   ))
   if (present(lev_list   )) &
       call slist_real_(ods%data%lev   , nobs, lev_list   , size(lev_list   ))
   if (present(kx_list    )) &
       call slist_int_ (ods%data%kx    , nobs, kx_list    , size(kx_list    ))
   if (present(kt_list    )) &
       call slist_int_ (ods%data%kt    , nobs, kt_list    , size(kt_list    ))
   if (present(ks_list    )) &
       call slist_int_ (ods%data%ks    , nobs, ks_list    , size(ks_list    ))
   if (present(xm_list    )) &
       call slist_real_(ods%data%xm    , nobs, xm_list    , size(xm_list    ))
   if (present(time_list  )) &
       call slist_int_ (ods%data%time  , nobs, time_list  , size(time_list  ))
   if (present(obs_list   )) &
       call slist_real_(ods%data%obs   , nobs, obs_list   , size(obs_list   ))
   if (present(OmF_list   )) &
       call slist_real_(ods%data%OmF   , nobs, OmF_list   , size(OmF_list   ))
   if (present(OmA_list   )) &
       call slist_real_(ods%data%OmA   , nobs, OmA_list   , size(OmA_list   ))
   if (present(Xvec_list  )) &
       call slist_real_(ods%data%Xvec  , nobs, Xvec_list  , size(Xvec_list  ))
   if (present(qcexcl_list)) &
       call slist_int_ (ods%data%qcexcl, nobs, qcexcl_list, size(qcexcl_list))
   if (present(qchist_list)) &
       call slist_int_ (ods%data%qchist, nobs, qchist_list, size(qchist_list))

!  Apply complement
!  ----------------
   if ( present(complement) ) then
        if ( complement ) selected = .not. selected
   end if

!  Check if any data remain
!  ------------------------
   nsel = count(selected)
   if ( nsel .eq. 0 ) then
        deallocate ( mask, selected )
        rc = 0
        return
   end if

!  Organize selected obs in ods structure
!  --------------------------------------
   if ( present(odss) ) then

!     allocate memory for output ods structure
!     ----------------------------------------
      nnkt = ods%meta%nkt
      nnkx = ods%meta%nkx
      nnqc = ods%meta%nqc
      nncr = ods%meta%ncr
      nnsyn = ods%data%nsyn

      allocate( num_in_lat( gsize ) )
      num_in_lat = 0
      num_in_lat(gid+1) = nobs
      call parcollective( commglobal, MAXOP, gsize, num_in_lat )
      call ODS_Init ( odss, nsel, ier, nnkt, nnkx, nnqc, nncr, nnsyn, &
                      dist=num_in_lat )
      if (ier .ne. 0) then
          deallocate ( mask, selected )
          rc = 1
          return
      end if
      deallocate( num_in_lat )

!     copy metadata from input ODS
!     ----------------------------
      odss%meta%kt_names  = ods%meta%kt_names
      odss%meta%kt_units  = ods%meta%kt_units
      odss%meta%kx_names  = ods%meta%kx_names
      odss%meta%kx_meta   = ods%meta%kx_meta
      odss%meta%qcx_names = ods%meta%qcx_names

!     set data attributes
!     -------------------
      odss%data%kid    = (/ (i,i=1,nsel) /)
      odss%data%gid    = pack(ods%data%gid   (1:nobs), mask=selected)
      odss%data%lat    = pack(ods%data%lat   (1:nobs), mask=selected)
      odss%data%lon    = pack(ods%data%lon   (1:nobs), mask=selected)
      odss%data%lev    = pack(ods%data%lev   (1:nobs), mask=selected)
      odss%data%kx     = pack(ods%data%kx    (1:nobs), mask=selected)
      odss%data%kt     = pack(ods%data%kt    (1:nobs), mask=selected)
      odss%data%ks     = pack(ods%data%ks    (1:nobs), mask=selected)
      odss%data%xm     = pack(ods%data%xm    (1:nobs), mask=selected)
      odss%data%time   = pack(ods%data%time  (1:nobs), mask=selected)
      odss%data%obs    = pack(ods%data%obs   (1:nobs), mask=selected)
      odss%data%OmF    = pack(ods%data%OmF   (1:nobs), mask=selected)
      odss%data%OmA    = pack(ods%data%OmA   (1:nobs), mask=selected)
      odss%data%Xvec   = pack(ods%data%Xvec  (1:nobs), mask=selected)
      odss%data%qcexcl = pack(ods%data%qcexcl(1:nobs), mask=selected)
      odss%data%qchist = pack(ods%data%qchist(1:nobs), mask=selected)

   else

!     move selected obs to front of input structure
!     ---------------------------------------------
      ods%data%gid   (1:nobs) = (/ pack(ods%data%gid   (1:nobs), mask=      selected),  &
                                   pack(ods%data%gid   (1:nobs), mask=.not. selected) /)
      ods%data%kid   (1:nobs) = (/ pack(ods%data%kid   (1:nobs), mask=      selected),  &
                                   pack(ods%data%kid   (1:nobs), mask=.not. selected) /)
      ods%data%lat   (1:nobs) = (/ pack(ods%data%lat   (1:nobs), mask=      selected),  &
                                   pack(ods%data%lat   (1:nobs), mask=.not. selected) /)
      ods%data%lon   (1:nobs) = (/ pack(ods%data%lon   (1:nobs), mask=      selected),  &
                                   pack(ods%data%lon   (1:nobs), mask=.not. selected) /)
      ods%data%lev   (1:nobs) = (/ pack(ods%data%lev   (1:nobs), mask=      selected),  &
                                   pack(ods%data%lev   (1:nobs), mask=.not. selected) /)
      ods%data%kx    (1:nobs) = (/ pack(ods%data%kx    (1:nobs), mask=      selected),  &
                                   pack(ods%data%kx    (1:nobs), mask=.not. selected) /)
      ods%data%kt    (1:nobs) = (/ pack(ods%data%kt    (1:nobs), mask=      selected),  &
                                   pack(ods%data%kt    (1:nobs), mask=.not. selected) /)
      ods%data%ks    (1:nobs) = (/ pack(ods%data%ks    (1:nobs), mask=      selected),  &
                                   pack(ods%data%ks    (1:nobs), mask=.not. selected) /)
      ods%data%xm    (1:nobs) = (/ pack(ods%data%xm    (1:nobs), mask=      selected),  &
                                   pack(ods%data%xm    (1:nobs), mask=.not. selected) /)
      ods%data%time  (1:nobs) = (/ pack(ods%data%time  (1:nobs), mask=      selected),  &
                                   pack(ods%data%time  (1:nobs), mask=.not. selected) /)
      ods%data%obs   (1:nobs) = (/ pack(ods%data%obs   (1:nobs), mask=      selected),  &
                                   pack(ods%data%obs   (1:nobs), mask=.not. selected) /)
      ods%data%OmF   (1:nobs) = (/ pack(ods%data%OmF   (1:nobs), mask=      selected),  &
                                   pack(ods%data%OmF   (1:nobs), mask=.not. selected) /)
      ods%data%OmA   (1:nobs) = (/ pack(ods%data%OmA   (1:nobs), mask=      selected),  &
                                   pack(ods%data%OmA   (1:nobs), mask=.not. selected) /)
      ods%data%Xvec  (1:nobs) = (/ pack(ods%data%Xvec  (1:nobs), mask=      selected),  &
                                   pack(ods%data%Xvec  (1:nobs), mask=.not. selected) /)
      ods%data%qcexcl(1:nobs) = (/ pack(ods%data%qcexcl(1:nobs), mask=      selected),  &
                                   pack(ods%data%qcexcl(1:nobs), mask=.not. selected) /)
      ods%data%qchist(1:nobs) = (/ pack(ods%data%qchist(1:nobs), mask=      selected),  &
                                   pack(ods%data%qchist(1:nobs), mask=.not. selected) /)

   end if

   deallocate ( mask, selected )

!  All done
!  --------
   rc = 0
   return


CONTAINS   ! ................. internals follow ......................

     subroutine slist_int_ ( attr, nobs, list, nlst )

!    select on list of integer attribute values

      implicit none
      integer, intent(in)     ::      nobs
      integer, intent(in)     ::      nlst
      integer, intent(in)     :: attr(nobs)
      integer, intent(in)     :: list(nlst)

      mask = .false.
      do i = 1,nlst
         where ( attr .eq. list(i) )
             mask = .true.
         end where
      end do
      selected = selected .and. mask

     end subroutine slist_int_

     subroutine slist_real_ ( attr, nobs, list, nlst )

!    select on list of real attribute values

      implicit none
      integer, intent(in)     ::      nobs
      integer, intent(in)     ::      nlst
      real,    intent(in)     :: attr(nobs)
      real,    intent(in)     :: list(nlst)

      mask = .false.
      do i = 1,nlst
         where ( attr .eq. list(i) )
             mask = .true.
         end where
      end do
      selected = selected .and. mask

     end subroutine slist_real_

   end subroutine ODS_Select

!-------------------------------------------------------------------------
!         NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE:  ODS_Transfer --- Transfers (redistributes) ODS vector
!
! !INTERFACE:
!
      subroutine ODS_Transfer ( transfer, ods_in, ods_out, rc )

! !USES:
      use parutilitiesmodule, only : commglobal, BCSTOP, parcollective,    &
                                     gid, parbegintransfer, parendtransfer
      use m_odstransfer, only      : transfertype, transfer_perform
      implicit none

! !INPUT/OUTPUT PARAMETERS:

  type(transfertype), intent(in) ::  transfer      ! Transfer pattern
  type(ods_vect), intent(in)     ::  ods_in        ! Many ODS vectors
  type(ods_vect), intent(inout)  ::  ods_out       ! Merge ODS vector
  integer,        intent(out)    ::  rc            ! Error return code:
                                                   !  0 - all is well
                                                   !  1 - could not allocate mem

! !DESCRIPTION:
! \label{MODS:Merge}
!    Transfer one ods vector to another.  This is fundamentally 
!    designed to support SPMD mode, and makes use of the PILGRIM
!    communication primitives.   Do not use this to copy one
!    vector to a new one in sequential mode -- {\tt OBS\_Merge} is much
!    better for this purpose.
!
! !REVISION HISTORY:
!
!  21mar2002   Sawyer    Initial code
!  25Jun2002   Zaslavsky Made changes related to relocation of nsyn from 
!                        metadata to data
!
!EOP
!-------------------------------------------------------------------------
  integer :: iv(5), i, j, ier
  character(len=NCR_MAX), allocatable :: buffer(:)

  rc = 1

! Replicate Metadata  -- assume correct copy on gid==0
! ----------------------------------------------------

  if ( gid == 0 ) then
    iv(1) = ods_in%meta%nkt
    iv(2) = ods_in%meta%nkx
    iv(3) = ods_in%meta%nqc
    iv(4) = ods_in%meta%ncr
    iv(5) = ods_in%data%nsyn  ! nsyn has been recently transferred from 
                              ! metadata to data
  end if
  call parcollective( commglobal, BCSTOP, 5, iv )

  ods_out%meta%nkt = iv(1)
  ods_out%meta%nkx = iv(2)
  ods_out%meta%nqc = iv(3) 
  ods_out%meta%ncr = iv(4)
  ods_out%data%nsyn= iv(5)    ! nsyn has been recently transferred from 
                              ! metadata to data

  allocate( buffer( 2*iv(1)+2*iv(2)+iv(3) ) )
  if (gid == 0 ) then
    buffer(1:iv(1)) = ods_in%meta%kt_names(1:iv(1))
    buffer(iv(1)+1:2*iv(1)) = ods_in%meta%kt_units(1:iv(1))
    buffer(2*iv(1)+1:2*iv(1)+iv(2)) = ods_in%meta%kx_names(1:iv(2))
    buffer(2*iv(1)+iv(2)+1:2*iv(1)+2*iv(2)) = ods_in%meta%kx_meta(1:iv(2))
    buffer(2*iv(1)+2*iv(2)+1:2*iv(1)+2*iv(2)+iv(3)) = ods_in%meta%qcx_names(1:iv(3))
  end if

  call parcollective( commglobal, BCSTOP, 2*iv(1)+2*iv(2)+iv(3), buffer )

!
! Now the meta strings
! --------------------

  ods_out%meta%kt_names(1:iv(1)) = buffer(1:iv(1))
  ods_out%meta%kt_units(1:iv(1)) = buffer(iv(1)+1:2*iv(1))
  ods_out%meta%kx_names(1:iv(2)) = buffer(2*iv(1)+1:2*iv(1)+iv(2))
  ods_out%meta%kx_meta(1:iv(2))  = buffer(2*iv(1)+iv(2)+1:2*iv(1)+2*iv(2))
  ods_out%meta%qcx_names(1:iv(3))= buffer(2*iv(1)+2*iv(2)+1:2*iv(1)+2*iv(2)+iv(3))
  deallocate( buffer )

  call transfer_perform( transfer, ods_in%data, ods_out%data, rc )

  return

 end subroutine ODS_Transfer


!-------------------------------------------------------------------------
!         NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE:  ODS_Partition_ --- Partition and Scatter ODS
!
! !INTERFACE:
  subroutine ODS_Partition_( ods_g, ods, rc, lat_max )

! !USES:
    use parutilitiesmodule, only : gid, gsize, commglobal, BCSTOP,   &
                                   parcollective
    use m_odstransfer, only : TransferType, Transfer_Init, Transfer_Clean
    use m_MPodsdata, only : odsd_sort

! !INPUT/OUTPUT PARAMETERS:
    type(ods_vect), intent(inout)    ::  ods_g      ! Global ODS vector (PE 0)
    type(ods_vect), intent(inout)    ::  ods        ! ODS vector (partitioned)
    real, intent(inout), optional    ::  lat_max(:) ! latitude maxima

! !OUTPUT PARAMETERS:
    integer, intent(out)             ::  rc         ! Error code

! !DESCRIPTION:
! \label{MODS:Partition_}
!  This routine determines the partition of the ods vector
!  which is globally defined on PE 0.  It first sorts the 
!  vector by latitude.  If {\tt lat\_max} is defined it accepts 
!  this as the definition of the latitude bands; if not, it
!  determines an optimal distribution of the observations
!  which evenly distributes the number of obs.  It then sets
!  the meta data in the output ods, determines the communication 
!  pattern for the scattering operation and scatters the data
!  to the output {\tt ods}.
!
!  This routine is internal to {\tt m\_ods} (not publically 
!  accessible).
!
! !REVISION HISTORY:
!
!  12Feb2003  Sawyer    Created portions of ods_get1_
!
!EOP
!-------------------------------------------------------------------------

! !LOCAL VARIABLES:
   type (transfertype)  :: scatter ! transfer global -> local data sets

   real    :: lat_tmp
   real    :: lat_range(2)
   integer :: iv(5)
   integer :: nobs_inh, nobs_gl, i, ipe, oldnsel, nsel
   integer :: cnt, nsup, lat_weight_sum, cnt_est
   real, allocatable    :: lat_sup(:), lat_weight(:)
   integer, allocatable :: num_in_lat(:), gl_dist(:)
   real, allocatable    :: latmax(:)

! Start of code
     allocate( latmax(gsize) )
     allocate( num_in_lat(gsize) )
     if ( gid == 0 ) then
        nobs_gl = ods_g%data%nobs

!
! Sort by latitude;  this is useful for the subsequent decomposition
! ------------------------------------------------------------------
        call odsd_sort( ods_g%data, nobs_gl, rc, attr1='ks',        &
                        attr2='kt', attr3='lat' )

!
! Reinitialize the global IDs for this ordering
! ---------------------------------------------
        ods_g%data%gid = (/ (i, i=1,nobs_gl) /)

!
! If lat_max is predefined
! ------------------------
        if ( present( lat_max ) .and. lat_max(gsize) == 90.0 ) then
           latmax = lat_max
!
! Distribute the observations by the specified latitudes
! ------------------------------------------------------
           oldnsel = 0
           lat_range(1) = -90.0
           do ipe=1,gsize
              lat_range(2) = lat_max(ipe)
              call ODS_select( ods_g, nobs_gl, nsel, rc, lat_range=lat_range )
              num_in_lat(ipe) = nsel - oldnsel
              oldnsel = nsel
           end do
        else
!
! Distribute approximately equal numbers of observations to the PEs
! -----------------------------------------------------------------
           allocate( lat_sup(nobs_gl) )
           allocate( lat_weight(nobs_gl) )

           nsup = 1
           lat_sup(1) = ods_g%data%lat(1)
           lat_tmp    = ods_g%data%lat(1)    ! added

           lat_weight(1) = 1
           do cnt = 2, nobs_gl
              if ( lat_tmp > ods_g%data%lat(cnt) ) then
                 print *, "ods_partition: TROUBLE lat(",cnt,") = ", &
                          ods_g%data%lat(cnt), " < previous", lat_tmp
              endif
              lat_tmp = ods_g%data%lat(cnt)
              if( lat_tmp .gt. lat_sup(nsup) ) then 
                 nsup = nsup + 1
                 lat_sup(nsup) = lat_tmp
                 lat_weight(nsup) = 1
              else
                 lat_weight(nsup) = lat_weight(nsup) + 1
              endif
           enddo

           oldnsel = 0
           cnt = 0
           lat_weight_sum = 0

           do ipe = 1, gsize-1
              cnt_est = (ipe*nobs_gl)/gsize    ! estimate for decomposition

              if ( cnt .lt. nsup) then

                 cnt = cnt + 1
                 lat_weight_sum = lat_weight_sum + lat_weight(cnt)
                   
                 do while ( lat_weight_sum+lat_weight (cnt+1) .le. cnt_est )
                    cnt = cnt + 1
                    lat_weight_sum = lat_weight_sum + lat_weight(cnt)
                 enddo
              endif

              num_in_lat(ipe) = lat_weight_sum - oldnsel
              oldnsel = lat_weight_sum

              latmax(ipe) = (lat_sup(cnt+1) + lat_sup(cnt)) / 2.0
           enddo

           num_in_lat(gsize) = nobs_gl - oldnsel

           latmax(gsize) = 90.0
           deallocate(lat_sup)
           deallocate(lat_weight)
        end if
     endif
     call parcollective( commglobal, BCSTOP, nobs_gl )
     call parcollective( commglobal, BCSTOP, gsize, num_in_lat )
     call parcollective( commglobal, BCSTOP, gsize, latmax )

     if ( gid == 0 ) then
        iv(1) = ods_g%meta%nkt
        iv(2) = ods_g%meta%nkx
        iv(3) = ods_g%meta%nqc
        iv(4) = ods_g%meta%ncr
        iv(5) = ods_g%data%nsyn
     end if

     call parcollective( commglobal, BCSTOP, 5, iv )

!
! Allocate a minimal ods vector on all other PEs
! ----------------------------------------------
     allocate( gl_dist(gsize) )
     gl_dist = 0
     gl_dist(1) = nobs_gl
     if (gid /= 0) then
        nsel = 0
        call ods_init( ods_g, nsel, rc, &
                       nkt=iv(1), nkx=iv(2), nqc=iv(3), ncr=iv(4),      &
                       nsyn=iv(5), dist=gl_dist )
     endif
     deallocate( gl_dist )


     nobs_inh = num_in_lat( gid+1 )    !  Number of obs on this PE
  
!
! Initialize the local ods vector  
! (User should ensure that ods has been cleaned already)
! ------------------------------------------------------

     lat_range(1) = -90.0
     if ( gid > 0 ) lat_range(1) = latmax(gid)
     lat_range(2) = latmax(gid+1)

     call ods_init( ods, nobs_inh, rc, &
                    nkt=iv(1), nkx=iv(2), nqc=iv(3), ncr=iv(4),    &
                    nsyn=iv(5), dist=num_in_lat,                   &
                    lmin=lat_range(1), lmax=lat_range(2) )
     if ( rc .ne. 0 ) then
        print *, 'Error initializing ods ', 'rc = ', rc
     end if

!
! Transfer the global ods vector to the local version
! ---------------------------------------------------
     call Transfer_Init( ods_g%decomp, ods%decomp, scatter, rc )
     call ODS_Transfer( scatter, ods_g, ods, rc )
     call Transfer_Clean(  scatter, rc )

     call ods_clean( ods_g, rc )


!
! Test the latitudes to make sure they are on the correct PE
! ----------------------------------------------------------
 
    cnt = 0
    do i=1,nobs_inh
       if ( ods%data%lat(i) > latmax(gid+1) .and. cnt > 100 ) then
          cnt = cnt + 1
          print *, "trouble on pe", gid, "index", i, "lat", ods%data%lat(i), &
               "not less than", latmax(gid+1)
       endif
    end do

    deallocate( num_in_lat )
!
!    Return lat_max if it is present but not defined
!    -----------------------------------------------
    if ( present( lat_max ) .and. lat_max(gsize) /= 90.0 ) then
       lat_max = latmax
    endif
    deallocate( latmax )


! All done
! --------
  return

end subroutine ODS_Partition_


!-------------------------------------------------------------------------
!         NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE:  ODS_get1pe --- Reads data from a single ODS file
! 
! !INTERFACE:

  subroutine ods_get1pe ( fname, nymd, nhms, ftype, ods, rc, &
                           ncid )                    ! Optional

!
! !USES
!
    use parutilitiesmodule, only : gsize, gid
    Implicit NONE

! !INPUT PARAMETERS: 
!
    character(len=*), intent(in)   :: fname   ! ODS file name
    integer, intent(in)            :: nymd    ! year-month-day, e.g., 19990701
    integer, intent(in)            :: nhms    ! hour-min-sec,   e.g., 120000

    integer, intent(in), optional  :: ncid    ! ODS file id as returned by
                                              ! ODS_Open(); in this case
                                              ! fname is not used
!
! !OUTPUT PARAMETERS:
!
    character(len=*), intent(out)  :: ftype   ! ODS file type: 'pre_anal'
                                              !  or 'post_anal'
    type(ods_vect), intent(inout)  :: ods     ! ODS vector

    integer, intent(out)           :: rc      ! Error return code:
                                              !  0 - all is well
                                              !  1 - could not open file
                                              !  2 - could not get table sizes
                                              !  3 - could not get no. obs
                                              !  4 - could not allocate memory
                                              !  5 - could not read tables
                                              !  6 - could not read attributes (pre)
                                              !  7 - could not read attributes (post)
                                              !  8 - could not obtain first julian day on file

! !DESCRIPTION:
! \label{MODS:Get1pe_}
!  This routine reads data from an ODS file, allocates the necessary
!  memory for the ODS vector, and loads the data for the synoptic time
! (nymd,nhms). Usually the ODS file is opened, the data is read in and
!  the file is closed. However, if the optional parameter {\tt ncid}
!  is specified, the file is assumed to have been externally opened
!  with routine {\tt ODS\_Open()}.
!
! !REVISION HISTORY: 
!
!  07Sep1999  da Silva   Initial code.
!  09Mar2000  da Silva   Now it checks if file exists before calling
!                        ODS_Open(). 
!  12Dec2000  C. Redder  Removed the optional input parameter, whms.
!                        Modified error handling procedures to clean
!                        up (that is close the file if ncid is not
!                        present) if an error occurs.  Added code to
!                        translate output sampling time (in minutes
!                        since the current synoptic date and time)
!                        from ODS sampling time (in minutes since
!                        0:00 GMT of the first day on file).  Added
!                        to extract and save the number of synoptic
!                        time per day.
!  16Feb2001   Todling   Changed to account for Xvec addition
!  19Feb2001   Todling   Bug fix; GetList(qcx_n) arg nkx instead of nqc
!  27Feb2001   Todling   Changed close to pass blank string.
!  28Mar2001   Todling   Bug fix for nobs=0 case(it now returns)
!  04Oct2001   RT/GLou   Bug fix for nobs=0; Added nsyn for RUC
!  11Jun2002   Sawyer    Procedure ODS_Get1_ was renamed to Ods_get1pe
!  14Jun2002   RT/Redder Revised exit status for nobs=0 case
!  18Jun2002   Zaslavsky Made changes related to relocation of nsyn from 
!                        metadata to data 
!  26Oct2002   Sawyer    Now sorted by latitude
!
!EOP
!-------------------------------------------------------------------------

  character(len=*), parameter :: myname = 'ods_get1pe'
  integer id, ier, nkt, nkt1, nkx, nqc, ncr, nsyn, nobs, khms
  integer jday, syn_hour, ods_nget
  integer first_jday, first_nymd, first_nhms
  integer, allocatable :: num_in_lat(:)
  logical fexists

  rc = 0

! Determine file ID
! -----------------
  if ( present(ncid) ) then
       id = ncid
  else
       inquire ( file=trim(fname), exist=fexists )
       if ( fexists ) then
          call ODS_Open ( id, fname, 'r', ier ) ! open the file
          if ( ier .ne. 0 ) then
             rc = 1
             return
          end if

       else
          rc = 1
          return
       end if
  end if

! Determine size of ODS tables
! ----------------------------
  call ODS_IGet ( id, 'nkt',    nkt,  ier ); if ( ier .ne. 0 ) rc = 2
  call ODS_IGet ( id, 'nkx',    nkx,  ier ); if ( ier .ne. 0 ) rc = 2
  call ODS_IGet ( id, 'nqcx',   nqc,  ier ); if ( ier .ne. 0 ) rc = 2
  call ODS_IGet ( id, 'strlen', ncr,  ier ); if ( ier .ne. 0 ) rc = 2
  call ODS_IGet ( id, 'nsyn',   nsyn, ier ); if ( ier .ne. 0 ) rc = 2
  ods%data%nsyn = nsyn
  if ( rc .ne. 0 ) then
     if ( .not. present(ncid) ) call ODS_Close ( id, ' ', ier )
     return
  end if

! For mem allocation, force nkt to be at least 7 (for recasting)
! --------------------------------------------------------------
  nkt1 = max ( nkt, 7 )

! Determine how many obs in this synoptic time
! --------------------------------------------
  jday = ODS_Julian ( nymd )
  syn_hour = nhms / 10000
  nobs = ODS_NGet ( id, jday, syn_hour, ier )
  if ( ier .ne. 0 ) then
     rc = 3
     if ( .not. present(ncid) ) call ODS_Close ( id, ' ', ier )
     return
  end if

  allocate( num_in_lat( gsize ) )
  num_in_lat = 0
  num_in_lat(gid+1) = nobs

! No observations, nothing left to do
! -----------------------------------
  if ( nobs == 0 ) then
     if ( .not. present(ncid) ) call ODS_Close ( id, ' ', ier )
     call ODS_Clean ( ods, ier ) 
     call ODS_Init  ( ods, nobs, ier, dist=num_in_lat )
     if ( ier .ne. 0 ) rc = 4
     ods%data%nymd = nymd
     ods%data%nhms = nhms
     return
  end if

! Allocate memory for ODS file
! ----------------------------
  call ODS_Clean ( ods, ier )      ! ignore ier
  call ODS_Init  ( ods, nobs, ier, &
                   nkt=nkt1, nkx=nkx, nqc=nqc, ncr=ncr, nsyn=nsyn, &
                   dist=num_in_lat )

  deallocate( num_in_lat )

  ods%data%nymd = nymd
  ods%data%nhms = nhms

  if ( ier .ne. 0 ) then
     rc = 4
     if ( .not. present(ncid) ) call ODS_Close ( id, ' ', ier )
     return
  end if


! Read in ODS file tables
! -----------------------
  call ODS_GetList ( id, 'kt_names',  nkt, ods%meta%kt_names,  ier )
  if ( ier .ne. 0 ) rc = 5
  call ODS_GetList ( id, 'kt_units',  nkt, ods%meta%kt_units,  ier )
  if ( ier .ne. 0 ) rc = 5
  call ODS_GetList ( id, 'kx_names',  nkx, ods%meta%kx_names,  ier )
  if ( ier .ne. 0 ) rc = 5
  call ODS_GetList ( id, 'kx_meta',   nkx, ods%meta%kx_meta,   ier )
  if ( ier .ne. 0 ) rc = 5
  call ODS_GetList ( id, 'qcx_names', nqc, ods%meta%qcx_names, ier )
  if ( ier .ne. 0 ) rc = 5
  if ( rc .ne. 0 ) then
     if ( .not. present(ncid) ) call ODS_Close ( id, ' ', ier )
     return
  end if

! Read in pre-analysis attributes
! -------------------------------
  call ODS_GetR ( id, 'lat',    jday, syn_hour, nobs, ods%data%lat,    ier )
  if ( ier .ne. 0 ) rc = 6
  call ODS_GetR ( id, 'lon',    jday, syn_hour, nobs, ods%data%lon,    ier )
  if ( ier .ne. 0 ) rc = 6
  call ODS_GetR ( id, 'lev',    jday, syn_hour, nobs, ods%data%lev,    ier )
  if ( ier .ne. 0 ) rc = 6
  call ODS_GetI ( id, 'kx',     jday, syn_hour, nobs, ods%data%kx,     ier )
  if ( ier .ne. 0 ) rc = 6
  call ODS_GetI ( id, 'kt',     jday, syn_hour, nobs, ods%data%kt,     ier )
  if ( ier .ne. 0 ) rc = 6
  call ODS_GetI ( id, 'ks',     jday, syn_hour, nobs, ods%data%ks,     ier )
  if ( ier .ne. 0 ) rc = 6
  call ODS_GetR ( id, 'xm',     jday, syn_hour, nobs, ods%data%xm,     ier )
  if ( ier .ne. 0 ) rc = 6
  call ODS_GetI ( id, 'time',   jday, syn_hour, nobs, ods%data%time,   ier )
  if ( ier .ne. 0 ) rc = 6
  call ODS_GetR ( id, 'obs',    jday, syn_hour, nobs, ods%data%obs,    ier )
  if ( ier .ne. 0 ) rc = 6
  call ODS_GetI ( id, 'qcexcl', jday, syn_hour, nobs, ods%data%qcexcl, ier )
  if ( ier .ne. 0 ) rc = 6
  call ODS_GetI ( id, 'qchist', jday, syn_hour, nobs, ods%data%qchist, ier )
  if ( ier .ne. 0 ) rc = 6
  if ( rc .ne. 0 ) then
     if ( .not. present(ncid) ) call ODS_Close ( id, ' ', ier )
     return
  end if


! Read in post-analysis attributes if this is the case
! ----------------------------------------------------
  call ODS_CGet ( id, ':type', ftype, ier )
  if ( ier .ne. 0 ) then
       ftype = 'pre_analysis'   ! why not?
  end if
  if ( ftype(1:4) .eq. 'post' ) then
    call ODS_GetR ( id, 'omf',  jday, syn_hour, nobs, ods%data%OmF,  ier )
    if ( ier .ne. 0 ) rc = 7
    call ODS_GetR ( id, 'oma',  jday, syn_hour, nobs, ods%data%OmA,  ier )
    if ( ier .ne. 0 ) rc = 7
    call ODS_GetR ( id, 'xvec', jday, syn_hour, nobs, ods%data%Xvec, ier )
    if ( ier .ne. 0 ) ods%data%Xvec(1:nobs) = 0.0  ! for backward compatibility
  else
    ods%data%OmF (1:nobs) = 0.0  
    ods%data%OmA (1:nobs) = 0.0
    ods%data%Xvec(1:nobs) = 0.0
  end if

! Output time is minutes since current synoptic date and hour
! -----------------------------------------------------------
  call ODS_IGet ( id, 'syn_beg:first_julian_day', first_jday, ier )
  if ( ier .ne. 0 ) then
     rc = 8
     if ( .not. present(ncid) ) call ODS_Close ( id, ' ', ier )
     return

  end if
  first_nymd = ODS_CalDat ( first_jday )
  first_nhms = 0
  call ODS_NewRefTime ( first_nymd, first_nhms,   &
                              nymd,       nhms, ods%data%time )

! Close the file if we have opened it
! -----------------------------------
  if ( .not. present(ncid) ) call ODS_Close ( id, ' ', ier )
  if ( rc .ne. 0 ) return

! All done
! --------
  rc = 0
  return

  end subroutine ods_get1pe

!-------------------------------------------------------------------------
!         NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE:  ODS_Put1pe --- Writes to an ODS file from 1 process
! 
! !INTERFACE:

  subroutine ods_put1pe ( fname, ftype,  nymd, nhms, ods, rc, &
                           ncid,  append, new )         ! Optional
!
! !USES
!
    Implicit NONE

! !INPUT PARAMETERS: 
!
    character(len=*), intent(in)   :: fname   ! ODS file name
    character(len=*), intent(in)   :: ftype   ! ODS file type: 'pre_anal'
                                              !  or 'post_anal'
    integer, intent(in)            :: nymd    ! year-month-day, e.g., 19990701
    integer, intent(in)            :: nhms    ! hour-min-sec,   e.g., 120000

    integer, intent(in), optional  :: ncid    ! ODS file id as returned by
                                              ! ODS_Open(); in this case
                                              ! fname is not used

    type(ods_vect), intent(in)     :: ods     ! ODS vector

    logical, intent(in), optional  :: append  ! If true, the data is appended
                                              ! to the ODS file at this time.
                                              ! Default: append = .false.

    logical, intent(in), optional  :: new     ! If true, file will be created
                                              ! whether it exists or not
                                              ! Default: new    = .false.
!
! !OUTPUT PARAMETERS:
!
    integer, intent(out)           :: rc      ! Error return code:
                                              !  0 - all is well
                                              !  1 - could not open file
                                              !  2 - could not create file
                                              !  3 - could not append file
                                              !  4 - could not write a
                                              !     'pre_anal' observation
                                              !      attribute
                                              !  5 - could not write a
                                              !     'post_anal' observation
                                              !      attribute
                                              !  6 - error in allocating
                                              !      scratch space.
                                              !  7 - could not obtain the
                                              !      number of observations
                                              !      already on file
                                              !  8 - could not obtain the
                                              !      the first Julian day
                                              !      on file.
!
! !DESCRIPTION:
! \label{MODS:Put1pe_}
!  This routine writes an ODS vector to an ODS file at synoptic time
! {\tt (nymd,nhms)}. Usually the ODS file is opened, the data is written
!  to and the file is closed. However, if the optional parameter
! {\tt ncid} is specified, the file is assumed to have been externally
!  opened with routine {\tt ODS\_Open()}.
!
! !REVISION HISTORY: 
!
!  07Sep1999  da Silva  Initial code.
!  15Nov1999  da Silva  Added optinal parameter "new".
!  13Dec2000  C. Redder Revisions in prologue, Modified error
!                       handling procedures to clean up (that is
!                       close the file if ncid is not present)
!                       if an error occurs.   Added code to
!                       translate input sampling time (in minutes
!                       since the current synoptic date and time)
!                       to the ODS sampling time (in minutes since
!                       0:00 GMT of the first day on file).  Added
!                       code to call ODS_Append only if there are
!                       already observations stored on file for
!                       the given synoptic day and time.  Added 
!                       code to set the number of synoptic times
!                       based on the value stored in the record
!                       ods%meta in the component nsyn.
!  16Feb2001   Todling  Changed to account for Xvec addition
!  27Feb2001   Todling  Changed close to pass blank string.
!  04Oct2001   Todling  Checking return code from Create
!  07Jun2002   Sawyer   Renamed from ODS_Put to Ods_put1pe
!  18Jun2002 Zaslavsky  Made changes related to relocation of nsyn from 
!                       metadata to data 
!EOP
!-------------------------------------------------------------------------

  character(len=*), parameter :: myname = 'ods_put1pe'

  integer :: ier, ier2, jday, nobs, syn_hour, nsyn, id
  integer :: first_jday, first_nymd, first_nhms
  integer :: ods_nget, NObs_Init
  logical :: appending, fexists, creating
  integer, dimension (:), allocatable :: time  ! scratch space

! Handle optional parameters
! --------------------------
  id = -1
  if ( present(append) ) then
     appending = append
  else
     appending = .false.
  end if
  if ( present(ncid) ) then
     id = ncid
  end if
  if ( present(new) ) then
       creating = new
  else
       creating = .false.
  end if

! Check whether file exists
! -------------------------
  inquire ( file=trim(fname), exist=fexists )

! Open ...
! --------
  jday = ods_julian ( nymd )
  syn_hour = nhms / 10000
  nobs = ods%data%nobs
  if ( (id .lt. 0 .and. fexists) .AND. (.not. creating) ) then
       call ODS_Open ( id, fname, 'w', ier ) ! open the file
       if ( ier .ne. 0 ) then
          rc = 1
          return
       end if

!      ... and determine the first date and time on file
!      -------------------------------------------------
       call ODS_IGet ( id, 'syn_beg:first_julian_day', first_jday, ier )
       if ( ier .ne. 0 ) then
          rc = 8
          if ( .not. present(ncid) ) call ODS_Close ( id, ' ', ier )
          return

       end if
       first_nymd = ODS_CalDat ( first_jday )
       first_nhms = 0

! or Create
! ---------
  else if ( id .lt. 0 .or. creating ) then
       nsyn = ods%data%nsyn
       call ODS_SetParmI      ( 'nsyn', nsyn, ier )
       if ( ier .ne. 0 ) then
          rc = 2
          return
       end if
       call ODS_Create ( id, fname, ftype, jday,                             &
                         ods%meta%nkt, ods%meta%kt_names, ods%meta%kt_units, &
                         ods%meta%nkx, ods%meta%kx_names, ods%meta%kx_meta,  &
                         ods%meta%nqc, ods%meta%qcx_names,                   &
                         ier )
       if ( ier /= 0 ) then
          rc = 8
          return
       end if
       first_nymd = nymd
       first_nhms = 0

   end if

!  Determine the number of ob reports already on file
!  --------------------------------------------------
   NObs_Init = ODS_NGet ( id, jday, syn_hour, ier )
   if ( ier .ne. 0 ) then
      rc = 7
      if ( .not. present(ncid) ) call ODS_Close ( id, ' ', ier )
      return
   end if

!   If appending to this synoptic time, do the
!   setup (if ob reports are already on file)
!   ------------------------------------------
    if ( appending .and. NObs_Init .gt. 0 ) then

       call ODS_Append ( id, nobs, ier )
       if ( ier .ne. 0 ) then
          rc = 3
          if ( .not. present(ncid) ) call ODS_Close ( id, ' ', ier )
          return
       end if
     end if

!   Reference for sampling time is 0:00 GMT on the first date on file
!   -----------------------------------------------------------------
    allocate ( time ( nobs ), stat = ier )
    if ( ier .ne. 0 ) then
       rc = 6
       if ( .not. present(ncid) ) call ODS_Close ( id, ' ', ier )
       return
    end if
    time = ods%data%time
    call ODS_NewRefTime (       nymd,       nhms,   &
                          first_nymd, first_nhms, time )

!   OK, now we are ready for the heavy duty writing
!   -----------------------------------------------
  call ODS_PutR ( id, 'lat',    jday, syn_hour, nobs, ods%data%lat,    ier )
  if ( ier .ne. 0 ) rc = 4
  call ODS_PutR ( id, 'lon',    jday, syn_hour, nobs, ods%data%lon,    ier )
  if ( ier .ne. 0 ) rc = 4
  call ODS_PutR ( id, 'lev',    jday, syn_hour, nobs, ods%data%lev,    ier )
  if ( ier .ne. 0 ) rc = 4
  call ODS_PutI ( id, 'kx',     jday, syn_hour, nobs, ods%data%kx,     ier )
  if ( ier .ne. 0 ) rc = 4
  call ODS_PutI ( id, 'kt',     jday, syn_hour, nobs, ods%data%kt,     ier )
  if ( ier .ne. 0 ) rc = 4
  call ODS_PutI ( id, 'ks',     jday, syn_hour, nobs, ods%data%ks,     ier )
  if ( ier .ne. 0 ) rc = 4
  call ODS_PutR ( id, 'xm',     jday, syn_hour, nobs, ods%data%xm,     ier )
  if ( ier .ne. 0 ) rc = 4
  call ODS_PutI ( id, 'time',   jday, syn_hour, nobs, time,            ier )
  if ( ier .ne. 0 ) rc = 4
  call ODS_PutR ( id, 'obs',    jday, syn_hour, nobs, ods%data%obs,    ier )
  if ( ier .ne. 0 ) rc = 4
  call ODS_PutI ( id, 'qcexcl', jday, syn_hour, nobs, ods%data%qcexcl, ier )
  if ( ier .ne. 0 ) rc = 4
  call ODS_PutI ( id, 'qchist', jday, syn_hour, nobs, ods%data%qchist, ier )
  if ( ier .ne. 0 ) rc = 4
  deallocate ( time )
  if ( rc .ne. 0 ) then
     if ( .not. present(ncid) ) call ODS_Close ( id, ' ', ier )
     return
  end if

  if ( ftype(1:4) .eq. 'post' ) then
    call ODS_PutR ( id, 'omf',  jday, syn_hour, nobs, ods%data%OmF,  ier )
    if ( ier .ne. 0 ) rc = 5
    call ODS_PutR ( id, 'oma',  jday, syn_hour, nobs, ods%data%OmA,  ier )
    if ( ier .ne. 0 ) rc = 5
    call ODS_PutR ( id, 'xvec', jday, syn_hour, nobs, ods%data%Xvec, ier )
    if ( ier .ne. 0 ) rc = 5
  end if
 
! Close the file if user is not doing the file management
! -------------------------------------------------------
  if ( .not. present(ncid) ) call ODS_Close ( id, ' ', ier )
  if ( rc .ne. 0 ) return

! All done
! --------
  end subroutine ods_put1pe

!-------------------------------------------------------------------------
!         NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE:  ROMS_Get1pe_ --- Reads data from a iRET OMS file on 1 PE
! 
! !INTERFACE:

  subroutine ROMS_Get1pe_ ( fname, nymd, nhms, ftype, ods, rc, &
                            only_Ts )

!
! !USES
!
    use parutilitiesmodule, only : gsize, gid
    use m_MPodsmeta
    implicit NONE

! !INPUT PARAMETERS: 
!
    character(len=*), intent(in)   :: fname   ! ODS file name
    integer, intent(in)            :: nymd    ! year-month-day, e.g., 19990701
    integer, intent(in)            :: nhms    ! hour-min-sec,   e.g., 120000

    logical, intent(in), OPTIONAL  :: only_Ts ! only skin temperature

!
! !OUTPUT PARAMETERS:
!
    character(len=*), intent(out)  :: ftype   ! ODS file type: 'pre_anal'
                                              !  or 'post_anal'
    type(ods_vect), intent(out)    :: ods     ! ODS vector

    integer, intent(out)           :: rc      ! Error return code:
                                              !  0 - all is well

! !DESCRIPTION:
! \label{ROMS:Get}
!  This routine reads data from an iRET OMS file, allocates the necessary
!  memory for the ODS vector, and loads the data for the synoptic time
! (nymd,nhms). Usually the OMS file is opened, the data is read in and
!  the file is closed. 
!
! !REVISION HISTORY: 
!
!  07Sep1999  da Silva   Initial code.
!  07Feb2003  Sawyer     Moved here from m_roms, renamed roms_get1pe_
!
!EOP
!-------------------------------------------------------------------------

  character(len=*), parameter :: myname = 'roms_get1pe_'
  integer id, ier, nsyn, nobs, ios, nqc
  logical fexists

  integer :: start(2), kount(2) ! hyperslab indicators
  integer, external :: NCOPN
  integer, allocatable :: qcl2(:), land(:)
  real,    allocatable :: emissmw(:)
  integer, allocatable :: num_in_lat(:)

  rc = 0

! Only T_s is supported for now
! -----------------------------
  if ( .not. present ( only_ts ) ) then
        rc = 1
        return
  end if

! Determine file ID
! -----------------
  inquire ( file=trim(fname), exist=fexists )
  if ( fexists ) then
     id = NCOPN ( fname, 0, ier ) ! open the file
     if ( ier .ne. 0 ) then
        rc = 2
        return
     end if
  else
     rc = 2
     return
  end if

! Determine hyperslab for this time interval
! ------------------------------------------
  nsyn = 4 ! 6 hourly data for now
  call hyperslab_ ( nymd, nhms, nsyn, start, kount, rc )
  if ( rc .ne. 0 ) return

! Determine how many obs in this synoptic time
! --------------------------------------------
  nobs = kount(1)
  nqc = nqcXmax

  allocate( num_in_lat( gsize ) )
  num_in_lat = 0
  num_in_lat(gid+1) = nobs

! Allocate memory for ODS file
! ----------------------------
  call ODS_Init  ( ods, nobs, ier, &
                   nqc=nqc, nsyn=nsyn, dist=num_in_lat )

  if ( ier .ne. 0 ) then
     rc = 4
     return
  end if

! Stop here if no sounding for this time
! --------------------------------------
  if ( nobs .eq. 0 ) then
       call NCCLOS ( id, ier ) ! ignore ier
       rc = 3
       return
  end if

! Prepare ODS file tables
! ------------------------------
  call setmeta_ ( ods%meta%kt_names, ods%meta%kt_units, ods%meta%kx_names, & 
                  ods%meta%kx_meta, ods%meta%qcx_names )

  allocate ( qcl2(nobs), land(nobs), emissmw(nobs), stat=ios )
  if ( ios .ne. 0 ) then
       rc = 4
       return
  end if

! Retrieve ODS attributes
! -----------------------
  call roms_gets_ ( nobs, start, kount, 'latitude', ods%data%lat )
  call roms_gets_ ( nobs, start, kount, 'longitude', ods%data%lon )
  call roms_gets_ ( nobs, start, kount, 'tsurf', ods%data%obs )
  call roms_gets_ ( nobs, start, kount, 'emissmw', emissmw )
  call roms_getl_ ( nobs, start, kount, 'qc', qcl2 )
  call roms_getb_ ( nobs, start, kount, 'land', land )

  if ( rc /= 0 ) then 
       call clean_()
       return
  end if

! Simple minded attribute assignment
! ----------------------------------
  ods%data%lev  = 0.0
  ods%data%kt   = KTSKINT
  ods%data%ks   = ods%data%kid
  ods%data%xm   = 0.0
  ods%data%time = 0.0

! Next reset QC flags/kx 
! ----------------------
  call setqc_ ( nobs, qcl2, land, emissmw,  ods%data%obs, &
                ods%data%qcexcl, ods%data%kx  )
  if ( rc .ne. 0 ) return
  ods%data%qchist = 0

  ftype = 'oms_tovs'

  call clean_()

! Close the file
! --------------
  call NCCLOS ( id, ier ) ! ignore ier
  if ( rc .ne. 0 ) return

! All done
! --------
  return

CONTAINS

  subroutine clean_()
  deallocate(qcl2,land)
  end subroutine clean_

!.........................................................................

  subroutine hyperslab_ ( nymd, nhms, nsyn, start, kount, rc )
  implicit NONE
  integer, intent(in)  :: nymd, nhms, nsyn
  integer, intent(out) :: start(2), kount(2)
  integer, intent(out) :: rc
!
!  Determines hyperslab parameters
!
  integer :: sounding
  integer, allocatable :: iymd(:), ihms(:), time(:)
  integer ier, did, vid, nymd1, nymd2, nhms1, nhms2, ios, jul1, jul2, i, jul
  integer time1, time2

  integer, external :: ncvid, ncdid, ods_julian, ods_caldat
  character(len=255) name

  rc = 0

! Total soundings on file
! ----------------------- 
  did = ncdid ( id, 'sounding', ier )
  call ncdinq ( id, did, name, sounding, ier ) 

! Get date/time for each sounding
! -------------------------------
  allocate(iymd(sounding), ihms(sounding), time(sounding), stat=ios )
  start(1) = 1; kount(1) = sounding
  vid = ncvid ( id, 'yyyymmdd', ier )
  call ncvgt ( id, vid, start, kount, iymd, ier )
  vid = ncvid ( id, 'hhmmss', ier )
  call ncvgt ( id, vid, start, kount, ihms, ier )

! Determine start/end dates for this period
! -----------------------------------------
  jul1 = ods_julian ( nymd )
  jul2 = jul1
  nhms1 = nhms - 120000/nsyn
  if ( nhms1 .lt. 0 ) then
       jul1 = jul1 - 1
       nhms1 = 240000 + nhms1
  endif
  nymd1 = ods_caldat ( jul1 )
  nhms2 = nhms + 120000/nsyn
  if ( nhms2 .gt. 240000 ) then
       jul2 = jul2 + 1
       nhms2 = nhms2 - 240000
  endif
  nymd2 = ods_caldat ( jul2 )

! Prepare combined date/time indices
! ----------------------------------
  time1 = (jul1-jul1)*1000000 + nhms1 ! = nhms1
  time2 = (jul2-jul1)*1000000 + nhms2
  do i = 1, sounding
     jul  = ods_julian ( iymd(i) )
     time(i) = (jul-jul1)*1000000 + ihms(i)
  end do

! Determine start for this synoptic time
! --------------------------------------
  start(1) = 0 ! no sounding in this range
  do i = 1, sounding
     if ( time(i) >= time1 .AND. time(i) < time2 ) then
          start(1) = i
          exit
     end if
  end do

! Determine how many soundings for this synoptic time
! ---------------------------------------------------
  if ( start(1) > 0 ) then
     kount(1) = 1
     do i = start(1)+1, sounding
        if ( time(i) >= time2 ) exit
        kount(1) = kount(1) + 1
     end do
  else
     kount(1) = 0
  end if

! All done
! --------
  deallocate ( iymd, ihms, time )

 end subroutine hyperslab_

!.........................................................................

      subroutine setmeta_ ( ktnames, ktunits, kxnames, kxmeta, qcxnames )

      implicit none

      character (len=*), intent(out), dimension(:) :: ktnames
      character (len=*), intent(out), dimension(:) :: ktunits
      character (len=*), intent(out), dimension(:) :: kxnames
      character (len=*), intent(out), dimension(:) :: kxmeta 
      character (len=*), intent(out), dimension(:) :: qcxnames

      qcxnames = ' '
      kxnames = ' '
      kxmeta = ' '
      ktnames = ' '
      ktunits = ' '

      ktnames   (8) = 'temperature'
      ktunits   (8) = 'k'
      ktnames   (23) = 'geopotential thickness'
      ktunits   (23) = 'm'
      ktnames   (7) = 'humidity mixing ratio'
      ktunits   (7) = 'g/kg'
      ktnames   (10) = 'humidity relative'
      ktunits   (10) = '%'
      ktnames   (11) = 'Specific humidity '
      ktunits   (11) = 'g/kg'
      ktnames   (5) = 'diff (obs-fg) specific humidity'
      ktunits   (5) = 'g/kg'
      ktnames   (4) = 'diff (obs-fg) temperature'
      ktunits   (4) = 'K'
      ktnames   (38) = 'Surface skin temperature'
      ktunits   (38) = 'K'

! Priority 2 SATELLITE
!========================================================
      kxnames   (93) = '93 DAOTOVS Land  P2 Type: Clr H/M/S'
      kxmeta    (93) = 'oms_file:daotovs_p2_lrhms.oms'
      kxnames   (94) = '94 DAOTOVS Land  P2 Type: Clr H/M'
      kxmeta    (94) = 'oms_file:daotovs_p2_lrhm.oms'
      kxnames   (95) = '95 DAOTOVS Land  P2 Type: MSU/SSU'
      kxmeta    (95) = 'oms_file:daotovs_p2_lms.oms'
      kxnames   (96) = '96 DAOTOVS Land  P2 Type: MSU'
      kxmeta    (96) = 'oms_file:daotovs_p2_lm.oms'
      kxnames   (97) = '97 DAOTOVS Land  P2 Type: SSU'
      kxmeta    (97) = 'oms_file:daotovs_p2_ls.oms'
      kxnames   (98) = '98 DAOTOVS Ocean P2 Type: Clr H/M/S'
      kxmeta    (98) = 'oms_file:daotovs_p2_orhms.oms'
      kxnames   (99) = '99 DAOTOVS Ocean P2 Type: Clr H/M'
      kxmeta    (99) = 'oms_file:daotovs_p2_orhm.oms'
      kxnames   (100) = '100 DAOTOVS Ocean P2 Type: MSU/SSU'
      kxmeta    (100) = 'oms_file:daotovs_p2_oms.oms'
      kxnames   (101) = '101 DAOTOVS Ocean P2 Type: MSU'
      kxmeta    (101) = 'oms_file:daotovs_p2_om.oms'
      kxnames   (102) = '102 DAOTOVS Ocean P2 Type: SSU'
      kxmeta    (102) = 'oms_file:daotovs_p2_os.oms'
      kxnames   (103) = '103 DAOTOVS Ice   P2 Type: Clr H/M/S'
      kxmeta    (103) = 'oms_file:daotovs_p2_irhms.oms'
      kxnames   (104) = '104 DAOTOVS Ice   P2 Type: Clr H/M'
      kxmeta    (104) = 'oms_file:daotovs_p2_irhm.oms'
      kxnames   (105) = '105 DAOTOVS Ice   P2 Type: MSU/SSU'
      kxmeta    (105) = 'oms_file:daotovs_p2_ims.oms'
      kxnames   (106) = '106 DAOTOVS Ice   P2 Type: MSU'
      kxmeta    (106) = 'oms_file:daotovs_p2_im.oms'
      kxnames   (107) = '107 DAOTOVS Ice   P2 Type: SSU'
      kxmeta    (107) = 'oms_file:daotovs_p2_is.oms'
      kxnames   (108) = '108 DAOTOVS Land  P2 Type: Cld Clr H/M/S'
      kxmeta    (108) = 'oms_file:daotovs_p2_ldhms.oms'
      kxnames   (109) = '109 DAOTOVS Land  P2 Type: Cld Clr H/M'
      kxmeta    (109) = 'oms_file:daotovs_p2_ldhm.oms'
      kxnames   (110) = '110 DAOTOVS Ocean P2 Type: Cld Clr H/M/S'
      kxmeta    (110) = 'oms_file:daotovs_p2_odhms.oms'
      kxnames   (111) = '111 DAOTOVS Ocean P2 Type: Cld Clr H/M'
      kxmeta    (111) = 'oms_file:daotovs_p2_odhm.oms'
      kxnames   (112) = '112 DAOTOVS Ice   P2 Type: Cld Clr H/M/S'
      kxmeta    (112) = 'oms_file:daotovs_p2_idhms.oms'
      kxnames   (113) = '113 DAOTOVS Ice   P2 Type: Cld Clr H/M'
      kxmeta    (113) = 'oms_file:daotovs_p2_idhm.oms'

! P1 SATELLITE
!========================================================
      kxnames   (125) = '125 DAOTOVS Land  P1 Type: Clr H/M/S'
      kxmeta    (125) = 'oms_file:daotovs_p1_lrhms.oms'
      kxnames   (126) = '126 DAOTOVS Land  P1 Type: Clr H/M'
      kxmeta    (126) = 'oms_file:daotovs_p1_lrhm.oms'
      kxnames   (127) = '127 DAOTOVS Land  P1 Type: MSU/SSU'
      kxmeta    (127) = 'oms_file:daotovs_p1_lms.oms'
      kxnames   (128) = '128 DAOTOVS Land  P1 Type: MSU'
      kxmeta    (128) = 'oms_file:daotovs_p1_lm.oms'
      kxnames   (129) = '129 DAOTOVS Land  P1 Type: SSU'
      kxmeta    (129) = 'oms_file:daotovs_p1_ls.oms'
      kxnames   (130) = '130 DAOTOVS Ocean P1 Type: Clr H/M/S'
      kxmeta    (130) = 'oms_file:daotovs_p1_orhms.oms'
      kxnames   (131) = '131 DAOTOVS Ocean P1 Type: Clr H/M'
      kxmeta    (131) = 'oms_file:daotovs_p1_orhm.oms'
      kxnames   (132) = '132 DAOTOVS Ocean P1 Type: MSU/SSU'
      kxmeta    (132) = 'oms_file:daotovs_p1_oms.oms'
      kxnames   (133) = '133 DAOTOVS Ocean P1 Type: MSU'
      kxmeta    (133) = 'oms_file:daotovs_p1_om.oms'
      kxnames   (134) = '134 DAOTOVS Ocean P1 Type: SSU'
      kxmeta    (134) = 'oms_file:daotovs_p1_os.oms'
      kxnames   (135) = '135 DAOTOVS Ice   P1 Type: Clr H/M/S'
      kxmeta    (135) = 'oms_file:daotovs_p1_orhms.oms'
      kxnames   (136) = '136 DAOTOVS Ice   P1 Type: Clr H/M'
      kxmeta    (136) = 'oms_file:daotovs_p1_irhm.oms'
      kxnames   (137) = '137 DAOTOVS Ice   P1 Type: MSU/SSU'
      kxmeta    (137) = 'oms_file:daotovs_p1_ims.oms'
      kxnames   (138) = '138 DAOTOVS Ice   P1 Type: MSU'
      kxmeta    (138) = 'oms_file:daotovs_p1_im.oms'
      kxnames   (139) = '139 DAOTOVS Ice   P1 Type: SSU'
      kxmeta    (139) = 'oms_file:daotovs_p1_is.oms'
      kxnames   (140) = '140 DAOTOVS Land  P1 Type: Cld Clr H/M/S'
      kxmeta    (140) = 'oms_file:daotovs_p1_ldhms.oms'
      kxnames   (141) = '141 DAOTOVS Land  P1 Type: Cld Clr H/M'
      kxmeta    (141) = 'oms_file:daotovs_p1_ldhm.oms'
      kxnames   (142) = '142 DAOTOVS Ocean P1 Type: Cld Clr H/M/S'
      kxmeta    (142) = 'oms_file:daotovs_p1_odhms.oms'
      kxnames   (143) = '143 DAOTOVS Ocean P1 Type: Cld Clr H/M'
      kxmeta    (143) = 'oms_file:daotovs_p1_odhm.oms'
      kxnames   (144) = '144 DAOTOVS Ice   P1 Type: Cld Clr H/M/S'
      kxmeta    (144) = 'oms_file:daotovs_p1_irhms.oms'
      kxnames   (145) = '145 DAOTOVS Ice   P1 Type: Cld Clr H/M'
      kxmeta    (145) = 'oms_file:daotovs_p1_idhms.oms'

!     ATOVS  SATELLITE AM
!========================================================
      kxnames   (186) = '186 DAOTOVS Land  P2 Type:Clr H/A/B'
      kxmeta    (186) = 'oms_file:daotovs_p2_lrhms.oms'
      kxnames   (187) = '187 DAOTOVS Land  P2 Type:Clr H/A'
      kxmeta    (187) = 'oms_file:daotovs_p2_lrhm.oms'
      kxnames   (188) = '188 DAOTOVS Land  P2 Type:AMSUA/B'
      kxmeta    (188) = 'oms_file:daotovs_p2_lms.oms'
      kxnames   (189) = '189 DAOTOVS Land  P2 Type:AMSUA'
      kxmeta    (189) = 'oms_file:daotovs_p2_lm.oms'
      kxnames   (190) = '190 DAOTOVS Land  P2 Type:AMSUB'
      kxmeta    (190) = 'oms_file:daotovs_p2_ls.oms'
      kxnames   (191) = '191 DAOTOVS Ocean P2 Type:Clr H/A/B'
      kxmeta    (191) = 'oms_file:daotovs_p2_orhms.oms'
      kxnames   (192) = '192 DAOTOVS Ocean P2 Type:Clr H/A'
      kxmeta    (192) = 'oms_file:daotovs_p2_orhm.oms'
      kxnames   (193) = '193 DAOTOVS Ocean P2 Type:AMSUA/B'
      kxmeta    (193) = 'oms_file:daotovs_p2_oms.oms'
      kxnames   (194) = '194 DAOTOVS Ocean P2 Type:AMSUA'
      kxmeta    (194) = 'oms_file:daotovs_p2_om.oms'
      kxnames   (195) = '195 DAOTOVS Ocean P2 Type: AMSUB'
      kxmeta    (195) = 'oms_file:daotovs_p2_os.oms'
      kxnames   (196) = '196 DAOTOVS Ice   P2 Type: Clr H/A/B'
      kxmeta    (196) = 'oms_file:daotovs_p2_irhms.oms'
      kxnames   (197) = '197 DAOTOVS Ice   P2 Type: Clr H/A'
      kxmeta    (197) = 'oms_file:daotovs_p2_irhm.oms'
      kxnames   (198) = '198 DAOTOVS Ice   P2 Type: AMSUA/B'
      kxmeta    (198) = 'oms_file:daotovs_p2_ims.oms'
      kxnames   (199) = '199 DAOTOVS Ice   P2 Type: AMSUA'
      kxmeta    (199) = 'oms_file:daotovs_p2_im.oms'
      kxnames   (200) = '200 DAOTOVS Ice   P2 Type: AMSUB'
      kxmeta    (200) = 'oms_file:daotovs_p2_is.oms'
      kxnames   (201) = '201 DAOTOVS Land  P2 Type:CldClr H/A/B'
      kxmeta    (201) = 'oms_file:daotovs_p2_ldhms.oms'
      kxnames   (202) = '202 DAOTOVS Land  P2 Type:CldClr H/A'
      kxmeta    (202) = 'oms_file:daotovs_p2_ldhm.oms'
      kxnames   (203) = '203 DAOTOVS Ocean P2 Type:CldClr H/A/B'
      kxmeta    (203) = 'oms_file:daotovs_p2_odhms.oms'
      kxnames   (204) = '204 DAOTOVS Ocean P2 Type:CldClr H/A'
      kxmeta    (204) = 'oms_file:daotovs_p2_odhm.oms'
      kxnames   (205) = '205 DAOTOVS Ice   P2 Type:CldClr H/A/B'
      kxmeta    (205) = 'oms_file:daotovs_p2_idhms.oms'
      kxnames   (206) = '206 DAOTOVS Ice   P2 Type:CldClr H/A'
      kxmeta    (206) = 'oms_file:daotovs_p2_idhm.oms'

!     ATOVS  SATELLITE AM
!========================================================
      kxnames   (207) = '207 DAOTOVS Land  P2 Type:Clr H/A/B'
      kxmeta    (207) = 'oms_file:daotovs_p2_lrhms.oms'
      kxnames   (208) = '208 DAOTOVS Land  P2 Type:Clr H/A'
      kxmeta    (208) = 'oms_file:daotovs_p2_lrhm.oms'
      kxnames   (209) = '209 DAOTOVS Land  P2 Type:AMSUA/B'
      kxmeta    (209) = 'oms_file:daotovs_p2_lms.oms'
      kxnames   (210) = '210 DAOTOVS Land  P2 Type:AMSUA'
      kxmeta    (210) = 'oms_file:daotovs_p2_lm.oms'
      kxnames   (211) = '211 DAOTOVS Land  P2 Type:AMSUB'
      kxmeta    (211) = 'oms_file:daotovs_p2_ls.oms'
      kxnames   (212) = '212 DAOTOVS Ocean P2 Type:Clr H/A/B'
      kxmeta    (212) = 'oms_file:daotovs_p2_orhms.oms'
      kxnames   (213) = '213 DAOTOVS Ocean P2 Type:Clr H/A'
      kxmeta    (213) = 'oms_file:daotovs_p2_orhm.oms'
      kxnames   (214) = '214 DAOTOVS Ocean P2 Type:AMSUA/B'
      kxmeta    (214) = 'oms_file:daotovs_p2_oms.oms'
      kxnames   (215) = '215 DAOTOVS Ocean P2 Type:AMSUA'
      kxmeta    (215) = 'oms_file:daotovs_p2_om.oms'
      kxnames   (216) = '216 DAOTOVS Ocean P2 Type: AMSUB'
      kxmeta    (216) = 'oms_file:daotovs_p2_os.oms'
      kxnames   (217) = '217 DAOTOVS Ice   P2 Type: Clr H/A/B'
      kxmeta    (217) = 'oms_file:daotovs_p2_irhms.oms'
      kxnames   (218) = '218 DAOTOVS Ice   P2 Type: Clr H/A'
      kxmeta    (218) = 'oms_file:daotovs_p2_irhm.oms'
      kxnames   (219) = '219 DAOTOVS Ice   P2 Type: AMSUA/B'
      kxmeta    (219) = 'oms_file:daotovs_p2_ims.oms'
      kxnames   (220) = '220 DAOTOVS Ice   P2 Type: AMSUA'
      kxmeta    (220) = 'oms_file:daotovs_p2_im.oms'
      kxnames   (221) = '221 DAOTOVS Ice   P2 Type: AMSUB'
      kxmeta    (221) = 'oms_file:daotovs_p2_is.oms'
      kxnames   (222) = '222 DAOTOVS Land  P2 Type:CldClr H/A/B'
      kxmeta    (222) = 'oms_file:daotovs_p2_ldhms.oms'
      kxnames   (223) = '223 DAOTOVS Land  P2 Type:CldClr H/A'
      kxmeta    (223) = 'oms_file:daotovs_p2_ldhm.oms'
      kxnames   (224) = '224 DAOTOVS Ocean P2 Type:CldClr H/A/B'
      kxmeta    (224) = 'oms_file:daotovs_p2_odhms.oms'
      kxnames   (225) = '225 DAOTOVS Ocean P2 Type:CldClr H/A'
      kxmeta    (225) = 'oms_file:daotovs_p2_odhm.oms'
      kxnames   (226) = '226 DAOTOVS Ice   P2 Type:CldClr H/A/B'
      kxmeta    (226) = 'oms_file:daotovs_p2_idhms.oms'
      kxnames   (227) = '227 DAOTOVS Ice   P2 Type:CldClr H/A'
      kxmeta    (227) = 'oms_file:daotovs_p2_idhm.oms'

      qcxnames = qcXnames ! from m_odsmeta

      return

  end subroutine setmeta_

!.........................................................................

  subroutine roms_gets_ ( nobs, start, kount, vname, var )
    integer, intent(in) :: nobs
    character(len=*), intent(in) :: vname
    integer, intent(out) :: start(2), kount(2)
    real, intent(out) :: var(nobs)
    integer, external :: ncvid
    integer*2 buffer(nobs)
    integer vid
    real*4 offset, scale  
    vid = ncvid ( id, trim(vname), ier )
    call ncvgt ( id, vid, start, kount, buffer, ier ); if (ier/=0) rc = 5
    call ncagt ( id, vid, 'add_offset', offset, ier ); if (ier/=0) rc = 5
    call ncagt ( id, vid, 'scale_factor', scale, ier ); if (ier/=0) rc = 5
    var = scale * buffer + offset
  end subroutine roms_gets_
  
  subroutine roms_getl_ ( nobs, start, kount, vname, var )
    integer, intent(in) :: nobs
    character(len=*), intent(in) :: vname
    integer, intent(out) :: start(2), kount(2)
    integer, intent(out) :: var(nobs)
    integer, external :: ncvid
    integer vid
    vid = ncvid ( id, trim(vname), ier )
    call ncvgt ( id, vid, start, kount, var, ier ); if (ier/=0) rc = 6
  end subroutine roms_getl_
  
  subroutine roms_getb_ ( nobs, start, kount, vname, var )
    integer, intent(in) :: nobs
    character(len=*), intent(in) :: vname
    integer, intent(out) :: start(2), kount(2)
    integer, intent(out) :: var(nobs)
    integer, external :: ncvid
    character*1 buffer(nobs)
    integer vid, i
    vid = ncvid ( id, trim(vname), ier )
    call ncvgt ( id, vid, start, kount, buffer, ier ); if (ier/=0) rc = 7
    do i = 1, nobs
       var(i) = ichar(buffer(i))
    end do
  end subroutine roms_getb_

!.........................................................................

  subroutine setqc_ ( nobs, qcl2, landflg, emissmw, tsurf, qcx, kx  )

    integer, intent(in) :: nobs
    integer, intent(in) :: qcl2(nobs)
    real,    intent(in) :: emissmw(nobs)
    real,    intent(in) :: tsurf(nobs)
    integer, intent(in) :: landflg(nobs)
    integer, intent(out) :: qcx(nobs)
    integer, intent(out) :: kx(nobs)

    integer, parameter :: onem=1000000
    integer, parameter :: oneb=1000000000
    integer, parameter :: tenm=10000000
    integer, parameter :: hunm=100000000

    integer :: quality(nobs), instrument(nobs)
    logical :: bad(nobs), clear(nobs), land(nobs), seaice(nobs), ocean(nobs)
    logical :: nohirs(nobs)


     instrument = mod(qcl2,tenm)/onem
     quality = mod(qcl2,oneb)/hunm
     nohirs  = ( (quality == 1) .OR. (quality == 2) )
     bad     = ( (qcl2 <= 0) )
     clear  = ( quality == 3 .or. quality == 5 )
     land   = landflg .eq. 1
     seaice = ( qcl2  >= oneb ) 
     ocean  = (landflg .eq. 0) .and. ( .not. seaice )

     where ( bad )  
             qcx = X_PRE_BAD
     elsewhere       
             qcx = 0 ! good observation
     end where

     kx = -1
     where ( land   .and.        clear   ) kx = 93
     where ( land   .and. (.not. clear ) ) kx = 108
     where ( ocean  .and.        clear   ) kx = 98
     where ( ocean  .and. (.not. clear ) ) kx = 110
     where ( seaice .and.        clear   ) kx = 103
     where ( seaice .and. (.not. clear ) ) kx = 112

     !don't use observations that didn't have HIRS
     !surface channels - JJ 
     !-------------------------------------------------------
     where ( land   .and.      nohirs    ) 
       kx = 95
       qcx = X_PRE_BAD
     endwhere
     where ( ocean  .and.      nohirs    ) 
       kx = 100
       qcx = X_PRE_BAD
     endwhere
     where ( seaice .and.      nohirs    ) 
       kx = 105
       qcx = X_PRE_BAD
     endwhere

     where ( kx == -1 ) 
            kx  = 93            ! why not?
            qcx = X_PRE_BAD
     end where

     !force to use sounding that had both IR and microwave
     !so had to pass emissivity checks - JJ 
     !-------------------------------------------------------
     where (instrument /= 4 .and. instrument /= 9) 
            qcx = X_PRE_BAD
     endwhere

     ! microwave emissivity check
     where (land .and. tsurf > 273. .and. emissmw < 0.85)
            qcx = X_PRE_BAD
     endwhere

     ! relax check over frozen surface, snow has lower emissivity
     where (land .and. tsurf <= 273. .and. emissmw < 0.65)
            qcx = X_PRE_BAD
     endwhere

end subroutine setqc_

  end subroutine roms_get1pe_

end module m_mpods


