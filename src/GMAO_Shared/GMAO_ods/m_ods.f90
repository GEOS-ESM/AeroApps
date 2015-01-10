!-------------------------------------------------------------------------
!         NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!-------------------------------------------------------------------------
!BOP
!
! !MODULE:  m_ods --- Implements observation data stream vector class.
!
! !INTERFACE:
!

   module  m_ods

! !USES:

   use m_odsmeta
   implicit none

! !PUBLIC TYPES:
!
   PRIVATE
   PUBLIC  ods_meta        ! ODS global metadata (kx names, etc.)
   PUBLIC  obs_vect        ! atomic attributes (lat, lon, etc.)
   PUBLIC  ods_vect        ! ODS vector consisting of ods_meta and obs_vect
!
! !PUBLIC MEMBER FUNCTIONS:
!
   PUBLIC  ods_init
   PUBLIC  ods_clean
   PUBLIC  ods_put
   PUBLIC  ods_get
   PUBLIC  ods_Merge
   PUBLIC  ods_MaskOut
   PUBLIC  ods_Tally
   PUBLIC  ods_Select
   PUBLIC  obs_Reorder
   PUBLIC  obs_MoveUp
!
! !PUBLIC DATA MEMBERS:
!
   PUBLIC  obs_missing      ! observation missing value

!
! !DESCRIPTION:
! \label{MODS:Mod}
!  This module defines the observation data stream vector class.
!  It relies on the ODS library and HDF.
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
!  10May2002 Dee       Added OBS_Reorder, OBS_MoveUp
!                         (based on GEOS-2 ODS_Reorder, ODS_Moveup)
!  17Oct2002 C. Redder Increased NCR_MAX (string length in kt_names, etc)
!                      from 32 to 255.
!  28Apr2004 Todling   Teached ODS_Get to read diag_conv(GSI-out) files
!  20Dec2004 Dee       Increased NKX_MAX
!  14Sep2007 Meta      Increased NKT_MAX to match odsmeta.h
!
!EOP
!-------------------------------------------------------------------------

!  Internal dimensions
!  -------------------
   integer, parameter :: NKT_MAX  =  128
   integer, parameter :: NKX_MAX  = 1000
   integer, parameter :: NQC_MAX  =  60
   integer, parameter :: NCR_MAX  = 255
   integer, parameter :: NSYN_DEF =   4
   integer, parameter :: NOBS_MAX = 450 * 1000

   real, parameter :: obs_missing = 1.0E+15


!  ODS global metadata (kx names, etc.)
!  ------------------------------------
   type ods_meta

      integer :: nkt  ! actual number of KT's
      integer :: nkx  ! actual number of KX's
      integer :: nqc  ! actual number of QC flags
      integer :: ncr  ! actual number of chars in tables
      integer :: nsyn ! actual number of synoptic times per day

      character(len=NCR_MAX), pointer :: kt_names(:)
      character(len=NCR_MAX), pointer :: kt_units(:)
      character(len=NCR_MAX), pointer :: kx_names(:)
      character(len=NCR_MAX), pointer :: kx_meta(:)
      character(len=NCR_MAX), pointer :: qcx_names(:)

   end type ods_meta


!  Observation atomic atttributes (lat, lon, lev, etc.)
!  ----------------------------------------------------
   type obs_vect

      integer          :: nobs      ! actual number of observations
      integer          :: nvct      ! vector size allocated (nobs .le. nvct)

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

!  ODS vector
!  ----------
   type  ods_vect

      type(ods_meta) meta           ! ODS file metadata (kx_names, etc.)
      type(obs_vect) data           ! Data vectors (lat, lon, obs, omf, etc.)

   end type ods_vect

!  Overloaded Interfaces
!  ---------------------
   Interface ODS_Get
       module procedure ODS_Get1_
       module procedure ODS_GetM_
   end Interface

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
                        nkt, nkx,  nqc, ncr, nsyn ) ! Optional

! !USES:

  implicit none

! !INPUT/OUTPUT PARAMETERS:

  type(ods_vect), intent(inout) ::  ods    ! ODS vector to be allocated

  integer, intent(inout)        :: nobs    ! number of observations; if
                                           ! nobs<0  nobs will be reset to
                                           !         a large internal value.

                                           ! Size of tables and atributes:
  integer, intent(inout), optional :: nkt  !  number of KT's
  integer, intent(inout), optional :: nkx  !  number of KX's
  integer, intent(inout), optional :: nqc  !  number of QC flags
  integer, intent(inout), optional :: ncr  !  number of chars in ODS tables
  integer, intent(inout), optional :: nsyn !  number of synoptic times per
                                           !    day
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
!
!EOP
!-------------------------------------------------------------------------

  character(len=*), parameter :: myname = 'ods_init'
  integer ios, nvct, i
  integer nnkt, nnkx, nnqc, nncr, nnsyn


! If optional parameters not specified, use default
! -------------------------------------------------
  if ( present(nkt ) ) then
       nnkt  = nkt
  else
       nnkt  = NKT_MAX
  end if
  if ( present(nkx ) ) then
       nnkx  = nkx
  else
       nnkx  = NKX_MAX
  end if
  if ( present(nqc ) ) then
       nnqc  = nqc
  else
       nnqc  = NQC_MAX
  end if
  if ( present(ncr ) ) then
       nncr  = ncr
  else
       nncr  = NCR_MAX
  end if
  if ( present(nsyn) ) then
       nnsyn = nsyn
  else
       nnsyn = NSYN_DEF
  end if

! Reset nobs if user specifies less than zero value
! -------------------------------------------------
  if ( nobs .lt. 0 ) nobs = NOBS_MAX


! Allocate memory for ODS global metadata
! ---------------------------------------
  ods%meta%nkt  = nnkt
  ods%meta%nkx  = nnkx
  ods%meta%nqc  = nnqc
  ods%meta%ncr  = nncr
  ods%meta%nsyn = nnsyn
  allocate ( ods%meta%kt_names(nnkt),   &
             ods%meta%kt_units(nnkt),   &
             ods%meta%kx_names(nnkx),   &
             ods%meta%kx_meta(nnkx),    &
             ods%meta%qcx_names(nnqc),  &
             stat = ios )
  if ( ios .ne. 0 ) then
       rc = 1
       return
  end if

! Allocated memory for ODS atomic attributes
! ------------------------------------------
  if ( nobs .gt. 0 ) then
       nvct = nobs
  else
       nvct = 1
  endif
  ods%data%nobs = nobs
  ods%data%nvct = nvct

  allocate    ( ods%data%kid(nvct),      &
                ods%data%lat(nvct),      &
                ods%data%lon(nvct),      &
                ods%data%lev(nvct),      &
                ods%data%kx(nvct),       &
                ods%data%kt(nvct),       &
                ods%data%ks(nvct),       &
                ods%data%xm(nvct),       &
                ods%data%time(nvct),     &
                ods%data%obs(nvct),      &
                ods%data%OmF(nvct),      &
                ods%data%OmA(nvct),      &
                ods%data%Xvec(nvct),     &
                ods%data%qcexcl(nvct),   &
                ods%data%qchist(nvct),   &
                stat = ios )
     if ( ios .ne. 0 ) then
        rc = 2
        return
     end if

! Set observation index
! ---------------------
  ods%data%kid = (/ (i, i=1,nvct) /)

! All done
! --------
  rc = 0
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

      implicit none

! !INPUT/OUTPUT PARAMETERS:

  type(ods_vect), intent(inout) ::  ods    ! ODS vector to be allocated

! !OUTPUT PARAMETERS:

  integer, intent(out)          ::  rc     ! Error return code:
                                           !  0 - all is well
                                           !  1 - error deallocating metadata
                                           !  2 - error deallocating attributes

! !DESCRIPTION:
! \label{MODS:Clean}
!  Frees memory used by an ODS vector.
!
! !REVISION HISTORY:
!
!  07sep1999   da Silva   Initial code
!  16Feb2001   Todling    Changed to account for Xvec addition
!  28Mar2001   Todling    Made it go through whole clean even in error.
!  25Mar2003   Todling    Added associated checks.
!
!EOP
!-------------------------------------------------------------------------

  character(len=*), parameter :: myname = 'ods_clean'
  logical alld
  integer ios

  rc = 0

! Deallocate memory for ODS global metadata
! -----------------------------------------
  ods%meta%nkt = 0
  ods%meta%nkx = 0
  ods%meta%nqc = 0
  ods%meta%ncr = 0
  alld = associated(ods%meta%kt_names ) .and. associated(ods%meta%kt_units) .and. &
         associated(ods%meta%kx_names ) .and. associated(ods%meta%kx_meta ) .and. &
         associated(ods%meta%qcx_names)
  if ( alld ) then
    deallocate ( ods%meta%kt_names,   &
                 ods%meta%kt_units,   &
                 ods%meta%kx_names,   &
                 ods%meta%kx_meta,    &
                 ods%meta%qcx_names,  &
                 stat = ios )
  else
    ios = 1
  end if
  if ( ios .ne. 0 ) rc = rc + 1

! Deallocated memory for ODS atomic attributes
! --------------------------------------------
  ods%data%nobs = 0
  ods%data%nvct = 0
  alld = associated(ods%data%kid)    .and. &
         associated(ods%data%lat)    .and. associated(ods%data%lon)    .and. &
         associated(ods%data%lev)    .and. associated(ods%data%kx )    .and. &
         associated(ods%data%kt )    .and. associated(ods%data%ks )    .and. &
         associated(ods%data%xm )    .and. associated(ods%data%time)   .and. &
         associated(ods%data%obs)    .and. associated(ods%data%OmF)    .and. &
         associated(ods%data%OmA)    .and. associated(ods%data%Xvec)   .and. &
         associated(ods%data%qcexcl) .and. associated(ods%data%qchist)
  if ( alld ) then
    deallocate ( ods%data%kid,      &
                 ods%data%lat,      &
                 ods%data%lon,      &
                 ods%data%lev,      &
                 ods%data%kx,       &
                 ods%data%kt,       &
                 ods%data%ks,       &
                 ods%data%xm,       &
                 ods%data%time,     &
                 ods%data%obs,      &
                 ods%data%OmF,      &
                 ods%data%OmA,      &
                 ods%data%Xvec,     &
                 ods%data%qcexcl,   &
                 ods%data%qchist,   &
                 stat = ios )
  else
    ios = 1
  end if
  if ( ios .ne. 0 ) rc = rc + 1

  return

  end subroutine  ods_clean


!-------------------------------------------------------------------------
!         NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE:  ODS_GetM_ --- Reads data from an ODS file
!
! !INTERFACE:

  subroutine ODS_GetM_ ( nfiles, fnames, nymd, nhms, ftypes, ods, rc )

!
! !USES
!
    Implicit NONE

! !INPUT PARAMETERS:
!
    integer, intent(in)            :: nfiles  ! number of files
    character(len=*), intent(in)   :: fnames(nfiles) ! ODS file names
    integer, intent(in)            :: nymd    ! year-month-day, e.g., 19990701
    integer, intent(in)            :: nhms    ! hour-min-sec,   e.g., 120000

! !OUTPUT PARAMETERS:
!
                                              ! ODS file type: 'pre_anal'
                                              !  or 'post_anal'. If 'unknown'
                                              !  the could not be read
    character(len=*), intent(out)   :: ftypes(nfiles)

    type(ods_vect), intent(inout)  ::  ods    ! ODS vector

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
! !REVISION HISTORY:
!
!  03Dec1999  da Silva  Initial code.
!  09Mar2000  da Silva  Now is OK if some of the files could not be read,
!                       the only requiment is that at least one file can
!                       read.
!  12Dec2999  Redder    Revised prologue to reflect the added flexibility
!                       of arbitrarily setting the ODS file parameter,
!                       nsyn.
!  11Dec2002  Herdies   Turned "ods" inout.
!
!EOP
!-------------------------------------------------------------------------

  character(len=*), parameter :: myname = 'ods_getm_'

  type(ods_vect) :: tods(nfiles)   ! to hold data for each file

  integer i, ier, n


! Read each of the ODS files
! --------------------------
  rc = 0
  n  = 1
  do i = 1, nfiles
     call ODS_Get1_ ( trim(fnames(i)), nymd, nhms, ftypes(i), tods(n), ier )
     if ( ier .eq. 0 ) then
          n = n + 1
     else
          ftypes(i) = 'unknown'
     end if
  end do
  n = n - 1  ! actual number of files that could be read

! At least one file could be read
! -------------------------------
  if ( n .le. 0 ) then
     rc = 1
     return
  end if


! Merge ODS files
! ---------------
  call ODS_Merge ( tods, n, ods, ier )
  if ( ier .ne. 0 ) then
       rc = 999 + ier
       return
  end if

! Free memory used to hold data for individual files
! --------------------------------------------------
  do i = 1, n
     call ODS_Clean ( tods(i), ier )  ! ignore ier
  end do

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

  subroutine ODS_Get1_ ( fname, nymd, nhms, ftype, ods, rc, &
                         ncid, ncf )               ! Optional

!
! !USES
!
    Implicit NONE

! !INPUT PARAMETERS:
!
    character(len=*), intent(in)   :: fname   ! ODS file name
    integer, intent(in)            :: nymd    ! year-month-day, e.g., 19990701
    integer, intent(in)            :: nhms    ! hour-min-sec,   e.g., 120000

    integer, intent(in), optional  :: ncid    ! ODS file id as returned by
                                              ! ODS_Open(); in this case
                                              ! fname is not used
    logical, intent(in), optional  :: ncf     ! Non-compliant format, when .t.
                                              !  assumes this input file not 
                                              !  ODS format; rather GSI-output
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
! \label{MODS:Get1}
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
!  14Jun2002   RT/Redder Revised exit status for nobs=0 case
!  11Dec2002   Herdies   Turned "ods" inout; removed clean for init.
!  25Jun2003   Todling   Defining nsyn only in case it's been read ok.
!  28Apr2004   Todling   Overload of ODS to read GSI-diag-out files
!  14Jun2004   Todling   Added ncf opt-arg to avoid errmsg from Open.
!                        
!
!EOP
!-------------------------------------------------------------------------

  character(len=*), parameter :: myname = 'ods_get1_'
  integer id, ier, nkt, nkt1, nkx, nqc, ncr, nsyn, nobs, khms
  integer jday, syn_hour, ods_nget
  integer first_jday, first_nymd, first_nhms
  logical fexists, dconv

  rc = 0
  dconv = .false.
  if ( present(ncf) ) then
       if(ncf) dconv=ncf
  endif

! Determine file ID
! -----------------
  if ( present(ncid) ) then
       id = ncid
  else
       inquire ( file=trim(fname), exist=fexists )
       if ( fexists ) then
          if ( .not. dconv ) then 
               call ODS_Open ( id, fname, 'r', ier ) ! open the file
               if ( ier .ne. 0 ) dconv = .true.
          endif

       else
          rc = 1
          return
       end if
  end if

! Try to read as if file were diag_conv instead of legitimate ODS
! ---------------------------------------------------------------
  if ( dconv ) then
       call ods_dcget ( fname, nymd, nhms, ods, ier )
         if ( ier/= 0 ) then
              rc = 1
              return
         end if
       ftype = 'post_anal'
       rc = 0
       return  ! if it gets here, it read diag_conv successfully
  end if

! Determine size of ODS tables
! ----------------------------
  call ODS_IGet ( id, 'nkt',    nkt,  ier ); if ( ier .ne. 0 ) rc = 2
  call ODS_IGet ( id, 'nkx',    nkx,  ier ); if ( ier .ne. 0 ) rc = 2
  call ODS_IGet ( id, 'nqcx',   nqc,  ier ); if ( ier .ne. 0 ) rc = 2
  call ODS_IGet ( id, 'strlen', ncr,  ier ); if ( ier .ne. 0 ) rc = 2
  call ODS_IGet ( id, 'nsyn',   nsyn, ier ); if ( ier .ne. 0 ) rc = 2
  if ( rc .ne. 0 ) then
     if ( .not. present(ncid) ) call ODS_Close ( id, ' ', ier )
     return
  end if
  ods%meta%nsyn = nsyn

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

! No observations, nothing left to do
! -----------------------------------
  if ( nobs == 0 ) then
     if ( .not. present(ncid) ) call ODS_Close ( id, ' ', ier )
!ams call ODS_Clean ( ods, ier )      ! ignore ier
     call ODS_Init  ( ods, nobs, ier, nsyn=nsyn )
     if ( ier .ne. 0 ) rc = 4
     return
  end if

! Allocate memory for ODS file
! ----------------------------
!Dirceu  call ODS_Clean ( ods, ier )      ! ignore ier
  call ODS_Init  ( ods, nobs, ier, &
                   nkt=nkt1, nkx=nkx, nqc=nqc, ncr=ncr, nsyn=nsyn )
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

  end subroutine ods_get1_

!-------------------------------------------------------------------------
!         NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE:  ODS_Put --- Writes to an ODS file
!
! !INTERFACE:

  subroutine ODS_Put ( fname, ftype,  nymd, nhms, ods, rc, &
                       ncid,  append, new )   ! Optional
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
! \label{MODS:Put}
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
!  10Jan2003   Redder   - Bug fixes: rc wasn't defined from start;
!                                    fix bounds in time array
!
!EOP
!-------------------------------------------------------------------------

  character(len=*), parameter :: myname = 'ods_put'

  integer :: ier, ier2, jday, nobs, syn_hour, nsyn, id
  integer :: first_jday, first_nymd, first_nhms
  integer :: ods_nget, NObs_Init
  logical :: appending, fexists, creating
  integer, dimension (:), allocatable :: time  ! scratch space

  rc = 0
  
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
       nsyn = ods%meta%nsyn
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
    time(1:nobs) = ods%data%time(1:nobs) 
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
  end subroutine ods_put


!-------------------------------------------------------------------------
!         NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE:  ODS_Merge --- Merges 2 or more ODS vectors
!
! !INTERFACE:
!
      subroutine ODS_Merge ( ods_in, nods, ods_out, rc )

! !USES:

      implicit none

! !INPUT PARAMETERS:

  integer,        intent(in)   ::  nods          ! number of input ODS vectors

   type(ods_vect), intent(in)  ::  ods_in(nods)  ! Many ODS vectors

! !OUTPUT PARAMETERS:

  type(ods_vect), intent(inout)::  ods_out       ! Merge ODS vector
  integer,        intent(out)  ::  rc            ! Error return code:
                                                 !  0 - all is well
                                                 !  1 - could not allocate mem

! !DESCRIPTION:
! \label{MODS:Merge}
!  Merges 2 or more ODS vectors.
!
! !BUGS:
!
!     The merged vector has the same metadata as the first vector, that is
!  the metadata is not really merged.
!
!
! !REVISION HISTORY:
!
!  03dec1999   da Silva   Initial code
!  23dec1999   da Silva   Fixed metadata size bug.
!  26jul2000   da  Silva  Removed kid from merging since kid is set
!                         to 1:nvct during ods_init().
!  16Feb2001   Todling    Changed to account for Xvec addition
!  04Oct2001   G Lou      Added nsyn for RUC
!  11Dec2002   Herdies    Turned "ods_out" inout.
!
!EOP
!-------------------------------------------------------------------------

  character(len=*), parameter :: myname = 'ods_merge'

  integer i, ier, n, n1, n2, ks_offset
  integer nkt, nkx, nqc, ncr, nobs, nnsyn

  rc = 0

! Get total number of observations
! --------------------------------
  nobs = 0
  do i = 1, nods
     nobs = nobs + ods_in(i)%data%nobs
  end do


! Metadata sizes
! --------------
  nkt=maxval(ods_in(1:nods)%meta%nkt)
  nkx=maxval(ods_in(1:nods)%meta%nkx)
  nqc=maxval(ods_in(1:nods)%meta%nqc)
  ncr=maxval(ods_in(1:nods)%meta%ncr)


! Allocate memory for merged ODS file
! -----------------------------------
  nnsyn = ods_in(1)%meta%nsyn
  call ODS_Clean ( ods_out, ier )      ! ignore ier
  call ODS_Init  ( ods_out, nobs, ier, nkt, nkx, nqc, ncr, nsyn=nnsyn )
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
     ods_out%data%obs(n1:n2)    = ods_in(i)%data%obs(1:n)
     ods_out%data%OmF(n1:n2)    = ods_in(i)%data%OmF(1:n)
     ods_out%data%OmA(n1:n2)    = ods_in(i)%data%OmA(1:n)
     ods_out%data%Xvec(n1:n2)   = ods_in(i)%data%Xvec(1:n)
     ods_out%data%qcexcl(n1:n2) = ods_in(i)%data%qcexcl(1:n)
     ods_out%data%qchist(n1:n2) = ods_in(i)%data%qchist(1:n)

!    Make sure souding index is unique
!    ---------------------------------
     ods_out%data%ks(n1:n2)     = ks_offset + ods_in(i)%data%ks(1:n)

     ks_offset = ks_offset + maxval ( ods_in(i)%data%ks(1:n) )

  end do

! All done
! --------
  return

end subroutine ODS_Merge



!-------------------------------------------------------------------------
!         NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE:  ODS_MaskOut - Mask out observations given
!
! !INTERFACE:
!
      subroutine ODS_MaskOut ( y, nobs, boxes, nboxes, X_FLAG, qcx, nset, &
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
! \label{MODS:MaskOut}
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
!
!EOP
!-------------------------------------------------------------------------

  character(len=*), parameter :: myname = 'ods_maskout'

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

end subroutine ODS_MaskOut

!..............................................................

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

  implicit none

! !INPUT PARAMETERS:

  type(ods_vect), intent(inout) :: ods        ! ods vector
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

  type(ods_vect), intent(inout), optional :: odss   ! selected obs vector
  integer, intent(out)                  :: nsel     ! number of selected obs
  integer, intent(out)                  :: rc       ! return code:
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
!  11Dec2002   Herdies    Turned "odss" inout.
!  14Jul2003   RT/AMS     Added nobs=0 bug fix that got lost from 1.1.2.17.2.1
!                         also extended bug fix to initial nobs=0 check
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

   if ( nobs .eq. 0 ) then
        nsel = 0
        if ( present(odss) ) then

!          allocate memory for output ods structure
!          ----------------------------------------
           nnkt = ods%meta%nkt
           nnkx = ods%meta%nkx
           nnqc = ods%meta%nqc
           nncr = ods%meta%ncr
           nnsyn = ods%meta%nsyn
           call ODS_Init ( odss, nsel, ier, nnkt, nnkx, nnqc, nncr, nsyn=nnsyn )
        end if
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
        if ( present(odss) ) then

!          allocate memory for output ods structure
!          ----------------------------------------
           nnkt = ods%meta%nkt
           nnkx = ods%meta%nkx
           nnqc = ods%meta%nqc
           nncr = ods%meta%ncr
           nnsyn = ods%meta%nsyn
           call ODS_Init ( odss, nsel, ier, nnkt, nnkx, nnqc, nncr, nsyn=nnsyn )
        end if
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
      nnsyn = ods%meta%nsyn
      call ODS_Init ( odss, nsel, ier, nnkt, nnkx, nnqc, nncr, nsyn=nnsyn )
      if (ier .ne. 0) then
          deallocate ( mask, selected )
          rc = 1
          return
      end if

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
! !IROUTINE:  OBS_MoveUp() --- Move selected observations to front.
!
! !INTERFACE:
!
      subroutine OBS_MoveUp ( y, nobs, nmoved, rc,  &
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
!  10May2002 (Dee)   Initial code based on GEOS-2 ODS_MoveUp
!  14May2002 (Dee)   Added optional output parameter idx
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
           call OBS_Reorder ( nobs, y, indx )
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

      end subroutine OBS_MoveUp

!-------------------------------------------------------------------------
!         NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE:  OBS_Reorder() --- Reorder observation vector.
!
! !INTERFACE:
!
      subroutine OBS_Reorder ( nobs, y, indx )

!
! !INPUT PARAMETERS:
!
      integer , intent(in)    ::  nobs          ! number of obs
      integer , intent(in)    ::  indx(nobs)    ! sorting index
!
! !INPUT/OUTPUT PARAMETERS:
!
      type ( obs_vect ) , intent (inout) :: y
!
! !DESCRIPTION:  Reorders the first nobs elements of each array
!                attribute of the observation vector y according
!                to the index vector indx
!
! !REVISION HISTORY:
!
!  10May2002 (Dee)   Initial code based on GEOS-2 ODS_Reorder
!
!EOP
!-------------------------------------------------------------------------

      character(len=*), parameter :: myname='OBS_Reorder'

      integer i

!     Nothing to do
!     -------------
      if ( nobs==0 ) return

!     Apply sorting permutation to data arrays
!     ----------------------------------------
      y%kid    (1:nobs) = y%kid    ( (/ (indx(i), i=1,nobs) /) )
      y%kt     (1:nobs) = y%kt     ( (/ (indx(i), i=1,nobs) /) )
      y%kx     (1:nobs) = y%kx     ( (/ (indx(i), i=1,nobs) /) )
      y%ks     (1:nobs) = y%ks     ( (/ (indx(i), i=1,nobs) /) )
      y%lon    (1:nobs) = y%lon    ( (/ (indx(i), i=1,nobs) /) )
      y%lat    (1:nobs) = y%lat    ( (/ (indx(i), i=1,nobs) /) )
      y%lev    (1:nobs) = y%lev    ( (/ (indx(i), i=1,nobs) /) )
      y%time   (1:nobs) = y%time   ( (/ (indx(i), i=1,nobs) /) )
      y%obs    (1:nobs) = y%obs    ( (/ (indx(i), i=1,nobs) /) )
      y%OmF    (1:nobs) = y%OmF    ( (/ (indx(i), i=1,nobs) /) )
      y%OmA    (1:nobs) = y%OmA    ( (/ (indx(i), i=1,nobs) /) )
      y%Xvec   (1:nobs) = y%Xvec   ( (/ (indx(i), i=1,nobs) /) )
      y%xm     (1:nobs) = y%xm     ( (/ (indx(i), i=1,nobs) /) )
      y%qcexcl (1:nobs) = y%qcexcl ( (/ (indx(i), i=1,nobs) /) )
      y%qchist (1:nobs) = y%qchist ( (/ (indx(i), i=1,nobs) /) )

      end subroutine OBS_Reorder

end module m_ods

