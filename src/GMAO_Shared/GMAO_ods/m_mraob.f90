!-------------------------------------------------------------------------
!         NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!-------------------------------------------------------------------------
!BOP
!
! !MODULE:  m_mraob --- Implements RAOB Observation Metadata Stream 
! 
! !INTERFACE:
!
  module  m_mraob

! !USES:
  implicit none   
!
! !PUBLIC TYPES:
      PRIVATE
      PUBLIC  raob_oms     ! ODS global metadata (kx names, etc.)
      PUBLIC  raob_rc      ! Returned status codes from read/write routines
!                          !   for each attribute
! !PUBLIC MEMBER FUNCTIONS:
      PUBLIC  mraob_init   ! allocates memory for OMS vector
      PUBLIC  mraob_clean  ! deallocates mmeory for OMS vector
      PUBLIC  mraob_get    ! reads OMS from file
      PUBLIC  mraob_put    ! writes OMS to file
      PUBLIC  mraob_map    ! maps ODS xm into OMS index
      PUBLIC  mraob_types  ! return list of WMO instrument types
      PUBLIC  mraob_radcor ! return list of WMO codes for radiation corr
!
! !DESCRIPTION:
!     This module implements the RAOB Observation Metadata Stream,
!
! !REVISION HISTORY: 
!     15Nov2001  da Silva   First crack.
!     18Apr2002  Lucchesi   Added integer station ID
!     18Sep2002  C. Redder  Added character station ID (Id), station
!                           elevation (elev) and launch time (ltime) and
!                           removed integer station ID (istatn).  Made
!                           module privite.
!     10Oct2002  C. Redder  Added it_names (descriptive name for each
!                           intrument type) and rc_desc (description of
!                           all radcor codes) to the type raob_oms.  Added
!                           the publicly callable routine, nraob_types 
!EOP
!-------------------------------------------------------------------------

! RAOB metadata type
! ------------------
  integer, parameter  :: DIdlen   =    8 ! String length of WMO station ID
  integer, parameter  :: NCR_MAX  =   63 ! Max no. of char in table entries 
  integer, parameter  :: nsyn_Def =    4 ! Default for no syn time/day
  integer, parameter  :: nsta_Def = 2000 ! Default for no stations

  type raob_rc                           ! Return status code from the ODS
     integer          :: xm, Id, lat, &  !   read routine for each attribute
                                 lon, elev, itype, radcor, ltime
  end type raob_rc

  type raob_oms
     integer          :: nsta            ! No. of stations/reports
     integer          :: nvct            ! Max size of vectors
     integer          :: nsyn            ! No. synoptic time on file
     integer          :: Idlen           ! Length of WMO station ID
     integer          :: ntypes          ! No of WMO instrument types
     integer          :: nradcor         ! No of WMO radcor code
     integer          :: ncr             ! Actual no. of char in table entries
     character ( len = NCR_MAX ), pointer &
                      :: it_names (:)    ! Names of WMO intrument types.  Last
                                         !   entry is reserved for WMO type 0.
     character ( len = NCR_MAX ), pointer &
                      :: rc_descr (:)    ! Description of the WMO radcor codes.
                                         !   Last entry is for WMO code 0.
     real,    pointer :: xm       (:)    ! Metadata index
     character ( len = DIdlen  ), &
              pointer :: Id       (:)    ! WMO station Id
     real,    pointer :: lat      (:)    ! Latitude  (degrees) [-90,+90]
     real,    pointer :: lon      (:)    ! Longitude (degrees) [-180,+180]
     real,    pointer :: elev     (:)    ! Station elevation (m)
     integer, pointer :: itype    (:)    ! WMO instrument type indices
     integer, pointer :: radcor   (:)    ! WMO solar/IR correction codes 
     integer, pointer :: ltime    (:)    ! Launch time (minutes since syn time)
     type (raob_rc)   :: rc              ! Returned status codes
  end type raob_oms
  character (len=*), parameter :: &
     AttList = 'xm,Id,lat,lon,elev,itype,radcor,ltime'

  integer, external   :: ods_julian, ods_caldat

! Interfaces
! ----------
  Interface Init
     module procedure MRAOB_Init
  end Interface
  Interface Clean
     module procedure MRAOB_Clean
  end Interface
  Interface Put
     module procedure MRAOB_Put
  end Interface

  Interface Get
     module procedure MRAOB_Get
  end Interface
  Interface Map
     module procedure MRAOB_Map
  end Interface

CONTAINS

!-------------------------------------------------------------------------
!         NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE:  MRAOB_Init --- Allocates necessary memory
! 
! !INTERFACE:
  subroutine MRAOB_Init ( this, nsta,   rc, &
                          nsyn, ntypes, nradcor )
!
! !USES
  implicit NONE

! !INPUT PARAMETERS: 
      integer, intent(in)            :: nsta    ! number of stations
      integer, intent(in), OPTIONAL  :: nsyn    ! no. syn times/day
      integer, intent(in), OPTIONAL  :: ntypes  ! no. instr type
      integer, intent(in), OPTIONAL  :: nradcor ! no. radcor codes

! !INPUT/OUTPUT PARAMETERS:
      type(raob_oms), intent(inout)  :: this    ! RAOB OMS

      integer, intent(out)           :: rc      ! Error return code:
                                                ! = 0    - all is well
!
! !DESCRIPTION:
!     This routine allocates memory for the data structure of type raob_oms.
!     This routine also initializes the file parameters and the file tables.
!
! !REVISION HISTORY: 
!     07Apr2002  da Silva   First crack.
!     18Sep2002  C. Redder  Added character station ID (Id), station
!                           elevation (elev) and launch time (ltime) and
!                           removed integer station ID (istatn).
!     10Oct2002  C. Redder  Added optional argument, ntypes.  Added code
!                           to allocate and initialize tables.
!     18Oct2002  C. Redder  Fixed bug that inappropriately accesses the
!                           optional input parameters, ntypes and nradcor.
!EOP
!-------------------------------------------------------------------------

  character(len=*), parameter :: myname = 'mraob_init'
  integer ios, nnsta, nvct, nnsyn, iIdlen, nntypes, nnradcor

  rc = 0

! If optional parameters not specified, use default
! -------------------------------------------------
  if ( present ( nsyn    )) then
     nnsyn    = nsyn
  else
     nnsyn    = nsyn_Def
  end if

  if ( present ( ntypes  )) then
     nntypes  = ntypes
  else
     call MRaob_Types  ( nntypes  )
  end if

  if ( present ( nradcor )) then
     nnradcor = nradcor
  else
     call MRaob_RadCor ( nnradcor )
  end if

! Store file parameters
! ---------------------
  this % nsyn    = nnsyn
  this % Idlen   = DIdlen
  this % ntypes  = nntypes
  this % nradcor = nnradcor
  this % ncr     = NCR_Max

! Allocate memory for tables
! --------------------------
  allocate ( this % it_names ( nntypes  ), &
             this % rc_descr ( nnradcor ), &
             stat = ios )
  if ( ios .ne. 0 ) then
     rc = 1
     return
  end if

! Initialize tables
! -----------------
  call MRaob_Types  ( nntypes,  this % it_names )
  call MRaob_RadCor ( nnradcor, this % rc_descr )

! Check number of stations
! ------------------------
  nnsta = nsta
  if ( nsta < 0 ) nnsta = nsta_Def

! Allocate memory for station attributes
! --------------------------------------
  this % nsta = nnsta
  this % nvct = max ( 1, nnsta )
  nvct = this % nvct 
  allocate ( this % xm       (nvct), &
             this % Id       (nvct), &
             this % lat      (nvct), &
             this % lon      (nvct), &
             this % elev     (nvct), &
             this % itype    (nvct), &
             this % radcor   (nvct), &
             this % ltime    (nvct), &
             stat = ios )
  if ( ios .ne. 0 ) then
       rc = 2
       deallocate ( this % it_names, &
                    this % rc_descr, &
                    stat = ios ) ! ignore status code
       return
  end if

! All done
! --------
  return

end subroutine MRAOB_Init

!-------------------------------------------------------------------------
!         NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE:  MRAOB_Clean --- Deallocates memory used by RAOB OMS
! 
! !INTERFACE:
  subroutine MRAOB_Clean ( this, rc )

!
! !USES
  implicit NONE

! !INPUT/OUTPUT PARAMETERS:
      type(raob_oms), intent(inout) ::  this   ! RAOB OMS

! !OUTPUT PARAMETERS:
      integer,       intent(out)    ::  rc     ! Error return code:
                                               ! = 0    - all is well
!
! !DESCRIPTION: This routine allocates memory for RAOB OMS.
!
! !REVISION HISTORY: 
!     07Apr2002  da Silva   First crack.
!     18Sep2002  C. Redder  Added character station ID (Id), station
!                           elevation (elev) and launch time (ltime) and
!                           removed integer station ID (istatn).
!     10Oct2002  C. Redder  Added code to clean up tables
!EOP
!-------------------------------------------------------------------------

  character(len=*), parameter :: myname = 'mraob_clean'
  integer ios

  rc = 0

! Reset file parameters
! ---------------------
  this % nsta    = 0
  this % nvct    = 0
  this % nsyn    = 0
  this % Idlen   = 0
  this % ntypes  = 0
  this % nradcor = 0
  this % ncr     = 0

! Deallocate memory for tables
! --------------------------
  deallocate ( this % it_names, &
               this % rc_descr, &
               stat = ios )
  if ( ios .ne. 0 ) rc = rc + 1

! Deallocate memory for station attributes
! ----------------------------------------
  deallocate ( this % xm,     &
               this % Id,     &
               this % lat,    &
               this % lon,    &
               this % elev,   &
               this % itype,  & 
               this % radcor, &
               this % ltime,  &
               stat = ios )
  if ( ios .ne. 0 ) rc = rc + 2

! All done
! --------
  return

end subroutine MRAOB_Clean

!-------------------------------------------------------------------------
!         NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE:  MRAOB_Get --- Reads RAOB OMS for a given synoptic time
! 
! !INTERFACE:
  subroutine MRAOB_Get ( fname, nymd, nhms, this, rc )

!
! !USES
  implicit NONE

! !INPUT PARAMETERS: 
      character(len=*), intent(in)  :: fname ! OMS file name
      integer,          intent(in)  :: nymd  ! year-month-day, e.g., 19990701
      integer,          intent(in)  :: nhms  ! hour-min-sec,   e.g., 120000
!
! !OUTPUT PARAMETERS:
      type(raob_oms),   intent(inout) :: this ! RAOB OMS
      integer,          intent(out) :: rc    ! Error return code:
                                             !  0 - all is well
!
! !DESCRIPTION:
!     This routine reads metadata from a RAOB OMS file, allocates the
!     necessary memory for the OMS vector, and loads the data for the
!     synoptic time (nymd,nhms). Usually the OMS file is opened, the data
!     is read in and the file is closed. 
!
! !REVISION HISTORY:
!     07apr2002  da Silva   Initial code.
!     18Sep2002  C. Redder  Added character station ID (Id), station
!                           elevation (elev) and launch time (ltime) and
!                           removed integer station ID (istatn).
!     10Oct2002  C. Redder  Added code to retrieve list of intrument types.
!     11Dec2002  Herdies    Turned "this" inout.
!     09Jan2003  C. Redder  Commented out the call to MRaob_Clean
!EOP
!-------------------------------------------------------------------------

  character(len=*), parameter :: myname = 'mraob_get'
  integer :: id, ier, nkt, nkt1, nkx, nqc, ncr, nsyn, nobs, khms
  integer :: jday, syn_hour, ods_nget, Idlen, ntypes, nradcor
  integer :: first_jday, first_nymd, first_nhms
  logical :: fexists, ODS_EGet, idexists, itexists, rcexists

  rc = 0

! Determine file ID
! -----------------
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

! Determine no synoptic times/day
! -------------------------------
  call ODS_IGet ( id, 'nsyn',   nsyn, ier ); if ( ier .ne. 0 ) rc = 2
  if ( rc .ne. 0 ) then
     call ODS_Close ( id, ' ', ier )
     return
  end if

! ... and string length of WMO station ID
! ---------------------------------------
  idexists = ODS_EGet  ( id, 'Idlen', ier ); if ( ier .ne. 0 ) rc = 2
  if ( rc .ne. 0 ) then
     call ODS_Close ( id, ' ', ier )
     return
  end if
  Idlen = DIdlen
  if ( idexists ) then
     call ODS_IGet ( id, 'Idlen', Idlen, ier ); if ( ier .ne. 0 ) rc = 2
     if ( rc .ne. 0 ) then
        call ODS_Close ( id, ' ', ier )
        return
     end if
  end if

! ... and the number of instrument types
! --------------------------------------
  itexists = ODS_EGet  ( id, 'ntypes', ier ); if ( ier .ne. 0 ) rc = 2
  if ( rc .ne. 0 ) then
     call ODS_Close ( id, ' ', ier )
     return
  end if
  call MRaob_Types ( ntypes )
  if ( itexists ) then
     call ODS_IGet ( id, 'ntypes', ntypes, ier ); if ( ier .ne. 0 ) rc = 2
     if ( rc .ne. 0 ) then
        call ODS_Close ( id, ' ', ier )
        return
     end if
  end if

! ... and the number of radcor codes
! ----------------------------------
  rcexists = ODS_EGet  ( id, 'nradcor', ier ); if ( ier .ne. 0 ) rc = 2
  if ( rc .ne. 0 ) then
     call ODS_Close ( id, ' ', ier )
     return
  end if
  call MRaob_RadCor ( nradcor )
  if ( rcexists ) then
     call ODS_IGet ( id, 'nradcor', nradcor, ier ); if ( ier .ne. 0 ) rc = 2
     if ( rc  .ne. 0 ) then
        call ODS_Close ( id, ' ', ier )
        return
     end if
  end if

! Determine how many obs in this synoptic time
! --------------------------------------------
  jday = ODS_Julian ( nymd )
  syn_hour = nhms / 10000
  nobs = ODS_NGet ( id, jday, syn_hour, ier )
  if ( ier .ne. 0 ) then
     rc = 3
     call ODS_Close ( id, ' ', ier )
     return
  end if

! No observations, nothing left to do
! -----------------------------------
  if ( nobs == 0 ) then
     call ODS_Close ( id, ' ', ier )
     call MRAOB_Init  ( this, nobs, ier )
     return
  end if

! Allocate memory for OMS file
! ----------------------------
!  call MRAOB_Clean ( this, ier )      ! ignore ier
  call MRAOB_Init  ( this, nobs, ier, nsyn    = nsyn,   &
                                      ntypes  = ntypes, &
                                      nradcor = nradcor )
  if ( ier .ne. 0 ) then
     rc = 4
     call ODS_Close ( id, ' ', ier )
     return
  end if

! Read in OMS file tables
! -----------------------
  if ( idexists ) then
     call ODS_GetList ( id, 'it_names', ntypes,  this % it_names, ier )
     if ( ier .ne. 0 ) rc = 5
  end if
  if ( rcexists ) then
     call ODS_GetList ( id, 'rc_descr', nradcor, this % rc_descr, ier )
     if ( ier .ne. 0 ) rc = 5
  end if
  if ( rc .ne. 0 ) then
     call ODS_Close ( id, ' ', ier )
     return
  end if

! Read in OMS attributes
! ----------------------
  call ODS_GetR ( id, 'xm',     jday, syn_hour, nobs, this % xm,     ier )
  this % rc % xm       = ier
  if ( ier .ne. 0 ) rc = 6
  call ODS_GetC ( id, 'Id',     jday, syn_hour, nobs, this % Id,     ier )
  this % rc % Id       = ier
  if ( ier .ne. 0 ) rc = 6
  call ODS_GetR ( id, 'lat',    jday, syn_hour, nobs, this % lat,    ier )
  this % rc % lat      = ier
  if ( ier .ne. 0 ) rc = 6
  call ODS_GetR ( id, 'lon',    jday, syn_hour, nobs, this % lon,    ier )
  this % rc % lon      = ier
  if ( ier .ne. 0 ) rc = 6
  call ODS_GetR ( id, 'elev',   jday, syn_hour, nobs, this % elev,   ier )
  this % rc % elev     = ier
  if ( ier .ne. 0 ) rc = 6
  call ODS_GetI ( id, 'itype',  jday, syn_hour, nobs, this % itype,  ier )
  this % rc % itype    = ier
  if ( ier .ne. 0 ) rc = 6
  call ODS_GetI ( id, 'radcor', jday, syn_hour, nobs, this % radcor, ier )
  this % rc % radcor   = ier
  if ( ier .ne. 0 ) rc = 6
  call ODS_GetI ( id, 'ltime',  jday, syn_hour, nobs, this % ltime,  ier )
  this % rc % ltime    = ier
  if ( ier .ne. 0 ) rc = 6

! Close the file if we have opened it
! -----------------------------------
  call ODS_Close ( id, ' ', ier )
  if ( rc .ne. 0 ) return

! All done
! --------
  rc = 0
  return

  end subroutine mraob_get

!-------------------------------------------------------------------------
!         NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE:  MRAOB_Put --- Writes RAOB OMS for a given synoptic time
! 
! !INTERFACE:
  subroutine MRAOB_Put ( fname, nymd, nhms, this, rc, &
                         append, new )

!
! !USES
  implicit NONE

! !INPUT PARAMETERS: 
      character(len=*),  intent(in)  :: fname  ! OMS file name
      integer,           intent(in)  :: nymd   ! year-month-day, e.g. 19990701
      integer,           intent(in)  :: nhms   ! hour-min-sec,   e.g. 120000

      type(raob_oms),    intent(in)  :: this   ! RAOB OMS
      logical, optional, intent(in)  :: append ! If true, the data is appended
                                               ! to the OMS file at this time.
                                               ! Default: append = .false.

      logical, optional, intent(in)  :: new    ! If true, file will be created
                                               ! whether it exists or not
                                               ! Default: new    = .false.
!
! !OUTPUT PARAMETERS:
       integer,          intent(out) :: rc     ! Error return code:
                                               !  0 - all is well

! !DESCRIPTION:
!     This routine writes metadata to a RAOB OMS file on a given synoptic
!     time. Usually the OMS file is opened, the metadata is written to
!     the file is closed. 
!
! !REVISION HISTORY:
!     07apr2002  da Silva   Initial code.
!     18Sep2002  C. Redder  Added character station ID (Id), station
!                           elevation (elev) and launch time (ltime) and
!                           removed integer station ID (istatn).
!     10Oct2002  C. Redder  Added code to process list of intrument types.
!     10Nov2002  C. Redder  Bug fix - initialized rc at the beginning.
!EOP
!-------------------------------------------------------------------------

  character(len=*), parameter :: myname = 'mraob_put'

  integer :: ier, ier2, id, jday, syn_hour
  integer :: nobs, nsyn, Idlen, ntypes, nradcor, ncr
  integer :: first_jday, first_nymd, first_nhms
  integer :: ods_nget, NObs_Init
  logical :: appending, fexists, creating

  rc = 0

! Handle optional parameters
! --------------------------
  id = -1
  if ( present(append) ) then
     appending = append
  else
     appending = .false.
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
  nobs = this%nsta
  if ( (id .lt. 0 .and. fexists) .AND. (.not. creating) ) then
     call ODS_Open ( id, fname, 'w', ier ) ! open the file
     if ( ier .ne. 0 ) then
        rc = 1
        return
     end if

!    ... and determine the first date and time on file
!    -------------------------------------------------
     call ODS_IGet ( id, 'syn_beg:first_julian_day', first_jday, ier )
     if ( ier .ne. 0 ) then
        rc = 8
        call ODS_Close ( id, ' ', ier )
        return

     end if
     first_nymd = ODS_CalDat ( first_jday )
     first_nhms = 0

! or Create
! ---------
  else if ( id .lt. 0 .or. creating ) then

     nsyn    = this % nsyn
     Idlen   = this % Idlen
     ntypes  = this % ntypes
     nradcor = this % nradcor
     ncr     = this % ncr
     call ODS_SetParmI   ( 'nsyn',    nsyn,    ier )
     if ( ier .ne. 0 ) then
        rc = 2
        return
     end if
     call ODS_SetParmI   ( 'Idlen',   Idlen,   ier )
     if ( ier .ne. 0 ) then
        rc = 2
        return
     end if
     call ODS_SetParmI   ( 'ntypes',  ntypes,  ier )
     if ( ier .ne. 0 ) then
        rc = 2
        return
     end if
     call ODS_SetParmI   ( 'nradcor', nradcor, ier )
     if ( ier .ne. 0 ) then
        rc = 2
        return
     end if
     call ODS_SetParmI   ( 'strlen',  ncr,     ier )
     if ( ier .ne. 0 ) then
        rc = 2
        return
     end if
     call MRAOB_Create   ( id, fname, jday,    ier )
     if ( ier /= 0 ) then
        rc = 8
        return
     end if
     call ODS_PutList    ( id, 'it_names', ntypes,  this % it_names, ier )
     if ( ier /= 0 ) then
        rc = 9
        return
     end if
     call ODS_PutList    ( id, 'rc_descr', nradcor, this % rc_descr, ier )
     if ( ier /= 0 ) then
        rc = 9
        return
     end if

     first_nymd = nymd
     first_nhms = 0

  end if
 

! Determine the number of ob reports already on file
! --------------------------------------------------
  NObs_Init = ODS_NGet ( id, jday, syn_hour, ier )
  if ( ier .ne. 0 ) then
     rc = 7
     call ODS_Close ( id, ' ', ier )
     return
  end if

! If appending to this synoptic time, do the
! setup (if ob reports are already on file)
! ------------------------------------------
  if ( appending .and. NObs_Init .gt. 0 ) then

     call ODS_Append ( id, nobs, ier )
     if ( ier .ne. 0 ) then
        rc = 3
        call ODS_Close ( id, ' ', ier )
        return
     end if
  end if

! OK, now we are ready for the heavy duty writing
! -----------------------------------------------
  call ODS_PutR ( id, 'xm',     jday, syn_hour, nobs, this % xm,     ier )
  if ( ier .ne. 0 ) rc = 4
  call ODS_PutC ( id, 'Id',     jday, syn_hour, nobs, this % Id,     ier )
  if ( ier .ne. 0 ) rc = 4
  call ODS_PutR ( id, 'lat',    jday, syn_hour, nobs, this % lat,    ier )
  if ( ier .ne. 0 ) rc = 4
  call ODS_PutR ( id, 'lon',    jday, syn_hour, nobs, this % lon,    ier )
  if ( ier .ne. 0 ) rc = 4
  call ODS_PutR ( id, 'elev',   jday, syn_hour, nobs, this % elev,   ier )
  if ( ier .ne. 0 ) rc = 4
  call ODS_PutI ( id, 'itype',  jday, syn_hour, nobs, this % itype,  ier )
  if ( ier .ne. 0 ) rc = 4
  call ODS_PutI ( id, 'radcor', jday, syn_hour, nobs, this % radcor, ier )
  if ( ier .ne. 0 ) rc = 4
  call ODS_PutI ( id, 'ltime',  jday, syn_hour, nobs, this % ltime,  ier )
  if ( ier .ne. 0 ) rc = 4

! Close the file if user is not doing the file management
! -------------------------------------------------------
  call ODS_Close ( id, ' ', ier )
  if ( rc .ne. 0 ) return

! All done
! --------
  end subroutine mraob_put

!-------------------------------------------------------------------------
!         NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE:  MRAOB_Map --- Maps ODS xm into OMS index
! 
! !INTERFACE:

  subroutine MRAOB_Map ( nobs, xm, this, ioms, rc )

!
! !USES
!
  implicit NONE

! !INPUT PARAMETERS: 
!
      integer,        intent(in)  :: nobs       ! size of XM vector
      real,           intent(in)  :: xm(nobs)   ! ODS metadata index
      type(raob_oms), intent(in)  :: this       ! RAOB OMS

!
! !OUTPUT PARAMETERS:
!
      integer,        intent(out) :: ioms(nobs) ! OMS index (same size as XM)
      integer,        intent(out) :: rc         ! Error return code:
                                                !  0 - all is well

! !DESCRIPTION:
!     This routine translates the ODS metadata index XM into the OMS 
!     metadata index IOMS. With this index, one can address the OMS
!     vector directory, say, this%itype(ioms(i)).
!
! !REVISION HISTORY: 
!     07Apr2002  da Silva   Initial code.
!     25Apr2002  C Redder   Fixed compilation error (idx renamed to indx).
!
!EOP
!-------------------------------------------------------------------------

  character(len=*), parameter :: myname = 'mraob_map'
  integer, allocatable :: indx(:)
  integer ios, n, i
  real, parameter :: TOL = 0.001

  rc = 0
  n = this%nvct
  allocate ( indx(n), stat = ios )
  if ( ios .ne. 0 ) then
     rc = 1
     return
  end if
  indx = (/ ( i, i = 1, this%nvct ) /) 

! Match XM with each station
! --------------------------
  do i = 1, this%nsta
     where ( abs ( xm - this%xm(i) ) < TOL ) ioms = indx(i)
  end do

  deallocate ( indx )

  end subroutine mraob_map

!-------------------------------------------------------------------------
!         NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE:  PrepList --- Chaeck and prepare list of attributes
! 
! !INTERFACE:

  function PrepList ( List )

!
! !USES
!

! !INPUT PARAMETERS: 
      implicit NONE
      character (len=*), intent (in)  :: List      ! input list 
!
! !OUTPUT PARAMETERS:
      character (len=len(List))       :: PrepList  ! Modified list

! !DESCRIPTION:
!     Checks and prepares list of attributes by stripping out intermediate
!     blanks, and removing duplicate entries and unnecessary commas.  If a
!     list entry is invalid, then a blank string is returned.  This
!     function is low-level
!
! !REVISION HISTORY: 
!     24Sep2002  C Redder   Initial code.
!
!EOP
!-------------------------------------------------------------------------
  integer :: iCh1, iCh2, iCh2_Next, NCh1, LVar
  character (len=len(List)) :: Var
  character (len=1) :: Ch1

  PrepList = ' '
  NCh1     = len_trim ( List )
  LVar     = 0
  iCh2     = 1
  do iCh1 = 1, NCh1 + 1                ! Add list entry delimeter to end of
     Ch1 = ','                         !   list to process last entry
     if ( iCh1 .le. NCh1 ) Ch1 = List ( iCh1 : iCh1 )
     if (  Ch1 .ne. ' '  ) then        ! Skip over blanks
     if (  Ch1 .ne. ','  ) then        ! If not list entry delimiter
        LVar = LVar + 1                !    then add character to variable
        Var ( LVar : LVar ) = Ch1      !    name stored in scratch space
     else                              ! Else ...
        if ( LVar .gt. 0 ) then        !    If entry is not valid ...
           if      ( .not. in_list ( Var      ( : LVar ), &
                                     AttList )) then
              PrepList  = ' '          !    ... return with blank string
              return
           else if ( .not. in_list ( Var      ( : LVar ), &
                                     PrepList ( : iCh2 - 1 ))) then
                                       ! Skip over duplicate entries
              if ( iCh2 .gt. 1 ) then  ! If not first entry then
                 PrepList ( iCh2 : iCh2 ) = ','
                 iCh2   = iCh2 + 1     ! ... add preceding comma. 
              end if
              iCh2_Next = iCh2 + LVar  ! ... to each entry in the output list.
              PrepList ( iCh2 : iCh2_Next - 1 ) = Var ( : LVar )
              iCh2      = iCh2_Next
           end if
        end if
        LVar      = 0

     end if
     end if 
  end do

  return

  contains 

!    Internal function for determining if given entry exists in the list
!    -------------------------------------------------------------------
     logical function in_list ( Var, List )
     implicit NONE
     character (len=*), intent(in) :: Var, List
     in_list = index ( ',' // List // ',', &
                       ',' // Var  // ',' ) .gt. 0
     return
     end function in_list
  end function PrepList

!-------------------------------------------------------------------------
!         NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE:  MRaob_Types - Define list of radiosonde instrument types
!
! !INTERFACE:
!
      subroutine MRaob_Types ( NTypes, List )

! !USES:
      implicit none

! !OUTPUT PARAMETERS:
      integer,           intent (out), optional :: &
         NTypes          ! Number of entries in the table
      character (len=*), intent (out), dimension (:), optional :: &
         List            ! Discription of each radiosonde instrument type

! !DESCRIPTION:
!%    \label{MODS:RaobTypes}
!     Returns the number of WMO radiosonde instrument types and, if, desired
!     the description of each type.  WMO code of 0 is assumed to be
!     NTypes, since the lower bounds of the output array is assumed to
!     be 1.
!
!    !REVISION HISTORY:
!     21Sep2001  C Redder  Orininal code
!     10Oct2002  C Redder  Adapted from module m_RaobMeta
!     27May2003  R Todling Bug fix: continuation line set incorrectly
!
!EOP
!-------------------------------------------------------------------------
  integer :: iType, NNTypes
  integer, parameter :: TypeMax = 256   ! Maximum value for raob type
  character (len=*), parameter, dimension (TypeMax) :: TypeList = (/ &
    ! Description                                                 ! Code
    ! -----------                                                 ! ----
     'Reserved                                                 ', & ! 1
     'No radiosonde - passive target (e.g. reflector)          ', & ! 2
     'No radiosonde - active target (e.g. transponder)         ', & ! 3
     'No radiosonde - passive temperature-humidity profiler    ', & ! 4
     'No radiosonde - active temperature-humidity profiler     ', & ! 5
     'No radiosonde - active radio-acoustic sounder            ', & ! 6
     'No radiosonde - ...(reserved)                            ', & ! 7
     'No radiosonde - ...(reserved)                            ', & ! 8
     'No radiosonde - system unknown or not specified          ', & ! 9
     'VIZ type A pressure-commutated (USA)                     ', & ! 10
     'VIZ type B time-commutated (USA)                         ', & ! 11
     'RS SDC (Space Data Corporation - USA)                    ', & ! 12
     'Astor (no longer made - Australia)                       ', & ! 13
     'VIZ Mark I MICROSONDE (USA)                              ', & ! 14
     'EEC Company type 23 (USA)                                ', & ! 15
     'Elin (Austria)                                           ', & ! 16
     'Graw G. (Germany)                                        ', & ! 17
     'Reserved for allocation of radiosonde                    ', & ! 18
     'Graw M60 (Germany)                                       ', & ! 19
     'Indian Meteorological Service MK3 (India)                ', & ! 20
     'VIZ/Jin Yang Mark I MICROSONDE (South Korea)             ', & ! 21
     'Meisei RS2-80 (Japan)                                    ', & ! 22
     'Mesural FMO 1950A (France)                               ', & ! 23
     'Mesural FMO 1945A (France)                               ', & ! 24
     'Mesural MH73A (France)                                   ', & ! 25
     'Meteolabor Basora (Switzerland)                          ', & ! 26
     'AVK-MRZ (Russian Federation)                             ', & ! 27
     'Meteorit Marz2-1 (Russian Federation)                    ', & ! 28
     'Meteorit Marz2-2 (Russian Federation)                    ', & ! 29
     'Oki RS2-80 (Japan)                                       ', & ! 30
     'VIZ/Valcom type A pressure-commutated (Canada)           ', & ! 31
     'Shanghai Radio (China)                                   ', & ! 32
     'UK Met Office MK3 (UK)                                   ', & ! 33
     'Vinohrady (Czechoslovakia)                               ', & ! 34
     'Vaisala RS18 (Finland)                                   ', & ! 35
     'Vaisala RS21 (Finland)                                   ', & ! 36
     'Vaisala RS80 (Finland)                                   ', & ! 37
     'VIZ LOCATE Loran-C (USA)                                 ', & ! 38
     'Sprenger E076 (Germany)                                  ', & ! 39
     'Sprenger E084 (Germany)                                  ', & ! 40
     'Sprenger E085 (Germany)                                  ', & ! 41
     'Sprenger E086 (Germany)                                  ', & ! 42
     'AIR IS - 4A - 1680 (USA)                                 ', & ! 43
     'AIR IS - 4A - 1680 X (USA)                               ', & ! 44
     'RS MSS (USA)                                             ', & ! 45
     'Air IS - 4A - 403 (USA)                                  ', & ! 46
     'Meisei RS2-91 (Japan)                                    ', & ! 47
     'VALCOM (Canada)                                          ', & ! 48
     'VIZ MARK II (USA)                                        ', & ! 49
     'GRAW DFM-90 (Germany)                                    ', & ! 50
     'VIZ-B2 (USA)                                             ', & ! 51
     'Vaisala RS80-57H                                         ', & ! 52
     'AVK-RF93 (Russian Federation)                            ', & ! 53
    ('Reserved for allocation of radiosondes                   ', & ! 54
      iType = 1, 6),                                              & !-59
     'Vaisala RS80/MicroCora (Finland)                         ', & ! 60
     'Vaisala RS80/DifiCora or Marwin (Finland)                ', & ! 61
     'Vaisala RS80/PCCora (Finland)                            ', & ! 62
     'Vaisala RS80/Star (Finland)                              ', & ! 63
     'Orbital Sciences Corporation, Space Data Division        ', & ! 64
     'VIZ transponder radiosonde, model number 1499-520 (USA)  ', & ! 65
    ('Reserved for additional automated sounding systems       ', & ! 66-
      iType = 1, 5),                                              & ! 70
     'RS90/MicroCora (Finland)                                 ', & ! 71
     'RS90/DigiCora or Marwin (Finland)                        ', & ! 72
     'RS90/PCCora (Finland)                                    ', & ! 73
     'RS90/Star (Finland)                                      ', & ! 74
     'AVK-MRZ-ARMA (Russian Federation)                        ', & ! 75
     'AVK-RF95-ARMA (Russian Federation)                       ', & ! 76
     'GEOLINK GPSonde GL98 (France)                            ', & ! 77
    ('Reserved for additional automated sounding systems       ', & ! 78-
      iType = 1, 12),                                             & ! 89
     'Radiosonde not specified or unknownd                     ', & ! 90
     'Pressure-only radiosonde                                 ', & ! 91
     'Pressure-only radiosonde plus transponder                ', & ! 92
     'Pressure-only radiosonde plus radar reflector            ', & ! 93
     'No-pressure-only radiosonde plus transponder             ', & ! 94
     'No-pressure-only radiosonde plus radar reflector         ', & ! 95
     'Descending radiosonde                                    ', & ! 96
    ('Reserved for sounding systems with incomplete sondes     ', & ! 97-
      iType = 1, 3),                                              & ! 99
    ('Reserved                                                 ', & ! 100-
      iType = 1, 155),                                            & ! 254
     'Missing value                                            ', & ! 255
     'Reserved (Same as WMO code 0)                            '/)  ! 256

  if ( present ( NTypes )) NTypes = size ( TypeList )
  if ( present ( List   )) then
     NNTypes = size ( List )
     NNTypes = min ( size ( TypeList ), size ( List ))
     List  ( : NNTypes ) &
             = TypeList ( : NNTypes )
     if ( NNTypes .lt. size ( List )) List ( NNTypes + 1 : ) = ' '
  end if

  return
  end subroutine MRaob_Types

!-------------------------------------------------------------------------
!         NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE:  MRaob_RadCor - Define list of WMO codes for radiation correction
!
! !INTERFACE:
!
      subroutine MRaob_RadCor ( NRadCor, Descr )

! !USES:
      implicit none

! !OUTPUT PARAMETERS:
      integer,           intent (out), optional :: &
         NRadCor         ! Number of entries in the table
      character (len=*), intent (out), dimension (:), optional :: &
         Descr           ! Discription of each radiosonde instrument type

! !DESCRIPTION:
!%    \label{MODS:RaobTypes}
!     Returns the number of WMO radiosonde bias codes and, if, desired
!     the description of each type.  WMO code of 0 is assumed to be
!     NRadCor, since the lower bounds of the output array is assumed to
!     be 1.
!
!    !REVISION HISTORY:
!     10Sep2002  C Redder  Orininal code
!     27May2003  R Todling Bug fix: continuation line set incorrectly
!
!EOP
!-------------------------------------------------------------------------
  integer :: iRadCor, NNRadCor
  integer, parameter :: RadCorMax = 16    ! Maximum value for radcor code
  character (len=*), parameter, dimension (RadCorMax) :: RadCorDescr = (/ &
  ! Description                                                         ! Code
  ! -----------                                                         ! ----
   'CIMO solar corrected and CIMO infrared corrected               ', & ! 1
   'CIMO solar corrected and infrared corrected                    ', & ! 2
   'CIMO solar corrected                                           ', & ! 3
   'Solar and infrared corrected automatically by radiosonde system', & ! 4
   'Solar corrected automatically by radiosonde system             ', & ! 5
   'Solar and infrared corrected as specified by country           ', & ! 6
   'Solar corrected as specified by country                        ', & ! 7
  ('Reserved                                                       ', & ! 8-
    iRadCor = 1, 7),                                                  & ! 14
   'Missing value                                                  ', & ! 15 
   'No correction (Same as WMO code 0)                             '/)  ! 16

  if ( present ( NRadCor )) NRadCor = size ( RadCorDescr )
  if ( present ( Descr   )) then
     NNRadCor = min ( size ( RadCorDescr ), size ( Descr ))
     Descr ( : NNRadCor ) &
              = RadCorDescr ( : NNRadCor )
     if ( NNRadCor .lt. size ( Descr )) Descr ( NNRadCor + 1 : ) = ' '
  end if

  return
  end subroutine MRaob_RadCor
end module m_mraob
