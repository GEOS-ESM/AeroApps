!-------------------------------------------------------------------------
!         NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!-------------------------------------------------------------------------
!BOP
!
! !MODULE:  m_roms --- Implements observation data stream vector class.
! 
! !INTERFACE:
!

    module  m_roms

! !USES:

   use m_odsmeta
   use m_ods

   implicit none
   
!
! !PUBLIC MEMBER FUNCTIONS:
!
   PUBLIC  roms_get
!
! !PUBLIC DATA MEMBERS:
!
   PUBLIC  obs_missing      ! observation missing value

!
! !DESCRIPTION:
! \label{ROMS:Mod}
!  This module ingests Retrieval OMS file into ODS vectors
!
! !REVISION HISTORY: 
!
! 15Nov2001  da Silva  First crack.
! 23Apr2002  da Silva  Fixed time index bug in hyperslab_(); the code now
!                      should give same results when reading OMS files
!                      produced with 16 or 32 PEs.
! 15Apr2003  Todling   elsewhere instead of else where as per WSawyer.
!
!EOP
!-------------------------------------------------------------------------

!  Overloaded Interfaces
!  ---------------------
   Interface ROMS_Get
       module procedure ROMS_Get1_
       module procedure ROMS_GetM_
   end Interface

CONTAINS

!-------------------------------------------------------------------------
!         NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE:  ODS_GetM_ --- Reads data from an ODS file
! 
! !INTERFACE:

  subroutine ROMS_GetM_ ( nfiles, fnames, nymd, nhms, ftypes, ods, rc, &
                          only_ts, nnsyn, iwindow, delay )

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

    integer, intent(in), OPTIONAL  :: only_Ts ! only skin temperature
    integer, intent(in), OPTIONAL  :: nnsyn   ! number of analysis in a day
    integer, intent(in), OPTIONAL  :: iwindow ! hhmmss time for obs windowing
    integer, intent(in), OPTIONAL  :: delay   ! hhmmss time-delay because fnames are 0/6/12/18

! !OUTPUT PARAMETERS:
!
                                              ! ODS file type: 'pre_anal'
                                              !  or 'post_anal'. If 'unknown'
                                              !  the could not be read
    character(len=*), intent(out)   :: ftypes(nfiles) 

    type(ods_vect), intent(inout)   ::  ods   ! ODS vector

    integer, intent(out)           ::  rc     ! Error return code:
                                              ! = 0    - all is well
                                              ! <1000  - could not read file # rc
                                              ! > 999  - rc-999 is rc from ODS_Merge 

! !DESCRIPTION:
! \label{ROMS:GetM}
!  This routine reads data from 1 or more TOVS OMS files, allocates the
!  necessary memory for the ODS vector, and loads the data for the
!  synoptic time (nymd,nhms). Check the contents of {\tt ftypes} to
!  check whether individual files could be read sucessfully. 
!
! !REVISION HISTORY: 
!
!  17Dec1999  da Silva  Initial code based on ods_getm_().
!  14Apr2003  G-P Lou   Added nnsyn for RUC support.
!  02Jul2003  Radakovich Added iwindow to allow windowing TOVS data.
!  16Sep2003  Todling    Added delay parameter to allow windowing at off-syn hours
!
!EOP
!-------------------------------------------------------------------------

  character(len=*), parameter :: myname = 'roms_getm_'

  type(ods_vect) :: tods(nfiles)   ! to hold data for each file

  integer i, ier, n 
  logical is_ods     ! whether input file is ODS

! Only T_s is supported for now
! -----------------------------
  if ( .not. present ( only_ts ) ) then
        rc = 1
        return
  end if

! Read each of the ODS files
! --------------------------
  rc = 0
  n  = 1
  do i = 1, nfiles
     call ROMS_Get1_ ( trim(fnames(i)), nymd, nhms, ftypes(i), tods(n), ier, &
                       only_ts = only_ts, nnsyn=nnsyn, iwindow=iwindow, delay=delay )
     if ( ier .eq. 0 ) then
!!!       print *, 'nobs = ', tods(n)%data%nobs, ' found on ', trim(fnames(i))
          n = n + 1
     else
!!!       print *, 'nobs = ', 0, ' found on ', trim(fnames(i)), ' rc = ',ier
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


end subroutine ROMS_GetM_

!-------------------------------------------------------------------------
!         NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE:  ROMS_Get1_ --- Reads data from a iRET OMS file
! 
! !INTERFACE:

  subroutine ROMS_Get1_ ( fname, nymd, nhms, ftype, ods, rc, &
                          only_Ts, nnsyn, iwindow, delay )

!
! !USES
!
    Implicit NONE

 !  include 'netcdf.inc'


! !INPUT PARAMETERS: 
!
    character(len=*), intent(in)   :: fname   ! ODS file name
    integer, intent(in)            :: nymd    ! year-month-day, e.g., 19990701
    integer, intent(in)            :: nhms    ! hour-min-sec,   e.g., 120000

    integer, intent(in), OPTIONAL  :: only_Ts ! only skin temperature
    integer, intent(in), OPTIONAL  :: nnsyn   ! number of analysis in a day
    integer, intent(in), OPTIONAL  :: iwindow ! hhmmss time for obs windowing 
    integer, intent(in), OPTIONAL  :: delay   ! hhmmss time-delay

!
! !OUTPUT PARAMETERS:
!
    character(len=*), intent(out)  :: ftype   ! ODS file type: 'pre_anal'
                                              !  or 'post_anal'
    type(ods_vect), intent(inout)  :: ods     ! ODS vector

    integer, intent(out)           :: rc      ! Error return code:
                                              !  0 - all is well

! !DESCRIPTION:
! \label{ROMS:Get}
!  This routine reads data from an iRET OMS file, allocates the necessary
!  memory for the ODS vector, and loads the data for the synoptic time
! (nymd,nhms). Usually the OMS file is opened, the data is read in and
!  the file is closed. 
!
! !TO DO: 
!    - revise algorithm to calculate minutes around current time.
!
! !REVISION HISTORY: 
!
!  07Sep1999  da Silva   Initial code.
!  14Apr2003  G.P. Lou   Added nnsyn to arg list; "algorithm" to calc
!                        minutes around current time. 
!  15Apr2003  Todling    Reintroduced qc changes that go lost from between
!                        versions 1.1.2.6 to 1.1.2.8
!  02Jul2003  Rakakovich Added iwindow to allowing windowing TOVS obs.
!  16Sep2003  Todling    Added delay parameter to allow windowing at off-syn hours
!
!EOP
!-------------------------------------------------------------------------

  character(len=*), parameter :: myname = 'roms_get'
  integer id, ier, nsyn, nobs, ios, nqc
  logical fexists

  integer :: start(2), kount(2) ! hyperslab indicators
  integer, external :: NCOPN
  integer, allocatable :: qcl2(:), land(:)
  real,    allocatable :: emissmw(:)

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
  if( present(nnsyn) ) then
      nsyn = nnsyn
  else
      nsyn = 4 ! 6 hourly data for now
  endif
  call hyperslab_ ( nymd, nhms, nsyn, start, kount, rc, iwindow=iwindow, delay=delay )
  if ( rc .ne. 0 ) return

! Determine how many obs in this synoptic time
! --------------------------------------------
  nobs = kount(1)
  nqc = nqcXmax

! Allocate memory for ODS file
! ----------------------------
!!!  call ODS_Clean ( ods, ier )      ! ignore ier
  call ODS_Init  ( ods, nobs, ier, &
                   nqc=nqc, nsyn=nsyn )

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

  subroutine hyperslab_ ( nymd, nhms, nsyn, start, kount, rc, iwindow, delay )
  implicit NONE
  integer, intent(in)  :: nymd, nhms, nsyn
  integer, intent(out) :: start(2), kount(2)
  integer, intent(out) :: rc
  integer, intent(in), OPTIONAL  :: iwindow ! window for TOVS obs used in analysis
  integer, intent(in), OPTIONAL  :: delay   ! delay  for TOVS obs used in analysis
!
!  Determines hyperslab parameters
!
  integer :: sounding
  integer, allocatable :: iymd(:), ihms(:), time(:)
  integer ier, did, vid, nymd1, nymd2, nhms1, nhms2, ios, jul1, jul2, i, jul
  integer time1, time2
  integer hh, hh2
  integer nsecf,nsec1,nsecwin,iwin,idel,mynhms

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
  iwin = 240000 / nsyn
  if(present(iwindow)) iwin = iwindow
  idel = 0
  if(present(delay)) idel = delay
  mynhms = nhms + idel

  jul1 = ods_julian ( nymd )
  jul2 = jul1
  nsec1 = mynhms/10000*3600 + mod(mynhms,10000)/100* 60 + mod(mynhms,100)
  nsecwin = iwin/10000*3600 + mod(iwin,10000)/100* 60 + mod(iwin,100)
  nsecf = nsec1 - nsecwin/2
  if ( nsecf < 0) then
     jul1 = jul1 - 1
     nsecf = nsecf + 86400
  endif
  nhms1 = int(nsecf/3600)*10000 + (mod(nsecf,3600 )/ 60)*100 + mod(nsecf, 60)
  nymd1 = ods_caldat ( jul1 )

  nsecf = nsec1 + nsecwin/2
  if(nsecf >= 86400) then
     nsecf=nsecf-86400
     jul2 = jul2 + 1
  endif
  nhms2 = int(nsecf/3600)*10000 + mod(nsecf,3600 )/ 60*100 + mod(nsecf, 60)
  nymd2 = ods_caldat ( jul2 )

!!!  print *, 'nymd, nhms 1: ', nymd1, nhms1
!!!  print *, 'nymd, nhms  : ', nymd,  mynhms
!!!  print *, 'nymd, nhms 2: ', nymd2, nhms2

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

!!!  if ( kount(1) == 0 ) then
!     print *, '                          nymd: ', nymd1, iymd(1), nymd2
!     print *, '                          nhms: ', nhms1, ihms(1), nhms2
!  else
!     print *, 'start, kount = ', start(1), kount(1)
!!!  end if

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
!!!    print *, trim(vname), ' = ', minval(var), maxval(var), offset, scale
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
!!!    print *, trim(vname), ' = ', minval(var), maxval(var)
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
!!!    print *, trim(vname), ' = ', minval(var), maxval(var)
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

  end subroutine roms_get1_


end module m_roms

