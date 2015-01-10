!-------------------------------------------------------------------------
!         NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!-------------------------------------------------------------------------
!BOP
!
! !MODULE: m_odsmeta --- Implements observation data stream metadata class.
!
! !INTERFACE:
!

   module  m_MPodsmeta

! !USES:

   implicit none

   PRIVATE
!
! !PUBLIC MEMBER FUNCTIONS:
!
   PUBLIC  Init,  ODSM_Init
   PUBLIC  Clean, ODSM_Clean

   Interface Init
      module procedure ODSM_Init
   end Interface

   Interface Clean
      module procedure ODSM_Clean
   end Interface


! !PUBLIC DATA MEMBERS:
                                                ! Internal dimensions
   integer, PUBLIC, parameter :: NKT_MAX  =  60
   integer, PUBLIC, parameter :: NKX_MAX  = 300
   integer, PUBLIC, parameter :: NQC_MAX  =  60
   integer, PUBLIC, parameter :: NCR_MAX  =  32
                                                ! surface variables
   integer, PUBLIC, parameter  :: ktus    =  1
   integer, PUBLIC, parameter  :: ktvs    =  2
   integer, PUBLIC, parameter  :: ktslp   =  3
   integer, PUBLIC, parameter  :: ktuu    =  4
   integer, PUBLIC, parameter  :: ktvv    =  5
   integer, PUBLIC, parameter  :: ktHH    =  6
   integer, PUBLIC, parameter  :: ktww    =  7  ! CAUTION: used to be ktqq
                                                !   won't work w/ PSAS
   integer, PUBLIC, parameter  :: ktmax   = 41
   integer, PUBLIC, parameter  :: ktskint = 38
   integer, PUBLIC, parameter  :: kxmax =  512  ! datasource (kx) param.
   integer, PUBLIC, parameter :: X_PRE_BAD = 1  ! preprocessing exclusion flag   
   integer, PUBLIC, parameter :: nqcXmax = 33   ! number of exclusion flags in use

!                             -------------------------
!                             Data Type (kt) Parameters
!                             --------------------------
!
                                                ! upper-air variables


        integer, PUBLIC, parameter  :: ktTT    =  8
        integer, PUBLIC, parameter  :: ktTd    =  9
        integer, PUBLIC, parameter  :: ktrh    = 10
        integer, PUBLIC, parameter  :: ktqq    = 11     ! CAUTION: see ktww above
                                                ! More surface variables
        integer, PUBLIC, parameter  :: ktus10  = 12
        integer, PUBLIC, parameter  :: ktTs10  = 13
        integer, PUBLIC, parameter  :: ktTds   = 14
        integer, PUBLIC, parameter  :: ktrhs   = 15
        integer, PUBLIC, parameter  :: ktqs10  = 16

        integer, PUBLIC, parameter  :: ktpr    = 17
        integer, PUBLIC, parameter  :: kttpw   = 18
        integer, PUBLIC, parameter  :: kttlw   = 19
        integer, PUBLIC, parameter  :: ktfcc   = 20
                                                ! Chemistry
        integer, PUBLIC, parameter  :: kttco3  = 21
        integer, PUBLIC, parameter  :: kto3    = 22
        integer, PUBLIC, parameter  :: ktuthk  = 23

        integer, PUBLIC, parameter  :: ktspd2m   = 24
        integer, PUBLIC, parameter  :: ktmxspd2m = 25
        integer, PUBLIC, parameter  :: ktwgust2m = 26
        integer, PUBLIC, parameter  :: ktt2m     = 27
        integer, PUBLIC, parameter  :: ktmxt2m   = 28
        integer, PUBLIC, parameter  :: ktmmnt2m  = 29
        integer, PUBLIC, parameter  :: ktdewt2m  = 30
        integer, PUBLIC, parameter  :: ktrh2m    = 31
        integer, PUBLIC, parameter  :: ktsphu2m  = 32
        integer, PUBLIC, parameter  :: ktps2m    = 33
        integer, PUBLIC, parameter  :: ktvis     = 34
        integer, PUBLIC, parameter  :: ktsdepth  = 35
        integer, PUBLIC, parameter  :: ktwcond   = 36
        integer, PUBLIC, parameter  :: ktth_UppA = 37
        integer, PUBLIC, parameter  :: ktSST     = 39
        integer, PUBLIC, parameter  :: ktTb      = 40
        integer, PUBLIC, parameter  :: ktTv      = 41

        integer, PUBLIC, dimension(22), parameter :: ktSurfAll = (/ ktus, ktvs, ktslp,     &
                                             ktus10, ktTs10, ktTds, ktrhs, ktqs10, &
                                             ktspd2m, ktmxspd2m, ktwgust2m, ktt2m, &
                                             ktmxt2m, ktmmnt2m, ktdewt2m, ktrh2m,  &
                                             ktsphu2m, ktps2m, ktskint, ktpr,      &
                                             kttpw, ktSST /)

        integer, PUBLIC, dimension(12), parameter :: ktUppaAll = (/ ktuu, ktvv, ktHH, ktww, &
                                                      ktTT, ktTd, ktrh, ktqq,       &
                                                      ktuthk, ktth_UppA,            &
                                                      ktTb, ktTv /)

        integer, PUBLIC, dimension( 2), parameter :: ktChemAll = (/ kttco3, kto3 /)

        integer, PUBLIC, dimension( 5), parameter :: ktOthrAll = (/ kttlw, ktfcc, ktvis, &
                                                      ktsdepth, ktwcond /)

!                             ------------------
!                             QC Flag parameters
!                             ------------------
!
!     History flags
!     -------------

      integer, PUBLIC, parameter :: H_PRE_SUSP   =   1  ! unspecified preprocessing history flag
      integer, PUBLIC, parameter :: H_UNDERG     =   3  ! suspect due to underground check

                                                ! NCEP quality marks that we treat as suspect
      integer, PUBLIC, parameter :: H_NCEP1      =   4  !   NCEP CQC - modified value
      integer, PUBLIC, parameter :: H_NCEP3      =   5  !   NCEP CQC - suspect value
      integer, PUBLIC, parameter :: H_NCEP89     =   6  !   NCEP PREVENTS - P suspect
                                                !     (underground or sfcP too high)
                                                !      or suspect specific humid.
      integer, PUBLIC, parameter :: H_NCEP14     =   7  !   NCEP SDM purged
      integer, PUBLIC, parameter :: H_NCEP_OTHER =   8  !   NCEP missing or unknown QM
                                                !     (QM != 0,1,2,3,8,9,13,14,15)
                                                !  Suspect mark derived from
                                                !     pressure value (obs maybe OK,
                                                !     pressure is suspect)
      integer, PUBLIC, parameter :: H_NCEP_PRES1 =   9  ! NCEP CQC modified P value
      integer, PUBLIC, parameter :: H_NCEP_PRES3 =  10  ! NCEP CQC suspect P value
      integer, PUBLIC, parameter :: H_NCEP_PRES14=  11  ! NCEP SDM purged P
                                                !  Suspect marks for moisture
                                                !  vars calcfrom suspect input values

      integer, PUBLIC, parameter :: H_SUSP_TEMP    = 12 ! moisture from suspect temp
      integer, PUBLIC, parameter :: H_SUSP_DEWTEMP = 13 ! moisture from suspect dewpt

      integer, PUBLIC, parameter :: H_BACKG   = 17      ! background check

      integer, PUBLIC, parameter :: H_YELLOW  = 20      ! obs marked as suspect by "Yellow List"
      integer, PUBLIC, parameter :: H_SIMUL   = 21      ! obs confidence level less than 1
      integer, PUBLIC, parameter :: H_BEAUFORT = 22     ! Beaufort corrected winds


!     Exclusion flags
!     ---------------
      integer, PUBLIC, parameter :: X_BAD_LOC     =   2 ! bad location flag
      integer, PUBLIC, parameter :: X_UNDERG      =   3 ! underground flag
      integer, PUBLIC, parameter :: X_OBS_FILL    =   4 ! observation fill flag
      integer, PUBLIC, parameter :: X_FCS_FILL    =   5 ! forecast at obs location fill flag
      integer, PUBLIC, parameter :: X_TOO_HIGH    =   6 ! obs above certain desired level
      integer, PUBLIC, parameter :: X_PASSIVE     =   7 ! passive data type exclusion flag
      integer, PUBLIC, parameter :: X_TIME_ACT    =   8 ! data outside active time window
      integer, PUBLIC, parameter :: X_NOT_ANAVAR  =   9 ! not an analysis variable


                                                ! NCEP quality marks treated as exclusion
      integer, PUBLIC, parameter :: X_NCEP13      =  10 ! NCEP CQC bad observation
      integer, PUBLIC, parameter :: X_NCEP15      =  11 ! NCEP PREPDATA bad observation
      integer, PUBLIC, parameter :: X_NCEP_PRES13 =  12 ! NCEP CQC bad pressure
      integer, PUBLIC, parameter :: X_NCEP_PRES15 =  13 ! NCEP PREPDATA bad pressure

                                                ! DAO preprocessing quality marks
      integer, PUBLIC, parameter :: X_RANGE       =  14 ! DAO range check failed
      integer, PUBLIC, parameter :: X_DUP         =  15 ! DAO duplicate obs. (>1 in 6 hr)
      integer, PUBLIC, parameter :: X_HYDRO       =  16 ! DAO failed hydrostatic check
                                                !  Rejection marks for moisture values
                                                !    calculated from rejected/dubious input


      integer, PUBLIC, parameter :: X_BUDDY       = 17  ! buddy check
      integer, PUBLIC, parameter :: X_WIND        = 18  ! wind check
      integer, PUBLIC, parameter :: X_NOSTATS     = 19  ! no error statistics
      integer, PUBLIC, parameter :: X_PROFILE     = 20  ! profile check
      integer, PUBLIC, parameter :: X_BACKG       = 21  ! background check
      integer, PUBLIC, parameter :: X_BADXM       = 22  ! bad xm for scaling

      integer, PUBLIC, parameter :: X_BADLAYER    = 25  ! layer too thin or too thick

      integer, PUBLIC, parameter :: X_BAD_TEMP    = 28  ! moisture from bad temperature
      integer, PUBLIC, parameter :: X_THIN        = 29  ! observation excl. by thinner

      integer, PUBLIC, parameter :: X_UNPHYSICAL  = 30  ! obs with unphysical value
      integer, PUBLIC, parameter :: X_RED         = 31  ! obs excluded by "Red List"
      integer, PUBLIC, parameter :: X_SIMUL       = 32  ! obs could not be simulated
      integer, PUBLIC, parameter :: X_PSAS        = 33  ! excluded by PSAS

!     Descriptions of history flags:
!     -----------------------------

                                ! number of history flags in use
      integer, PUBLIC, parameter :: nqcHmax = 22

      character(len=32), PUBLIC, parameter :: qcHnames(nqcHmax)=(/ &
                      'unspecified preprocessing flag  ', &
                      '                                ', & 
                      'gcm slightly underground        ', &
                      'NCEP CQC - modified value       ', &
                      'NCEP CQC - suspect value        ', &
                      'NCEP PREVENTS undergnd/bad P, q ', &
                      'NCEP SDM purged                 ', &
                      'NCEP unknown or missing QM      ', &
                      'NCEP CQC modified P             ', &
                      'NCEP CQC suspect P value        ', &
                      'NCEP SDM purged P               ', &
                      'moisture from suspect temp.     ', &
                      'moisture from suspect dewpt.    ', &
                      '                                ', &
                      '                                ', &
                      '                                ', &
                      'outlier wrt background          ', &
                      '                                ', &
                      '                                ', &
                      'marked suspect by Yellow list   ', & 
                      'obs simulated w/ confidence < 1 ', &
                      'Beaufort Corrected Winds        '     /)


!     Descriptions of exclusion flags:
!     -------------------------------

      character(len=33), PUBLIC, parameter :: qcXnames(nqcXmax)=(/ &
                      'unspecified preprocessing flag  ', &
                      'impossible location             ', &
                      'gcm deep underground            ', &
                      'observation value undefined     ', &
                      'forecast value undefined        ', &
                      'observation level too high      ', &
                      'passive data type               ', & 
                      'outside active time window      ', &
                      'not an analysis variable        ', & 
                      'NCEP CQC bad observation        ', &
                      'NCEP PREPDATA bad observation   ', &
                      'NCEP CQC bad pressure           ', &
                      'NCEP PREPDATA bad pressure      ', &
                      'DAO range check failed          ', &
                      'DAO duplicate obs. (>1 in 6 hr) ', & 
                      'DAO failed hydrostatic check    ', &
                      'failed buddy check              ', &
                      'incomplete wind vector          ', &
                      'improper error statistics       ', &
                      'incomplete vertical profile     ', &
                      'extreme outlier wrt background  ', & 
                      'bad xm value for omf scaling    ', &
                      '                                ', &
                      '                                ', &
                      'layer too thin or too thick     ', &
                      '                                ', &
                      '                                ', & 
                      'moisture from bad temp.         ', &
                      'observation removed by thinner  ', &
                      'unphysical value                ', &
                      'excluded by "Red List"          ', &
                      'obs could not be simulated      ', &
                      'excluded by PSAS                '/)

! !PUBLIC TYPES:
!
   PUBLIC  ods_meta            ! ODS global metadata (kx names, etc.)
!  ODS global metadata (kx names, etc.)
!  ------------------------------------
   type ods_meta

      integer :: nkt  ! actual number of KT's
      integer :: nkx  ! actual number of KX's
      integer :: nqc  ! actual number of QC flags
      integer :: ncr  ! actual number of chars in tables

      character(len=NCR_MAX), pointer :: kt_names(:)
      character(len=NCR_MAX), pointer :: kt_units(:)
      character(len=NCR_MAX), pointer :: kx_names(:)
      character(len=NCR_MAX), pointer :: kx_meta(:)
      character(len=NCR_MAX), pointer :: qcx_names(:)

   end type ods_meta

!
! !DESCRIPTION:
!
!  This module defines the observation data stream metadata class.
!  It relies on the ODS library and HDF.
!
! !REVISION HISTORY:
!
!  10Apr93 - Jing G.       - Initial code (kt_max.h)
!  24Mar95 - Jing G.       - New version for text based data table.
!  09Mar98 - Todling       - Added all kt definition from ON-95.01
!  10Apr98 - Todling       - Updated kt according to Y.Kondratyeva
!  29Dec99 - Todling       - Embedded in ODSmeta module.
!
!  10Jun2002 Zaslavsky Created by splitting m_ods module into m_odsmeta, 
!                      m_odsdata, and m_ods modules. Two separate classes 
!                      (for metadata and data) were created. The content of
!                      m_odsmeta.f file (metadata parameters) was included
!                      in this file.
!  17Jun2002 Zaslavsky Generic name Init for ODSM_Init and Clean for ODSM_Clean
!                      were introduced.
!  17Jun2002 Zaslavsky Public Data Members described.
!  18Jun2002 Zaslavsky nsyn moved from this module to m_odsdata module.
!EOP
!-------------------------------------------------------------------------

CONTAINS

!-------------------------------------------------------------------------
!         NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE:  ODSM_Init --- Allocate memory for ODS metadata.
!
! !INTERFACE:
!
  subroutine ODSM_Init ( m, rc,       &
                         nkt, nkx,  nqc, ncr ) ! Optional

! !USES:

  implicit none

! !INPUT PARAMETERS:
                                           ! Size of tables and atributes:
  integer, intent(in), optional :: nkt     !  number of KT's
  integer, intent(in), optional :: nkx     !  number of KX's
  integer, intent(in), optional :: nqc     !  number of QC flags
  integer, intent(in), optional :: ncr     !  number of chars in ODS tables


! !INPUT/OUTPUT PARAMETERS:

  type(ods_meta), intent(inout) :: m       ! ODS metadata to be allocated

! !OUTPUT PARAMETERS:

  integer, intent(out)          ::  rc     ! Error return code:
                                           !  0 - all is well
                                           !  1 - error allocating metadata

! !DESCRIPTION:
!
!  Allocates memory for an ODS metadata.
!
! !REVISION HISTORY:
!
!  10Jun2002   Zaslavsky  Created using existing ODS_Init procedure.
!  18Jun2002   Zaslavsky  Added checks to avoid memory leaks.
!
!EOP
!-------------------------------------------------------------------------

  character(len=*), parameter :: myname = 'odsm_init'
  integer ios
  integer nnkt, nnkx, nnqc, nncr


! Avoid memory leaks
! ------------------
  if ( associated ( m%kt_names ) .or. & 
       associated ( m%kt_units ) .or. &
       associated ( m%kx_names ) .or. &
       associated ( m%kx_meta  ) .or. &
       associated ( m%qcx_names ) )   &
  then
       rc = 1 
       return 
  end if


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


! Allocate memory for ODS global metadata
! ---------------------------------------
  m%nkt  = nnkt
  m%nkx  = nnkx
  m%nqc  = nnqc
  m%ncr  = nncr

  allocate ( m%kt_names(nnkt),   &
             m%kt_units(nnkt),   &
             m%kx_names(nnkx),   &
             m%kx_meta(nnkx),    &
             m%qcx_names(nnqc),  &
             stat = ios )
  if ( ios .ne. 0 ) then
       rc = 1
  else
       rc = 0
  end if

  return

  end subroutine odsm_init



!-------------------------------------------------------------------------
!         NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE:  ODSM_Clean - Deallocates ODS metadata
!
! !INTERFACE:
!
      subroutine ODSM_Clean ( m, rc )

! !USES:

      implicit none

! !INPUT/OUTPUT PARAMETERS:

  type(ods_meta), intent(inout) ::  m      ! ODS metadata to be deallocated

! !OUTPUT PARAMETERS:

  integer, intent(out)          ::  rc     ! Error return code:
                                           !  0 - all is well
                                           !  1 - error deallocating metadata

! !DESCRIPTION:
!
!  Frees memory used by an ODS vector.
!
! !REVISION HISTORY:
!
!  10Jun2002   Zaslavsky  Created using existing ODS_Clean procedure
!  18Jun2002   Zaslavsky  Added checks to avoid memory leaks.
!  23Oct2002   Sawyer     Bug fix:  added deallocate( kt_units )
!
!EOP
!-------------------------------------------------------------------------

  character(len=*), parameter :: myname = 'odsm_clean'
  integer :: ios1 = 0, ios2 = 0, ios3 = 0, ios4 = 0, ios5 = 0

  rc = 0

! Deallocate memory for ODS global metadata
! -----------------------------------------
  m%nkt = 0
  m%nkx = 0
  m%nqc = 0
  m%ncr = 0
  if ( associated ( m%kt_names ) ) deallocate (  m%kt_names, stat = ios1 )
  if ( associated ( m%kx_names ) ) deallocate (  m%kx_names, stat = ios2 )
  if ( associated ( m%kx_meta  ) ) deallocate (  m%kx_meta, stat = ios3  )
  if ( associated (m%qcx_names ) ) deallocate (  m%qcx_names, stat = ios4)
  if ( associated ( m%kt_units ) ) deallocate (  m%kt_units, stat = ios5 )

  if ( ( ios1 .ne. 0 ) .or. &
       ( ios2 .ne. 0 ) .or. &
       ( ios3 .ne. 0 ) .or. &
       ( ios4 .ne. 0 ) .or. &
       ( ios5 .ne. 0 ) )    &
  then
         rc = rc + 1
  end if

  return

  end subroutine  odsm_clean

end module m_MPodsmeta

