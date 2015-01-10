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
       
!       PRIVATE
!
! !PUBLIC MEMBER FUNCTIONS:
!
       PUBLIC  Init,   ODSM_Init
       PUBLIC  Clean,  ODSM_Clean

       Interface Init
          module procedure ODSM_Init
       end Interface

       Interface Clean
          module procedure ODSM_Clean
       end Interface

! !PUBLIC TYPES:
!
       PUBLIC  ods_meta            ! ODS global metadata (kx names, etc.)

! !PUBLIC DATA MEMBERS:

       include 'odsmeta.h'

!
!  ODS global metadata (kx names, etc.)
!  ------------------------------------
       integer, PUBLIC, parameter :: NKT_MAX  =  60
       integer, PUBLIC, parameter :: NKX_MAX  = 300
       integer, PUBLIC, parameter :: NQC_MAX  =  60
       integer, PUBLIC, parameter :: NCR_MAX  =  32

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
         subroutine ODSM_Init ( m, rc,       
     &                          nkt, nkx,  nqc, ncr ) ! Optional
         
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
         if ( associated ( m%kt_names ) .or.   
     &        associated ( m%kt_units ) .or.  
     &        associated ( m%kx_names ) .or.  
     &        associated ( m%kx_meta  ) .or.  
     &        associated ( m%qcx_names ) )  
     &   then
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

         allocate ( m%kt_names(nnkt),   
     &        m%kt_units(nnkt),   
     &        m%kx_names(nnkx),   
     &        m%kx_meta(nnkx),    
     &        m%qcx_names(nnqc),  
     &        stat = ios )
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

         if ( ( ios1 .ne. 0 ) .or.  
     &        ( ios2 .ne. 0 ) .or.  
     &        ( ios3 .ne. 0 ) .or.  
     &        ( ios4 .ne. 0 ) .or.  
     &        ( ios5 .ne. 0 ) )     
     &   then
            rc = rc + 1
         end if

         return

       end subroutine  odsm_clean

      end module m_MPodsmeta

