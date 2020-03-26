!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !MODULE: m_kt_corF - correlation matrix generators for Fcst.Err.
!
! !DESCRIPTION:
!
! !INTERFACE:

    module m_kt_corF
      use m_mvcorF_blox,only : fHHcor1
      use m_mvcorF_blox,only : fHHcorx
      use m_mvcorF_blox,only : fDDcor1
      use m_mvcorF_blox,only : fDDcorx
      use m_mvcorF_blox,only : fHDcorx

      use m_kt_uvcorF,only : fQQcor1
      use m_kt_uvcorF,only : fQQcorx

      implicit none
      private	! except

      public :: fHHcor1
      public :: fHHcorx
      public :: fDDcor1
      public :: fDDcorx
      public :: fHDcorx

      public :: fQQcor1
      public :: fQQcorx

! !REVISION HISTORY:
!	29Aug00	- Jing Guo
!		- initial prototype/prolog/code
!		- combined fHHcor1.F etc. into this module
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname='m_kt_corF'

end module m_kt_corF

