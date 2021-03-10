!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !MODULE: m_mvcorF_blox - block matrix constructors of mvcorF
!
! !DESCRIPTION:
!
! !INTERFACE:

    module m_mvcorF_blox

      use m_mvcorF_nsep_bmop,only : fHHcorx, fHHcor1
      use m_mvcorF_nsep_bmop,only : fHDcorx
      use m_mvcorF_nsep_bmop,only : fDDcorx, fDDcor1

      use m_mvcorF_nsep,only : norm_HH => normHH
      use m_mvcorF_nsep,only : norm_DD => normDD

      implicit none
      private	! except

      public :: fHHcorx, fHHcor1
      public :: fHDcorx
      public :: fDDcorx, fDDcor1
      public :: norm_HH, norm_DD

! !REVISION HISTORY:
! 	06Dec01	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname='m_mvcorF_blox'

end module m_mvcorF_blox
