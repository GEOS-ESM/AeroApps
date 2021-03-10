!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !MODULE: m_mvcorF_bldr - Builder of mvcorF operators
!
! !DESCRIPTION:
!
!   This module builds mvcorF block matrix operators.
!
! !INTERFACE:

    module m_mvcorF_bldr
      use m_mvcorF_nsep_bldr,only : mvcorF_init  => mvcorF_nsep_init
      use m_mvcorF_nsep_bldr,only : mvcorF_clean => mvcorF_nsep_clean
      implicit none
      private	! except

      public :: mvcorF_init
      public :: mvcorF_clean

! !REVISION HISTORY:
! 	06Dec01	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname='m_mvcorF_bldr'

end module m_mvcorF_bldr
