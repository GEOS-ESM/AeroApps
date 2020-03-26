!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !MODULE: m_ObsSmry - List a summary of observations
!
! !DESCRIPTION:
!
! !INTERFACE:

    module m_ObsSmry
      implicit none
      private	! except

      public :: ObsSmry		! The class data structure

      interface ObsSmry
	subroutine ObsSmry(lu,nobs,kx,kt)
	  implicit none
	  integer,intent(in) :: lu
	  integer,intent(in) :: nobs
	  integer,dimension(nobs),intent(in) :: kx
	  integer,dimension(nobs),intent(in) :: kt
	end subroutine ObsSmry
      end interface

! !REVISION HISTORY:
! 	17Dec99	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname='m_ObsSmry'

end module m_ObsSmry
