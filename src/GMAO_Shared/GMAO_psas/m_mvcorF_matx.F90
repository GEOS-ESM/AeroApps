!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !MODULE: m_mvcorF_matx - block mat-vec operators of mvcorF
!
! !DESCRIPTION:
!
! !INTERFACE:

    module m_mvcorF_matx

      use m_mvcorF_nsep_bmop,only : sHH1cxpy, xHH1cxpy, rHH1cxpy
      use m_mvcorF_nsep_bmop,only : xHD1cxpy, rHD1cxpy, rDH1cxpy
      use m_mvcorF_nsep_bmop,only : sDD1cxpy, xDD1cxpy, rDD1cxpy

      use m_mvcorF_nsep_bmop,only : sHHmcxpy, xHHmcxpy, rHHmcxpy
      use m_mvcorF_nsep_bmop,only : xHDmcxpy, rHDmcxpy, rDHmcxpy
      use m_mvcorF_nsep_bmop,only : sDDmcxpy, xDDmcxpy, rDDmcxpy

      implicit none
      private	! except

      public :: sHH1cxpy, xHH1cxpy, rHH1cxpy
      public :: xHD1cxpy, rHD1cxpy, rDH1cxpy
      public :: sDD1cxpy, xDD1cxpy, rDD1cxpy

      public :: sHHmcxpy, xHHmcxpy, rHHmcxpy
      public :: xHDmcxpy, rHDmcxpy, rDHmcxpy
      public :: sDDmcxpy, xDDmcxpy, rDDmcxpy

! !REVISION HISTORY:
! 	06Dec01	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname='m_mvcorF_matx'

end module m_mvcorF_matx
