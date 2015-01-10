!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !MODULE: rlev_imat - a module of level entry table
!
! !DESCRIPTION:
!
! !INTERFACE:

    module rlev_imat

      use config, only : MXveclev
      implicit none
      private	! except

      public :: MXveclev
      public :: nveclev
      public :: pveclev
      public :: rlev_imat0

! !REVISION HISTORY:
! 	07Jan99 - Jing Guo <guo@thunder> - initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname='rlev_imat'

  integer,save :: nveclev
  real	 ,save :: pveclev(MXveclev)

  CONTAINS

    subroutine rlev_imat0
    end subroutine rlev_imat0

!.
end module rlev_imat
