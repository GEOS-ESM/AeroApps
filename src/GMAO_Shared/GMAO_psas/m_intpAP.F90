!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !MODULE: m_intpAP - An interface module of intpAPmiss and intpAPelem
!
! !DESCRIPTION:
!
! !INTERFACE:

    module m_intpAP
      use m_intpAPmiss,only : intpAP
      use m_intpAPelem,only : intpAP
      implicit none
      private	! except

      public :: intpAP

! !REVISION HISTORY:
!	25Apr00	- Jing Guo
!		. initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname='m_intpAP'

end module m_intpAP
