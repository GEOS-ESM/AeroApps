!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !MODULE: m_FEsigFi_tabl - the data module of sigFi resource
!
! !DESCRIPTION:
!
! !INTERFACE:

    module m_FEsigFi_tabl
      implicit none
      private	! except

      public :: FEsigFi_Name		! The class data structure
      public :: FEsigFi_clean

      interface FEsigFi_clean; module procedure clean_; end interface

! !REVISION HISTORY:
! 	07Jan99 - Jing Guo <guo@thunder> - initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname='m_FEsigFi_tabl'

  character(len=120),save :: FEsigFi_Name=' '

contains
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: clean_ - clean the contents of the data module
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine clean_()
      implicit none

! !REVISION HISTORY:
! 	07Jan99 - Jing Guo <guo@thunder> - initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::clean_'

  FEsigFi_Name=' '

end subroutine clean_
end module m_FEsigFi_tabl
