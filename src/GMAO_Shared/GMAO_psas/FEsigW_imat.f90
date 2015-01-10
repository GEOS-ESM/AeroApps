module FEsigW_imat
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!-----------------------------------------------------------------------
!
! !ROUTINE: FEsigW_imat - F90 module of FEsigW_imat type definition
!
! !INTERFACE: (to do)
! !DESCRIPTION: (to do)
! !EXAMPLES: (to do)
! !BUGS: (to do)
! !SEE ALSO: (to do)
! !SYSTEM ROUTINES: (to do)
!
! !REVISION HISTORY:
! 	23Apr96 - J. Guo	- (to do)
!	19Sep96 - J. Guo	- add `save' to all public variables
!_______________________________________________________________________
use config, only : MXveclev,MXveclat
implicit none
private	! except
public	:: FEsigS_imat, FEsigV_imat, FEsigW_imat0

  real, save :: FEsigS_imat(MXveclev,MXveclat)	! sigma_(u,v) rotational
  real, save :: FEsigV_imat(MXveclev,Mxveclat)	! sigma_(u,v) divegent

contains

  subroutine FEsigW_imat0
  end subroutine FEsigW_imat0

end module FEsigW_imat
!.
