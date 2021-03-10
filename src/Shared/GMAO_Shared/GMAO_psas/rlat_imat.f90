module rlat_imat
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!-----------------------------------------------------------------------
!
! !ROUTINE: rlat_imat - a lookup table of latitude for indexed matrices
!
! !INTERFACE:
!	use rlat_imat[, only : nveclat, veclats]
!
! !DESCRIPTION:
!	Module rlat_imat defines a lookup table of latitude values for
!	indexed matrix.
!
! !EXAMPLES: (to do)
! !BUGS: (to do)
! !SEE ALSO: (to do)
! !REVISION HISTORY:
! 	04Mar96 - J. Guo	- (to do)
!	19Sep96 - J. Guo	- add `save' to all public variables
!_______________________________________________________________________

  use config, only : MXveclat
  implicit none

  private	!  except ...
  public	:: nveclat, veclats
  public	:: MXveclat
  public        :: rlat_imat0

  integer,save :: nveclat		! actual size of veclats(:)
  real,   save :: veclats(MXveclat)	! a lookup table

Contains

  subroutine rlat_imat0
  end subroutine rlat_imat0

end module rlat_imat
!.
