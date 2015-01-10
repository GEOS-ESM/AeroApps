module voecH_tbl
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!-----------------------------------------------------------------------
!
! !ROUTINE: voecH_tbl - a F90 module of voecH input tables
!
! !INTERFACE:	(to do)
!
! !DESCRIPTION:
!
! !EXAMPLES:
!
! !BUGS:
!
! !SEE ALSO:
!
! !SYSTEM ROUTINES:
!
! !FILES USED:
!
! !REVISION HISTORY:
! 	10Jan96 - J. Guo	- programmed and added the prolog
!	19Sep96 - J. Guo	- add `save' to all public variables
!_______________________________________________________________________
use config, only : lvmax_vc,MX_voecH
implicit none
public

	! observation error correlation functions

  integer,save :: n_voecH

  character*16,save :: name_voecH(MX_voecH)
  character*16,save :: type_voecH(MX_voecH)
  character*64,save :: desc_voecH(MX_voecH)

  integer,save :: nlev_voecH(MX_voecH)
  real,   save :: plev_voecH(lvmax_vc,MX_voecH)
  real,   save :: corr_voecH(lvmax_vc,lvmax_vc,MX_voecH)

contains

  subroutine voecH_tbl0()
  end subroutine voecH_tbl0

end module voecH_tbl
