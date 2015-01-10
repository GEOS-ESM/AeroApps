module hoecH_tbl
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!-----------------------------------------------------------------------
!
! !ROUTINE: hoecH_tbl - a F90 module of hoecH input tables
!
! !INTERFACE: (to do)
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
use config, only : lvmax_hc, MXpar_hc, MX_hoecH
implicit none
public

  integer,save :: n_hoecH	! a counter upto MX_hoecH

  character*16,save :: name_hoecH(MX_hoecH)	! classes
  character*16,save :: type_hoecH(MX_hoecH)	! functional form
  character*64,save :: desc_hoecH(MX_hoecH)	! descriptions

  integer,save :: nlev_hoecH(MX_hoecH)
  integer,save :: npar_hoecH(MX_hoecH)

  real,save :: plev_hoecH(lvmax_hc,MX_hoecH)
  real,save :: pars_hoecH(MXpar_hc,lvmax_hc,MX_hoecH)

contains

  subroutine hoecH_tbl0()
  end subroutine hoecH_tbl0

end module hoecH_tbl
