module vfecH_tbl
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!-----------------------------------------------------------------------
!
! !ROUTINE: vfecH_tbl - a F90 module of vfecH input tables
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
use config, only : lvmax_vc,MX_fecH
implicit none
private			! except
public	:: name_vfecH
public	:: type_vfecH,desc_vfecH
public	:: nlev_vfecH,plev_vfecH
public	:: corr_vfecH

  character(len=17),parameter	:: name_vfecH(MX_fecH) =&
    (/	'FcstErr*vCor_HH::',				&
	'FcstErr*vCor_SS::',				&
	'FcstErr*vCor_VV::' 	/)

  character*16,save :: type_vfecH(MX_fecH)
  character*64,save :: desc_vfecH(MX_fecH)

  integer,save :: nlev_vfecH(MX_fecH)
  real,   save :: plev_vfecH(lvmax_vc,MX_fecH)
  real,   save :: corr_vfecH(lvmax_vc,lvmax_vc,MX_fecH)
!-----------------------------------------------------------------------
end module vfecH_tbl
!.
