module hfecH_tbl
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!-----------------------------------------------------------------------
!
! !ROUTINE: hfecH_tbl - a F90 module of hfecH input tables
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
use config, only : lvmax_hc,MXpar_hc,MX_fecH
implicit none
private	! except

public	:: name_hfecH
public	:: type_hfecH,desc_hfecH
public	:: nlev_hfecH,plev_hfecH
public	:: npar_hfecH,pars_hfecH

  character(len=17),parameter	:: name_hfecH(MX_fecH) =&
    (/	'FcstErr*hCor_HH::',				&
	'FcstErr*hCor_SS::',				&
	'FcstErr*hCor_VV::' 	/)

  character*16,save :: type_hfecH(MX_fecH)	! a class name
  character*64,save :: desc_hfecH(MX_fecH)	! a description

  integer,save :: nlev_hfecH(MX_fecH)		! number of levels
  real,   save :: plev_hfecH(lvmax_hc,MX_fecH)	! levels defined

  integer,save :: npar_hfecH(MX_fecH)		! number of parameters
  real,   save :: pars_hfecH(MXpar_hc,lvmax_hc,MX_fecH)	! parameters
!-----------------------------------------------------------------------
end module hfecH_tbl
!.
