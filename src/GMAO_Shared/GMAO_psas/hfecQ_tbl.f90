module hfecQ_tbl
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!-----------------------------------------------------------------------
!
! !ROUTINE: hfecQ_tbl - a F90 module of hfecQ input tables
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
use config, only : lvmax_hc,MXpar_hc
implicit none
public

  character(len=*),parameter ::	name_hfecQ='FcstErr*hCor_QQ::'

  character*16,save :: type_hfecQ	! a class name
  character*64,save :: desc_hfecQ	! a description

  integer,save :: nlev_hfecQ		! number of levels
  real,   save :: plev_hfecQ(lvmax_hc)	! levels defined

  integer,save :: npar_hfecQ			! number of parameters
  real,   save :: pars_hfecQ(MXpar_hc,lvmax_hc)	! parameters defined

Contains

  subroutine hfecQ_tbl0()
  end subroutine hfecQ_tbl0

end module hfecQ_tbl
