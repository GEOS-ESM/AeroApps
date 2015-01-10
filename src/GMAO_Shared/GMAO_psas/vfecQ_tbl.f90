module vfecQ_tbl
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!-----------------------------------------------------------------------
!
! !ROUTINE: vfecQ_tbl - a F90 module of vfecQ input tables
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
use config, only : lvmax_vc
implicit none
public

  character(len=*),parameter :: name_vfecQ='FcstErr*vCor_QQ::'

  character*16,save :: type_vfecQ
  character*64,save :: desc_vfecQ

  integer,save :: nlev_vfecQ
  real,   save :: plev_vfecQ(lvmax_vc)
  real,   save :: corr_vfecQ(lvmax_vc,lvmax_vc)

Contains

  subroutine vfecQ_tbl0()
  end subroutine vfecQ_tbl0

end module vfecQ_tbl
