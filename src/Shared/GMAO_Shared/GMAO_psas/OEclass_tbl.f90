!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!-----------------------------------------------------------------------
!
! !ROUTINE: OEclass_tbl - table of ObsErr classes
!
! !INTERFACE:

! !DESCRIPTION: (to do)
! !EXAMPLES: (to do)
! !BUGS: (to do)
! !SEE ALSO: (to do)
! !SYSTEM ROUTINES: (to do)
!
! !REVISION HISTORY:
!	19Sep96 - J. Guo	- add `save' to all public variables
!_______________________________________________________________________

module OEclass_tbl
use config, only : lvmax_oe,kxmax,ktmax
implicit none
public

	! The resource name of OE levels
  character*14, parameter :: RC_OEplev='ObsErr*Levels:'

  character*8,  save ::	levtype
  integer,      save ::	nlev_oe		! number of levels, <=lvmax
  real,         save :: plev_oe(lvmax_oe)! pressure levels (for small tables)

  integer,      save ::	nOEclas		! nOEclas <= kxmax
  character*8,  save ::	OEclas(kxmax)	! error classes
  logical,      save :: KTclas(ktmax,kxmax)	! permitted datatype(kt)

!-------------------------------
  character*16, save ::	hoecH_c(kxmax)	! correlation types of hoecH_c
  character*16, save ::	voecH_c(kxmax)	! correlation types of voecH_c

  integer, save :: loc_sigOc(ktmax,kxmax)
  integer, save :: len_sigOc(ktmax,kxmax)
  real,    save :: sigOc(lvmax_oe,ktmax,kxmax)

!-------------------------------
  character*16, save ::	voecH_u(kxmax)	! correlation types of voecH_u

  integer, save :: loc_sigOu(ktmax,kxmax)
  integer, save :: len_sigOu(ktmax,kxmax)
  real,    save :: sigOu(lvmax_oe,ktmax,kxmax)

!-------------------------------

Contains

  subroutine OEclass_tbl0()
  end subroutine OEclass_tbl0

end module OEclass_tbl
!.
