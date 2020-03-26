module config
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!-----------------------------------------------------------------------
!
! !ROUTINE: config - parameters shared by both F77 and F90 programs
!
! !INTERFACE:
!	use config
!
! !DESCRIPTION:
!	Module config defines constants normally defined in include
!	files (.h) in F77.  Future version may remove include, but
!	use only parameter settings, when all files are converted to
!	Fortran 90 and take the advantage of "use <module>" statement.
!
! !EXAMPLES:
!	use config, only : lvmax
!
! !REVISION HISTORY:
!	25Aug00	- Jing Guo
!		. Removed "rade.h" where a wrong Radius_of_Earth value
!		  has been used.
! 	25Jan96 - J. Guo	- prolog added
!_______________________________________________________________________
  include "lvmax.h"
  include "kxmax.h"
  include "ktmax.h"
  include "latmax.h"

  include "MX_hfecH.h"
  include "MX_hoecH.h"
  include "MXpar_hc.h"

  include "MX_voecH.h"

  integer, parameter	:: lvmax_vc=lvmax
  integer, parameter	:: lvmax_hc=lvmax
  integer, parameter	:: lvmax_oe=lvmax

  include "kind_covs.h"
  include "kind_mats.h"
  include "pres4slp.h"

  integer, parameter	:: MXveclat=5*180+1

	! Default GEOS "FILL" value
  real,    parameter    :: FILL=1.e+15
end module config
