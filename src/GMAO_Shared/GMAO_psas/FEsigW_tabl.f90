module FEsigW_tabl
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!-----------------------------------------------------------------------
!
! !ROUTINE: FEsigW_tabl - F90 module of FEsigW_tabl type definition
!
! !INTERFACE: (to do)
! !DESCRIPTION: (to do)
! !EXAMPLES: (to do)
! !BUGS: (to do)
! !SEE ALSO: (to do)
! !SYSTEM ROUTINES: (to do)
!
! !REVISION HISTORY:
!	08Oct99	- Jing Guo
!		. Added _defined and _setbyRC to support m_sigDCWindErr
!		. Defined FEsigW_Mlev, hoping to remove the dependency
!		  on lvmax and MXpar_hc (config.mod) in the future.
! 	23Apr96 - J. Guo	- (to do)
!	19Sep96 - J. Guo	- add `save' to all public variables
!_______________________________________________________________________
use config, only : lvmax,MXpar => MXpar_hc
implicit none
private	! except
public	:: FEsigW_rsrc, FEsigW_type, FEsigW_desc
public	:: FEsigW_Mlev, FEsigW_nlev, FEsigW_plev
public	:: FEsigW_Mpar, FEsigW_npar, FEsigW_pars
public  :: FEsigW_defined
public  :: FEsigW_setbyRC

  character(len=*), parameter	:: FEsigW_rsrc = 'FcstErr*Sigma_Wind::'

  logical,save :: FEsigW_defined=.false.
  logical,save :: FEsigW_setbyRC=.false.

  character(len= 32),save :: FEsigW_type	! a class name
  character(len=128),save :: FEsigW_desc	! a description

  integer,      save :: FEsigW_nlev		! number of levels
  real,         save :: FEsigW_plev(lvmax)	! levels defined

  integer,      save :: FEsigW_npar		! number of parameters
  real,         save :: FEsigW_pars(MXpar,lvmax)! parameters defined

  integer, parameter :: FEsigW_Mpar = MXpar
  integer, parameter :: FEsigW_Mlev = lvmax

end module FEsigW_tabl
!.
