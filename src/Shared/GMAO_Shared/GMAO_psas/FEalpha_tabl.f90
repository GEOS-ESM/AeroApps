module FEalpha_tabl
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!-----------------------------------------------------------------------
!
! !ROUTINE: FEalpha_tabl - a parameter tables of the ALPHA operator
!
! !INTERFACE:
!	use FEalpha_tabl
!
! !DESCRIPTION:
!	Alpha_tabl defines interpolation tables for ALPHA operators.
!	The tables are used to produce indexed matrix/transpose ALPHA
!	operators.  Tables are read from the resource file with a table
!	titled by RC_Aref for a table Aref_tabl(:) vs. rlev_A.
!
! !EXAMPLES: (to do)
! !BUGS: (to do)
!
! !SEE ALSO: tabl_FEalpha() and imat_alpha().
!
! !REVISION HISTORY:
! 	09Feb96 - J. Guo	- (to do)
!	13Feb97 - Jing Guo	- fixed the foggotten SAVEs
!_______________________________________________________________________

  use config, only : lvmax
  use config, only : MXpar => MXpar_hc
  implicit none

  private	! except:

  public :: FEalpha_rsrc	! Name of the parameter table
  public :: FEalpha_type	! Model the parameter table
  public :: FEalpha_desc	! A one-line description of the model
  public :: FEalpha_nlev	! number of actual table levels
  public :: FEalpha_plev	! pressure levels
  public :: FEalpha_npar	! number of parameters for a level
  public :: FEalpha_pars	! the parameter table
  public :: FEalpha_Mpar
  public :: FEalpha_tabl0

	! The name of the table

  character(len=*), parameter	:: FEalpha_rsrc    = "FcstErr*Aref::"

  character(len=32), save :: FEalpha_type ! model of the parameter table
  character(len=128),save :: FEalpha_desc ! description of the model

  integer,save		:: FEalpha_nlev	! number of levels
  real,   save		:: FEalpha_plev(lvmax)	! a level table


  integer,save		:: FEalpha_npar	! actual number of parameters
  real,   save		:: FEalpha_pars(MXpar,lvmax)	! parameters

  integer, parameter	:: FEalpha_Mpar = MXpar

Contains

  subroutine FEalpha_tabl0
  end subroutine FEalpha_tabl0

end module FEalpha_tabl
!.
