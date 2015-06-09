
!	----------------------
!	This is an include file to define some constants and variables 
! used in the geo_vlidort program and subroutines
!-------------------------

	! VLIDORT array sizes
	integer, parameter                    :: nobs = 1               ! number of profiles VLIDORT will work on
  real, parameter                       :: ptop = 1.0             ! top (edge) pressure [Pa]
  integer, parameter                    :: nq   = 14                ! number of tracers

  ! Vector calculation constants
  integer, parameter 										:: nPol = 6                ! number of components of the scattering matrix
  integer, parameter 										:: nMom = 300 						 ! number of phase function moments - 300 is max for dust

  ! Verbose
  integer, parameter                    :: verbose    = 0
  integer, parameter                    :: verbose_mp = 0

  ! MODIS kernel stuff
  integer, parameter                    :: nkernel = 3
  integer, parameter                    :: nparam  = 2
  integer, parameter                    :: nbands  = 7                ! number of modis bands

  ! Missing Value 
  real, parameter                       :: MISSING = 1e15

  ! RC file
  character(len=*), parameter           :: rcfile            = 'Aod_EOS.rc'  ! resource file
  character(len=16), parameter          :: vnames_string(nq) = (/'du001', 'du002', 'du003', 'du004', 'du005', &
                                                             'ss001', 'ss002', 'ss003', 'ss004', 'ss005', &
                                                             'BCphobic', 'BCphilic',                      &
                                                             'OCphobic', 'OCphilic'/) ! array of variable name strings
  character                             :: vnames(nq,16)          ! character array of variable names

  ! Physical constants
  real, parameter       								:: grav = 9.81               ! gravity

  ! Cloud Filtering
  real, parameter                     	:: cldmax = 0.01  

  ! Geomtry filtering
  real, parameter                       :: szamax = 90.0

