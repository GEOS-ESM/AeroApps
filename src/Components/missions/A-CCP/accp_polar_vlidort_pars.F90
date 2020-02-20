
! ----------------------
! This is an include file to define some constants and variables 
! used in the geo_vlidort program and subroutines
!-------------------------

  ! VLIDORT array sizes
  integer, parameter                    :: nobs = 1               ! number of profiles VLIDORT will work on
  real, parameter                       :: ptop = 1.0             ! top (edge) pressure [Pa]
  integer, parameter                    :: nq   = 15                ! number of tracers

  ! Vector calculation constants
  integer, parameter                    :: nPol = 6                ! number of components of the scattering matrix
  integer, parameter                    :: nMom = 300              ! number of phase function moments - 300 is max for dust

  ! Verbose
  integer, parameter                    :: verbose    = 0
  integer, parameter                    :: verbose_mp = 0

  ! MODIS kernel stuff
  integer, parameter                    :: nkernel = 3
  integer, parameter                    :: nparam  = 2

  ! Missing Value 
  real, parameter                       :: MISSING = 1e15

  ! RC file
  character(len=16), parameter          :: vnames_string(nq) = (/'du001', 'du002', 'du003', 'du004', 'du005', &
                                                             'ss001', 'ss002', 'ss003', 'ss004', 'ss005', &
                                                             'BCphobic', 'BCphilic',                      &
                                                             'OCphobic', 'OCphilic','SO4'/) ! array of variable name strings
  character(len=16), parameter          :: vnames_stringDU(5) = (/'du001', 'du002', 'du003', 'du004', 'du005'/)
  character(len=16), parameter          :: vnames_stringSS(5) = (/'ss001', 'ss002', 'ss003', 'ss004', 'ss005'/)
  character(len=16), parameter          :: vnames_stringBC(2) = (/'BCphobic', 'BCphilic'/)
  character(len=16), parameter          :: vnames_stringOC(2) = (/'OCphobic', 'OCphilic'/) 
  character(len=16), parameter          :: vnames_stringSU(1) = (/'SO4'/)   
  character                             :: vnames(nq,16)          ! character array of variable names
  character                             :: vnames_du(5,16)          ! character array of variable names  
  character                             :: vnames_ss(5,16)          ! character array of variable names  
  character                             :: vnames_bc(2,16)          ! character array of variable names  
  character                             :: vnames_oc(2,16)          ! character array of variable names 
  character                             :: vnames_su(1,16)          ! character array of variable names    


  ! Physical constants
  real, parameter                       :: grav = 9.81               ! gravity



