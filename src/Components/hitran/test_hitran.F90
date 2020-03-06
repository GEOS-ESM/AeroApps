! Simple program to test the HITRAN module
! ----------------------------------------
program test_hitran

  use HITRAN

  implicit none

  character(len=256)         :: hitran_path
  character(len=6)           :: the_molecule
  real*8                     :: fwhm
  integer                    :: nlambda
  real*8,allocatable         :: lambda(:)
  real*8,allocatable         :: cross_section(:)  ![1/(molecule cm-2 )]
  logical                    :: is_wavenum
  integer                    :: nz
  real*8,allocatable         :: ps(:)  ! atm
  real*8,allocatable         :: ts(:)  ! K
  logical                    :: do_dP, do_dT
  integer                    :: errstat
  character(len=256)         :: errmessage


  hitran_path  = './'
  the_molecule = 'NO2'
  fwhm         = 0.6 ! nm
  nlambda      = 2
  allocate( lambda(nlambda) )
  allocate( cross_section(nlambda))
  lambda(1)    = 200 ! nm
  lambda(2)    = 405 ! nm
  is_wavenum   = .FALSE.
  nz           = 1
  allocate( ps(nz) )
  allocate( ts(nz) )
  ps(1) = 1.0
  ts(1) = 298
  do_dP = .FALSE.
  do_dT = .FALSE.
  
  write(*,*) 'inputs set',trim(hitran_path), the_molecule, nlambda, lambda, &
           is_wavenum, nz, ps, ts, fwhm, do_dP, do_dT
  call get_HITRAN_cross_section &
         ( hitran_path, the_molecule, nlambda, lambda, &
           is_wavenum, nz, ps, ts, fwhm, do_dP, do_dT, &
           cross_section, errstat, errmessage )  

  write(*,*) the_molecule,' cross section',cross_section
  write(*,*) 'STATUS',errstat,errmessage
 
end program test_HITRAN
