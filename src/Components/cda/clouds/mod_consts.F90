! file: mod_consts.f90
! physical constants

module mod_consts

  implicit none

  ! dont change, from m_const.f90
  real*4, parameter :: cpm   = 1004.64
  real*4, parameter :: rgas  = 287.04
  real*4, parameter :: kappa = rgas/cpm
  real*4, parameter :: zvir  = 4.61e2/rgas - 1.  ! should be rvap/rgas-1
  real*4, parameter :: eps   = 0.622             ! should be rgas/rvap
  real*4, parameter :: grav  = 9.80616
  real*4, parameter :: alhl  = 2.499e6  ! Latent heat consensation
  real*4, parameter :: alhs  = 2.845e6  ! Latent heat sublimation

  ! derived
  real*4, parameter :: reps = 1. / eps
  real*4, parameter :: omreps = 1. - reps
  real*4, parameter :: UNDEF = 1.0e15

  ! from physparams
  real*4, parameter :: MIN_RL = 10.e-6
  real*4, parameter :: MAX_RL = 21.e-6
  real*4, parameter :: MIN_RI = 20.e-6
  real*4, parameter :: MAX_RI = 40.e-6
  real*4, parameter :: RI_ANV = 30.e-6

end module mod_consts
