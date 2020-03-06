! Simple program to test the SOLAR module
! ----------------------------------------
program test_solar

  use SOLAR

  implicit none

  character(len=256)         :: solar_path
  integer                    :: nheader
  integer                    :: nwavenum
  real*8, allocatable        :: wavenum(:)
  real*8, allocatable        :: irradiance(:)
  real*8                     :: fwhm
  real*8                     :: clambda(1)
  real*8                     :: cwavenum(1)
  real*8                     :: cirradiance(1)

  nheader = 2
  solar_path = 'newkur.dat'

  nwavenum = get_SOLAR_nwavenum(solar_path,nheader)
  write(*,*) 'nwavenum',nwavenum

  allocate( wavenum(nwavenum) )
  allocate( irradiance(nwavenum) )
  call read_SOLAR_reference(solar_path,nheader,nwavenum,wavenum,irradiance)

  write(*,*) 'wavenum',wavenum(1)
  write(*,*) 'irradiance',irradiance(1)

  fwhm = 0.6 ! nm
  cwavenum(1) = 1.0D7/400
  call convolve_SOLAR_instrument(wavenum,irradiance,1.0D7/fwhm,cwavenum,cirradiance)
  write(*,*) '400 nm irradiance',cirradiance

 
end program test_solar