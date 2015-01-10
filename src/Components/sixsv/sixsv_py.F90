  subroutine ocnAlbedo (speed, wavelength, del_azim, salinity, pigment, nobs, albedo)
  implicit NONE

  integer, intent(in) :: nobs
  real*8, intent(in) :: speed(nobs)       ! speed of wind (in m/s)
  real*8, intent(in) :: del_azim(nobs)    ! azim. of sun - azim. of wind (in deg.)
  real*8, intent(in) :: salinity(nobs)    ! salinity (in ppt) 
  real*8, intent(in) :: pigment(nobs)     ! pigment concentration (in mg.m-3)
  real*8, intent(in) :: wavelength        ! wavelength of the computation (in nanometer)

  real*8, intent(out) :: albedo(nobs)     ! the spherical albedo of the ocean 

!                            ----

  integer :: i
  real :: pws,paw,xsal,pcl,pwl,brdfalbe

  do i = 1, nobs
     pws = speed(i)
     paw = del_azim(i)
     xsal = salinity(i)
     pcl = pigment(i)
     pwl = wavelength/1000.
     call oceaalbe(pws,paw,xsal,pcl,pwl,brdfalbe)
     albedo(i) = brdfalbe
  end do

end subroutine ocnAlbedo
