!
! Compilation note:
!
! With gfortran 4.1 or older, specify -frecord-marker=4 
!
!     gfortran -frecord-marker=4 g5pace_py.F90
!
! With ifort, simply
!
!     ifort g5pace_py.F90
!
! If using python, you can build an extension with:
!
!     f2py -c -m g5pace_ -I. g5pace_py.F90
!
! and then
!
! >>> import g5pace_ as g
! >>> g.read('g5aero_p01.bin')
! >>> print g.g5pace.tau
!
!...........................................................................

module G5PACE

   integer(4), parameter :: Version=1      ! File format version number
   integer(4), parameter :: nLayers=72     ! number of vertical layers
   integer(4), parameter :: nScat=83       ! number of scattering angles
   integer(4), parameter :: nChannels=26   ! number of channels

   integer(4) :: year, month, day, hour, min, sec ! date & time (UTC)
   real(4) :: lon, lat        ! longitude, latitude [degrees]

   real(4) :: u10m, v10m      ! 10 meter wind components

   real(4) :: channels(nChannels) ! wavelength in nm
   real(4) :: scat_angles(nScat)  ! scattering angles

   real(4) :: pe(nLayers+1)   ! pressure at layer edges [Pa]
   real(4) :: he(nLayers+1)   ! height at layer edges [m]
   real(4) :: te(nLayers+1)   ! temperature at layer edges [K]

   real(4) :: qv(nLayers)     ! water vapor specific humidity [kg/kg]
   real(4) :: o3(nLayers)     ! ozone mass mixing ratio [kg/kg]

   real(4) :: tau(nLayers,nChannels)  ! AOD of layer
   real(4) :: ssa(nLayers,nChannels)  ! aerosol single scattering albedo
   real(4) :: g(nLayers,nChannels)    ! asymmetry parameter

                                   ! Phase matrix
   real(4) :: P11(nLayers,nScat,nChannels) 
   real(4) :: P12(nLayers,nScat,nChannels)
   real(4) :: P33(nLayers,nScat,nChannels)
   real(4) :: P34(nLayers,nScat,nChannels)

end module G5PACE

!---------------------------
subroutine write_(filename)

   use G5PACE
   implicit NONE

   character(len=*), intent(in) :: filename

   print *
   print *, '[] Writing '//trim(filename)
   call print_()

   open(10,file=filename,access='sequential',form='unformatted')
   write(10) Version, nLayers, nScat, nChannels
   write(10) year, month, day, hour, min, sec, &
             lon, lat, channels, scat_angles, &
	     u10m, v10m, qv, o3, &
             pe, he, te,  &
             tau, ssa, g, P11, P12, P33, P34
   close(10)

end subroutine write_

!------------------
subroutine print_()
  use G5PACE
  implicit NONE

  print *, '   -       When: ', year, month, day, hour, min, sec
  print *, '   -       Where: ', lon, lat
  print *, '   -        Wind: ', u10m, v10m
  print *, '   -    Channels: ', channels
  print *, '   - Scat Angles: ', scat_angles
  print *
  print *, '      pe = ', minval(pe),  maxval(pe)
  print *, '      he = ', minval(he),  maxval(he)
  print *, '      te = ', minval(te),  maxval(te)
  print *, '      qv = ', minval(qv),  maxval(qv)
  print *, '      o3 = ', minval(o3),  maxval(o3)
  print *, '     tau = ', minval(tau), maxval(tau)
  print *, '     ssa = ', minval(ssa), maxval(ssa)
  print *, '       g = ', minval(g),   maxval(g)
  print *, '     p11 = ', minval(p11), maxval(p11)
  print *, '     p12 = ', minval(p12), maxval(p12)
  print *, '     p33 = ', minval(p33), maxval(p33)
  print *, '     p34 = ', minval(p34), maxval(p34)

end subroutine print_

!-------------------------
subroutine read_(filename)

   use G5PACE
   implicit NONE

   character(len=*), intent(in) :: filename

   integer(4) :: Version_, nLayers_, nScat_, nChannels_

   print *
   print *, '[] Reading '//trim(filename)

   open(10,file=filename,access='sequential',form='unformatted')
   read(10) Version_, nLayers_, nScat_, nChannels_

   if ( Version    /= Version_ .OR. &
        nLayers_   /= nLayers  .OR. &
        nScat_     /= nScat    .OR. &
        nChannels_ /= nChannels ) then
      print *, '*** ERRROR ***'
      print *, 'Expecting Version, nLayers, nScat, nChannels = ', Version, nLayers, nScat, nChannels
      print *, 'But found Version, nLayers, nScat, nChannels = ', Version_, nLayers_, nScat_, nChannels_
      call exit(1)
   end if

   read(10) year, month, day, hour, min, sec, &
            lon, lat, channels, scat_angles, &
	    u10m, v10m, qv, o3, &
            pe, he, te,  &
            tau, ssa, g, P11, P12, P33, P34

   close(10)

   call print_()

end subroutine read_

!...............................
program Reader

 call read_('g5aero_p01.bin')
 call read_('g5aero_p02.bin')

end program Reader
