   subroutine ThreshVelocity (n, uthresh0, radius_microns)

! !USES:

  implicit NONE

! !INPUT PARAMETERS:

   real,    intent(in) :: radius_microns(n)        ! particle radius [micron]
   integer, intent(in) :: n

! !OUTPUT PARAMETERS:

   real,    intent(out) :: uthresh0(n)             ! Local emission

!EOP
!-------------------------------------------------------------------------

! !Local Variables
   integer :: i
   real, parameter ::  air_dens = 1.25       ! Air density = 1.25 kg m-3
   real, parameter ::  soil_density = 2650. ! km m-3
   real, parameter ::  grav = 9.80           ! this is what GEOS-5 uses, weird!
   real            ::  diameter              ! dust effective diameter [m]
 
!  Calculate the threshold velocity of wind erosion [m/s] for each radius
!  for a dry soil, as in Marticorena et al. [1997].
!  The parameterization includes the air density which is assumed 
!  = 1.25 kg m-3 to speed the calculation.  The error in air density is
!  small compared to errors in other parameters.
!  ----------------------------------------------------------------------
   do i = 1, n
      diameter = 2. * radius_microns(i) * 1e-6 ! in meters
      uthresh0(i) = 0.13 * sqrt(soil_density*grav*diameter/air_dens) &
                    * sqrt(1.+6.e-7/(soil_density*grav*diameter**2.5)) &
           / sqrt(1.928*(1331.*(100.*diameter)**1.56+0.38)**0.092 - 1.)
   end do

   return 
   end subroutine ThreshVelocity

!
! !IROUTINE:  DustEmissionGOCART - Compute the dust emissions
!
! !INTERFACE:
!

   subroutine DustEmissionGOCART( nobs, radius, fracland, gwettop, w10m, &
                                  emissions )

! !USES:

  implicit NONE

! !INPUT PARAMETERS:

   integer, intent(in) :: nobs
   real,    intent(in) :: radius            ! particle radius [m]
   real,    intent(in) :: fracland(nobs)    ! land fraction
   real,    intent(in) :: gwettop(nobs)     ! ground wetness
   real,    intent(in) :: w10m(nobs)        ! 10m wind speed

! !OUTPUT PARAMETERS:

   real,    intent(out) :: emissions(nobs)         ! Local emission

! !DESCRIPTION: Computes the dust emissions dependence on wind speed
!               and soil wetness for one time step.
!               GOCART-based dust emission scheme (modified from
! Ginoux et al., JGR, 2001) .Pass in a grid of wind speeds, surface
! characteristics, and an aerosol particle radius and returns the
! total dust emission mass flux [kg m-2 s-1].  Emissions need to be
! scaled by source function (ie., fraction of grid cell emitting),
! tuning factor, and fractional soil content of particle size.
!
!
! !REVISION HISTORY:
!
!  24Jan2016  da Silva - Based on code from Chem_Shared. 
!  29Dec2009  Colarco  - Modifications to change calling
!  06Nov2003  Colarco  - Based on Ginoux
!
!EOP
!-------------------------------------------------------------------------

! !Local Variables
   real, parameter ::  air_dens = 1.25       ! Air density = 1.25 kg m-3
   real, parameter ::  soil_density  = 2650. ! km m-3
   real, parameter ::  grav = 9.80           ! this is what GEOS-5 uses, weird!
   integer         ::  i
   real            ::  diameter              ! dust effective diameter [m]
   real            ::  u_thresh0
   real            ::  u_thresh
 
!  Initialize local variables
!  --------------------------
   emissions(:) = 0.

!  Calculate the threshold velocity of wind erosion [m/s] for each radius
!  for a dry soil, as in Marticorena et al. [1997].
!  The parameterization includes the air density which is assumed 
!  = 1.25 kg m-3 to speed the calculation.  The error in air density is
!  small compared to errors in other parameters.
!  ----------------------------------------------------------------------
   diameter = 2. * radius
   u_thresh0 = 0.13 * sqrt(soil_density*grav*diameter/air_dens) &
                    * sqrt(1.+6.e-7/(soil_density*grav*diameter**2.5)) &
           / sqrt(1.928*(1331.*(100.*diameter)**1.56+0.38)**0.092 - 1.)
   
   print *, 'u_thresh0 = ',  u_thresh0 

!  Spatially dependent part of calculation
!  ---------------------------------------
   do i = 1, nobs

!    Modify the threshold depending on soil moisture 
!    as in Ginoux et al. [2001]
!    -----------------------------------------------
     if(gwettop(i) .lt. 0.5) then

       u_thresh = amax1(0.,u_thresh0* &
                 (1.2+0.2*alog10(max(1.e-3,gwettop(i)))))

!      Emission of dust [kg m-2 s-1]
!      -----------------------------
       if(w10m(i) .gt. u_thresh) then     
          emissions(i) = fracland(i) * w10m(i)**2. * (w10m(i)-u_thresh)
       end if

      end if

    end do ! i

  end subroutine DustEmissionGOCART

