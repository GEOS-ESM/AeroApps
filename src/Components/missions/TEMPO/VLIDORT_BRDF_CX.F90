!
!  Simple f77 wrapper for the Python interface to VLIDORT for OMI aerosol
!  channels.
!
!.............................................................................

subroutine Scalar (km, nch, nobs,channels,        &
                   tau, ssa, g, pe, he, te, albedo, U10m,V10m, &
                   mr, solar_zenith, relat_azymuth, sensor_zenith, &
                   MISSING,verbose,radiance_VL,reflectance_VL,reflectance_cx, BRDF,rc)
!
! Uses VLIDORT in scalar mode to compute OMI aerosol TOA radiances.
!
  use VLIDORT_ScatMod

  implicit NONE

! !INPUT PARAMETERS:

  integer,          intent(in)  :: km    ! number of levels on file
  integer,          intent(in)  :: nch   ! number of channels
  integer,          intent(in)  :: nobs  ! number of observations
                                      
                  
  real*8, target,   intent(in)  :: channels(nch)    ! wavelengths [nm]

!                                                   ! --- Mie Parameters ---
  real*8, target,   intent(in)  :: tau(km,nch,nobs) ! aerosol optical depth
  real*8, target,   intent(in)  :: ssa(km,nch,nobs) ! single scattering albedo
  real*8, target,   intent(in)  :: g(km,nch,nobs)   ! asymmetry factor


  real*8, target,   intent(in)  :: pe(km+1,nobs)    ! pressure at layer edges [Pa]
  real*8, target,   intent(in)  :: he(km+1,nobs)    ! height above sea-level  [m]
  real*8, target,   intent(in)  :: te(km+1,nobs)    ! temperature at layer edges [K]

  real*8, target,   intent(in)  :: U10m(nobs)   ! Wind speed components [m/s]
  real*8, target,   intent(in)  :: V10m(nobs)    
  real*8, target,   intent(in)  :: mr(nch)       ! refractive index
  
  real*8, target,   intent(in)  :: MISSING          ! MISSING VALUE                                      
  real*8, target,   intent(in)  :: solar_zenith(nobs)  
  real*8, target,   intent(in)  :: relat_azymuth(nobs) 
  real*8, target,   intent(in)  :: sensor_zenith(nobs) 

  real*8, target,   intent(in)  :: albedo(nobs,nch)       ! surface albedo
  
  integer,          intent(in)  :: verbose

! !OUTPUT PARAMETERS:

  real*8,           intent(out) :: radiance_VL(nobs,nch)       ! TOA normalized radiance from VLIDORT
  integer,          intent(out) :: rc                          ! return code
  real*8,           intent(out) :: reflectance_VL(nobs, nch)   ! TOA reflectance from VLIDORT
  real*8,           intent(out) :: reflectance_cx(nobs, nch)   ! TOA reflectance from VLIDORT Cox Munk
  real*8,           intent(out) :: BRDF(nobs, nch)  
!                         ---  
  integer             :: i,j, ier

  real*8              :: ROT(km,nobs,nch)            ! rayleigh optical thickness 
  real*8              :: Q(nobs, nch)   ! Stokes parameter Q
  real*8              :: U(nobs, nch)   ! Stokes parameter U     
  type(VLIDORT_scat) :: SCAT

#define IS_MISSING(x) (abs(x/MISSING-1)<0.001)
  rc = 0
  ier = 0

  call VLIDORT_Init( SCAT%Surface%Base, km, rc)
  if ( rc /= 0 ) return

  do j = 1,nobs

     ! Make sure albedo and angles are available
     ! -----------------------------------------
     if ( IS_MISSING(solar_zenith(j))  .OR. & 
          IS_MISSING(sensor_zenith(j)) .OR. &
          IS_MISSING(relat_azymuth(j))  )  then

        radiance_VL(j,:) = MISSING
        reflectance_VL(j,:) = MISSING
        cycle

      end if
      
      SCAT%pe => pe(:,j)
      SCAT%ze => he(:,j)
      SCAT%te => te(:,j) 
 
      ! Loop over channels
      ! ------------------
      do i = 1, nch 
       
           ! Mare sure albedo is defined
           ! ---------------------------
           if ( IS_MISSING(albedo(j,i)) ) then
              radiance_VL(j,i) = MISSING
              reflectance_VL(j,i) = MISSING
              cycle
           end if
           print*, 'DO SURFACE LAMB'
           call VLIDORT_SurfaceLamb(SCAT%Surface,albedo(j,i),solar_zenith (j),sensor_zenith(j),&
                               relat_azymuth(j),.true.)

           SCAT%wavelength = channels(i)
           SCAT%tau => tau(:,i,j)
           SCAT%ssa => ssa(:,i,j)
           SCAT%g => g(:,i,j)
         
           call VLIDORT_Run (SCAT, radiance_VL(j,i), reflectance_VL(j,i), &
                                 ROT(:,j,i),Q(j,i),U(j,i), .true., .true., ier)
!           print *, 'radiance albedo',albedo(j,i),radiance_VL(j,i), reflectance_VL(j,i) 
           if ( ier /= 0 ) then
              radiance_VL(j,i) = MISSING
              reflectance_VL(j,i) = MISSING
              cycle
           end if
           print*, 'DO COX MUNK'
           call VLIDORT_GissCoxMunk(SCAT%Surface,U10m(j),V10m(j),mr(i),solar_zenith (j),&
                                    sensor_zenith(j),relat_azymuth(j),.true.,BRDF(j,i),rc)
           if ( rc /= 0 ) return

           SCAT%wavelength = channels(i)
           SCAT%tau => tau(:,i,j)
           SCAT%ssa => ssa(:,i,j)
           SCAT%g => g(:,i,j)
         
           call VLIDORT_Run (SCAT, radiance_VL(j,i), reflectance_cx(j,i), &
                                 ROT(:,j,i),Q(j,i),U(j,i), .true., .true., ier)

!           print *, 'radinace cox munk',radiance_VL(j,i), reflectance_cx(j,i) 
           if ( ier /= 0 ) then
              radiance_VL(j,i) = MISSING
              reflectance_VL(j,i) = MISSING
              cycle
           end if


        end do ! end loop over channels
     
        if ( verbose > 0 ) then
           if ( mod(j-1,1000) == 0 ) then
              print *, '<> VLIDORT Scalar: ', nint(j*100./nobs), '%'
           end if
        end if

  end do ! Loop over obs

end subroutine Scalar

!..........................................................................

subroutine Vector (km, nch, nobs, channels, nMom,  &
                   nPol,tau, ssa, g, pmom, pe, he, te, albedo, &
                   U10m, V10m, mr, &
                   solar_zenith, relat_azymuth, sensor_zenith, &
                   MISSING,verbose, radiance_VL,reflectance_VL,reflectance_cx, BRDF,rc)
!
! Place holder.
!
   use VLIDORT_ScatMod
 
   implicit NONE

! !INPUT PARAMETERS:

  integer,          intent(in)  :: km    ! number of levels on file
  integer,          intent(in)  :: nch   ! number of channels
  integer,          intent(in)  :: nobs  ! number of observations

  integer, target,  intent(in)  :: nMom  ! number of moments 
  integer, target,  intent(in)  :: nPol  ! number of components                               
                  
  real*8, target,   intent(in)  :: channels(nch)    ! wavelengths [nm]

!                                                   ! --- Mie Parameters ---
  real*8, target,   intent(in)  :: tau(km,nch,nobs) ! aerosol optical depth
  real*8, target,   intent(in)  :: ssa(km,nch,nobs) ! single scattering albedo
  real*8, target,   intent(in)  :: g(km,nch,nobs)   ! asymmetry factor
  
  real*8, target,   intent(in)  :: pmom(km,nMom,nPol,nch,nobs) !components of the scat phase matrix

  real*8, target,   intent(in)  :: MISSING          ! MISSING VALUE
  real*8, target,   intent(in)  :: pe(km+1,nobs)    ! pressure at layer edges [Pa]
  real*8, target,   intent(in)  :: he(km+1,nobs)    ! height above sea-level  [m]
  real*8, target,   intent(in)  :: te(km+1,nobs)    ! temperature at layer edges [K]

  real*8, target,   intent(in)  :: U10m(nobs)   ! Wind speed components [m/s]
  real*8, target,   intent(in)  :: V10m(nobs)    
  real*8, target,   intent(in)  :: mr(nch)       ! refractive index
                       
  real*8, target,   intent(in)  :: solar_zenith(nobs)  
  real*8, target,   intent(in)  :: relat_azymuth(nobs) 
  real*8, target,   intent(in)  :: sensor_zenith(nobs) 

  real*8, target,   intent(in)  :: albedo(nobs,nch)       ! surface albedo

  integer,          intent(in)  :: verbose

! !OUTPUT PARAMETERS:

   real*8,           intent(out) :: radiance_VL(nobs,nch)       ! TOA normalized radiance from VLIDORT
  integer,          intent(out) :: rc                          ! return code
  real*8,           intent(out) :: reflectance_VL(nobs, nch)   ! TOA reflectance from VLIDORT
  real*8,           intent(out) :: reflectance_cx(nobs, nch)   ! TOA reflectance from VLIDORT Cox Munk
  real*8,           intent(out) :: BRDF(nobs, nch)  
  real*8            :: Q(nobs, nch)   ! Stokes parameter Q
  real*8            :: U(nobs, nch)   ! Stokes parameter U
!                               ---
  
  integer             :: i,j, ier 
  real*8              :: ROT(km,nobs,nch)            ! rayleigh optical thickness 
  
  type(VLIDORT_scat) :: SCAT

#define IS_MISSING(x) (abs(x/MISSING-1)<0.001)
  rc = 0
  ier = 0
 
  call VLIDORT_Init( SCAT%Surface%Base, km, rc)
  if ( rc /= 0 ) return

  SCAT%nMom = nMom
  SCAT%nPol = nPol
  do j = 1, nobs
     
     ! Make sure albedo and angles are available
     ! -----------------------------------------
     if ( IS_MISSING(solar_zenith(j))  .OR. & 
          IS_MISSING(sensor_zenith(j)) .OR. &
          IS_MISSING(relat_azymuth(j))  )  then

        radiance_VL(j,:) = MISSING
        reflectance_VL(j,:) = MISSING
      
        cycle

      end if
     
     SCAT%pe => pe(:,j)
     SCAT%ze => he(:,j)
     SCAT%te => te(:,j) 

     do i = 1, nch
       if ( IS_MISSING(albedo(j,i)) ) then
              radiance_VL(j,i) = MISSING
              reflectance_VL(j,i) = MISSING
            
              cycle
       end if

       call VLIDORT_SurfaceLamb(SCAT%Surface,albedo(j,i),solar_zenith (j),sensor_zenith(j),&
                               relat_azymuth(j),.false.)

        SCAT%wavelength = channels(i)        
        SCAT%tau => tau(:,i,j)
        SCAT%ssa => ssa(:,i,j)
        SCAT%g => g(:,i,j)
        SCAT%pmom => pmom(:,:,:,i,j)

        call VLIDORT_Run (SCAT, radiance_VL(j,i),reflectance_VL(j,i),&
                          ROT(:,j,i), Q(j,i),U(j,i),.false., .true., ier)
        if ( ier /= 0 ) then
              radiance_VL(j,i) = MISSING
              reflectance_VL(j,i) = MISSING
             
              cycle
           end if
!        print*, 'DO COX MUNK'
           call VLIDORT_GissCoxMunk(SCAT%Surface,U10m(j),V10m(j),mr(i),solar_zenith (j),&
                                    sensor_zenith(j),relat_azymuth(j),.false.,BRDF(j,i),rc)
           if ( rc /= 0 ) return

           SCAT%wavelength = channels(i)
           SCAT%tau => tau(:,i,j)
           SCAT%ssa => ssa(:,i,j)
           SCAT%g => g(:,i,j)
         
           call VLIDORT_Run (SCAT, radiance_VL(j,i), reflectance_cx(j,i), &
                                 ROT(:,j,i) ,Q(j,i),U(j,i), .false., .true., ier)

!           print *, 'radinace cox munk',radiance_VL(j,i), reflectance_cx(j,i) 
           if ( ier /= 0 ) then
              radiance_VL(j,i) = MISSING
              reflectance_VL(j,i) = MISSING
             
              cycle
           end if
        end do ! end loop over channels
     
        if ( verbose > 0 ) then
           if ( mod(j-1,1000) == 0 ) then
              print *, '<> VLIDORT Vector: ', nint(j*100./nobs), '%'
           end if
        end if

  end do ! Loop over obs

end subroutine Vector

