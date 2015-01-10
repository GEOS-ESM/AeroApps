!
!  Simple f77 wrapper for the Python interface to VLIDORT for OMI aerosol
!  channels.
!
!.............................................................................

subroutine Scalar (km, nch, nobs,channels, nMom_cd,  &
                   tau, ssa, g,cod,ssac,cdsmom,      &
                   pe, he, te, albedo,               &
                   solar_zenith, relat_azymuth, sensor_zenith, &
                   MISSING,verbose,radiance_VL,reflectance_VL, rc)
!
! Uses VLIDORT in scalar mode to compute OMI aerosol TOA radiances.
!
  use VLIDORT_ScatMod_cd

  implicit NONE

! !INPUT PARAMETERS:

  integer,          intent(in)  :: km    ! number of levels on file
  integer,          intent(in)  :: nch   ! number of channels
  integer,          intent(in)  :: nobs  ! number of observations
                                                      
  real*8, target,   intent(in)  :: channels(nch)    ! wavelengths [nm]
  integer, target,  intent(in)  :: nMom_cd          ! number of moments P11 cloud 

!                                                   ! --- Mie Parameters ---
  real*8, target,   intent(in)  :: tau(km,nch,nobs) ! aerosol optical depth
  real*8, target,   intent(in)  :: ssa(km,nch,nobs) ! single scattering albedo
  real*8, target,   intent(in)  :: g(km,nch,nobs)   ! asymmetry factor

  real*8, target,   intent(in)  :: cod(km,nch,nobs)    ! cloud optical depth
  real*8, target,   intent(in)  :: ssac(km,nch,nobs)   ! cloud single scattering albedo
  real*8, target,   intent(in)  :: cdsmom(nMom_cd) ! cloud phase function P11

  real*8, target,   intent(in)  :: pe(km+1,nobs)    ! pressure at layer edges [Pa]
  real*8, target,   intent(in)  :: he(km+1,nobs)    ! height above sea-level  [m]
  real*8, target,   intent(in)  :: te(km+1,nobs)    ! temperature at layer edges [K]
  
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
!                               ---  
  integer             :: i,j, ier

  real*8              :: AOT(nobs,nch)            ! total aerosol optical thickness 
  real*8              :: Q(nobs,nch),U(nobs,nch)    
  type(VLIDORT_scat) :: SCAT

#define IS_MISSING(x) (abs(x/MISSING-1)<0.001)
  rc = 0
  ier = 0

  call VLIDORT_Init( SCAT%Surface%Base, km, rc)
  if ( rc /= 0 ) return


  SCAT%cdsmom => cdsmom(:)
  SCAT%nMom_cd = nMom_cd
  
  print*, 'nmom,npol,nmom_cd', nMom_cd

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
           
           call VLIDORT_SurfaceLamb(SCAT%Surface,albedo(j,i),solar_zenith (j),sensor_zenith(j),&
                               relat_azymuth(j),.true.)
           
           SCAT%wavelength = channels(i)
           SCAT%tau => tau(:,i,j)
           SCAT%ssa => ssa(:,i,j)
           SCAT%g => g(:,i,j)
           SCAT%cod => cod(:,i,j)
           SCAT%ssac => ssac(:,i,j)
                print*, 'nmom,npol,nmom_cd', nMom_cd
           
           call VLIDORT_Run (SCAT, radiance_VL(j,i), reflectance_VL(j,i), &
                                 AOT(j,i),ier,Q(j,i),U(j,i), .true., .true.)

           if ( ier /= 0 ) then
              if ( solar_zenith (j) > 89.9 ) then ! if SZA > 90
              radiance_VL(j,i) = 0.0
              reflectance_VL(j,i) = 0.0
              else
              radiance_VL(j,i) = MISSING
              reflectance_VL(j,i) = MISSING
              endif
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


