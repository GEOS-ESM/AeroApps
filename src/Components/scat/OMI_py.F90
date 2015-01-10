!
!  Simple f77 wrapper for the Python interface to VLIDORT for OMI aerosol
!  channels.
!
!.............................................................................

subroutine Scalar (km, nch, nobs, nymd, nhms, channels,        &
                   tau, ssa, g, pe, he, te,                    & 
                   lon, lat, ps, albedo,                       &
                   solar_zenith, relat_azymuth, sensor_zenith, &
                   MISSING,radiance, AI, verbose,radiance_VL,          &
                   AI_VL, refl, rc)
!
! Uses VLIDORT in scalar mode to compute OMI aerosol TOA radiances.
!
  use VLIDORT_LambMod

  implicit NONE

! !INPUT PARAMETERS:

  integer,          intent(in)  :: km    ! number of levels on file
  integer,          intent(in)  :: nch   ! number of channels
  integer,          intent(in)  :: nobs  ! number of observations
  integer,          intent(in)  :: nymd  ! date, eg., 20080629
  integer,          intent(in)  :: nhms  ! time, eg., 120000
                                         
                  
  real*8, target,   intent(in)  :: channels(nch)    ! wavelengths [nm]

!                                                   ! --- Mie Parameters ---
  real*8, target,   intent(in)  :: tau(km,nch,nobs) ! aerosol optical depth
  real*8, target,   intent(in)  :: ssa(km,nch,nobs) ! single scattering albedo
  real*8, target,   intent(in)  :: g(km,nch,nobs)   ! asymmetry factor


  real*8, target,   intent(in)  :: pe(km+1,nobs)    ! pressure at layer edges [Pa]
  real*8, target,   intent(in)  :: he(km+1,nobs)    ! height above sea-level  [m]
  real*8, target,   intent(in)  :: te(km+1,nobs)    ! temperature at layer edges [K]

                                                    ! --- OMI Metadata ---
  real*8,           intent(in)  :: lon(nobs)        ! longitudes [-180,180]
  real*8,           intent(in)  :: lat(nobs)        ! latitudes [-90,90]
  real*8,           intent(in)  :: ps(nobs)         ! surface pressure [Pa]
                       
  real*8, target,   intent(in)  :: solar_zenith(nobs)  
  real*8, target,   intent(in)  :: relat_azymuth(nobs) 
  real*8, target,   intent(in)  :: sensor_zenith(nobs) 
  real*8,           intent(in)  :: MISSING 

  real*8, target,   intent(in)  :: albedo(nobs,nch)       ! surface albedo
  real*8,           intent(in)  :: radiance(nobs,nch)     ! radiance observed OMI
!  real*8,           intent(in)  :: reflectivity(nobs,nch) ! reflectivity fichier OMI
 
  
  integer,          intent(in)  :: verbose
  real*8,           intent(in)  :: AI(nobs)                    ! Aerosol Index OMI file
! !OUTPUT PARAMETERS:

  real*8,           intent(out) :: radiance_VL(nobs,nch)       ! TOA radiance from VLIDORT
  integer,          intent(out) :: rc                          ! return code
  real*8,           intent(out) :: AI_VL(nobs)                 ! Aerosol Index model VLIDORT
  real*8,           intent(out) :: refl(nobs, nch)          ! surface reflectivity R*
  
!                               ---
  
  integer             :: n, i,j, k
  real*8              :: spher_alb(nobs,nch)      ! spherical albedo
  real*8              :: trans(nobs, nch)         ! transmission
  real*8              :: reflectivity_388nm       ! Reflectivity Vlidort 388 nm
  real*8              :: reflectivity_354nm       ! Reflectivity Vlidort 354 nm
  real*8              :: reflectivity_388nm_OMI   ! Reflectivity OMI file
  real*8              :: reflectivity_354nm_OMI   ! Reflectivity OMI file
  real*8              :: reflectivity_388nm_adj   ! Reflectivity OMI file
  real*8              :: delta_R_model            ! Difference R
  real*8              :: delta_R_OMI              ! Difference R
  real*8              :: AOT(nobs,nch)            ! total aerosol optical thickness 
  real*8              :: g_tot(nobs,nch)          ! total asymetry factor (column - one profile) 
  real*8              :: sfc_albedo_354nm         ! surface albedo 354 nm OMI file (climato TOMS)
  real*8              :: sfc_albedo_388nm         ! surface albedo 388 nm OMI file (climato TOMS)
  real*8              :: raymodel(nobs)           ! radiance Ray 354 nm from VLIDORT
  real*8              :: ray_OMI                  ! radiance Ray 354 nm OMI file
  real*8              :: reflectance_VL(nobs,nch)    ! TOA reflectance from VLIDORT
  type(VLIDORT_Lamb) :: lamb

  
  rc = 0
  radiance_VL = 0
  AOT =0

  call VLIDORT_LambInit( lamb, km, rc, .true. )

  do j = 1, nobs
     if (abs(solar_zenith(j) - MISSING) < 0.1) then
         do i = 1, nch
         radiance_VL(j,i) = MISSING
         AI_VL(j) = MISSING
         refl(j,i) = MISSING
         end do
     else
     lamb%solar_zenith = solar_zenith (j)
     lamb%relative_azymuth = relat_azymuth(j)
     lamb%sensor_zenith = sensor_zenith(j) 
     lamb%pe => pe(:,j)
     lamb%ze => he(:,j)
     lamb%te => te(:,j) 
    
     do i = 1, nch 
        if (abs(albedo(j,i) - MISSING) < 0.1) then
        radiance_VL(j,i) = MISSING
        refl(j,i) = MISSING
        else
        lamb%wavelength = channels(i)
        lamb%albedo = albedo(j,i)
        lamb%tau => tau(:,i,j)
        lamb%ssa => ssa(:,i,j)
        lamb%g => g(:,i,j)
 
        call VLIDORT_LambRun (lamb, radiance_VL(j,i),reflectance_VL(j,i), AOT(j,i), rc, .true., .true.)
        
        if (abs(lamb%wavelength - 354.00) < 0.01)  sfc_albedo_354nm = albedo(j,i) ! surface albedo  at 354 nm
        
        if (abs(lamb%wavelength - 388.00) < 0.01) then  ! calcul of S, T and R* at 388 nm
        call VLIDORT_LER(lamb, radiance_VL(j,i), refl(j,i), rc)
        sfc_albedo_388nm = albedo(j,i)
        reflectivity_388nm = refl(j,i)
        reflectivity_388nm_adj = refl(j,i)-(sfc_albedo_388nm -sfc_albedo_354nm)    
        end if 
        end if                
     end do
      
     do i = 1, nch 
        if (abs(albedo(j,i) - MISSING) < 0.1) then
        AI_VL(j) = MISSING
        else
        lamb%wavelength = channels(i)
        lamb%albedo = albedo(j,i)
        lamb%tau => tau(:,i,j)
        lamb%ssa => ssa(:,i,j)
        lamb%g => g(:,i,j) 
     
        if (abs(lamb%wavelength - 354.00) < 0.01) then   ! calcul AI only at 354 nm using refl_388nm  
        call VLIDORT_AI(lamb, radiance_VL(j,i),reflectivity_388nm_adj, AI_VL(j), rc)
        end if
        end if     
     end do 
  
     end if
  end do 

end subroutine Scalar

!..........................................................................

subroutine Vector (km, nch, nobs, nymd, nhms, channels, nMom,  &
                   nPol,tau, ssa, g, pmom, pe, he, te,         & 
                   lon, lat, ps, albedo,                       &
                   solar_zenith, relat_azymuth, sensor_zenith, &
                   MISSING,radiance,AI, verbose, radiance_VL,  &
                   AI_VL, refl, rc)
!
! Place holder.
!
   use VLIDORT_LambMod
 
   implicit NONE

! !INPUT PARAMETERS:

  integer,          intent(in)  :: km    ! number of levels on file
  integer,          intent(in)  :: nch   ! number of channels
  integer,          intent(in)  :: nobs  ! number of observations
  integer,          intent(in)  :: nymd  ! date, eg., 20080629
  integer,          intent(in)  :: nhms  ! time, eg., 120000
  integer, target,  intent(in)  :: nMom  ! number of moments 
  integer, target,  intent(in)  :: nPol  ! number of components                               
                  
  real*8, target,   intent(in)  :: channels(nch)    ! wavelengths [nm]

!                                                   ! --- Mie Parameters ---
  real*8, target,   intent(in)  :: tau(km,nch,nobs) ! aerosol optical depth
  real*8, target,   intent(in)  :: ssa(km,nch,nobs) ! single scattering albedo
  real*8, target,   intent(in)  :: g(km,nch,nobs)   ! asymmetry factor
  
  real*8, target,   intent(in)  :: pmom(km,nch,nobs,nMom,nPol) !components of the scat phase matrix


  real*8, target,   intent(in)  :: pe(km+1,nobs)    ! pressure at layer edges [Pa]
  real*8, target,   intent(in)  :: he(km+1,nobs)    ! height above sea-level  [m]
  real*8, target,   intent(in)  :: te(km+1,nobs)    ! temperature at layer edges [K]

                                                    ! --- OMI Metadata ---
  real*8,           intent(in)  :: lon(nobs)        ! longitudes [-180,180]
  real*8,           intent(in)  :: lat(nobs)        ! latitudes [-90,90]
  real*8,           intent(in)  :: ps(nobs)         ! surface pressure [Pa]
                       
  real*8, target,   intent(in)  :: solar_zenith(nobs)  
  real*8, target,   intent(in)  :: relat_azymuth(nobs) 
  real*8, target,   intent(in)  :: sensor_zenith(nobs) 
  real*8,           intent(in)  :: MISSING 


  real*8, target,   intent(in)  :: albedo(nobs,nch)       ! surface albedo
  real*8,           intent(in)  :: radiance(nobs,nch)     ! radiance observed OMI
!  real*8,           intent(in)  :: reflectivity(nobs,nch)! reflectivity fichier OMI
  real*8,           intent(in)  :: AI(nobs)               ! Aerosol Index OMI file
  
  integer,          intent(in)  :: verbose

! !OUTPUT PARAMETERS:

  real*8,           intent(out) :: radiance_VL(nobs,nch)       ! TOA radiance from VLIDORT
  integer,          intent(out) :: rc                          ! return code
  real*8,           intent(out) :: AI_VL(nobs)                 ! Aerosol Index model VLIDORT
  real*8,           intent(out) :: refl(nobs, nch)          ! surface reflectivity R*
!                               ---
  
  integer             :: n, i,j, k, imom, ipol, i_354
  real*8              :: spher_alb(nobs,nch)      ! spherical albedo
  real*8              :: trans(nobs, nch)         ! transmission
  real*8              :: reflectivity_388nm       ! Reflectivity Vlidort 388 nm
  real*8              :: reflectivity_354nm       ! Reflectivity Vlidort 354 nm
  real*8              :: reflectivity_388nm_OMI   ! Reflectivity OMI file
  real*8              :: reflectivity_354nm_OMI   ! Reflectivity OMI file
  real*8              :: reflectivity_388nm_adj   ! Reflectivity OMI file
  real*8              :: delta_R_model            ! Difference R
  real*8              :: delta_R_OMI              ! Difference R
  real*8              :: AOT(nobs,nch)            ! total aerosol optical thickness 
  real*8              :: g_tot(nobs,nch)            ! total asymetry factor 
  real*8              :: sfc_albedo_354nm         ! surface albedo 354 nm OMI file (climato TOMS)
  real*8              :: sfc_albedo_388nm         ! surface albedo 388 nm OMI file (climato TOMS)
  real*8              :: raymodel(nobs)           ! radiance Ray 354 nm from VLIDORT
  real*8              :: ray_OMI                  ! radiance Ray 354 nm OMI file
  real*8              :: reflectance_VL(nobs,nch)    ! TOA reflectance from VLIDORT

  
  type(VLIDORT_Lamb) :: lamb

  rc = 0
  radiance_VL = 0
  AOT =0
 
  call VLIDORT_LambInit( lamb, km, rc, .false. )

  lamb%nmom = nMom
 
  do j = 1, nobs
     if (abs(solar_zenith(j) - MISSING) < 0.1) then
         do i = 1, nch
         radiance_VL(j,i) = MISSING
         AI_VL(j) = MISSING
         refl(j,i) = MISSING
         end do
     else
     lamb%solar_zenith = solar_zenith (j)
     lamb%relative_azymuth = relat_azymuth(j)
     lamb%sensor_zenith = sensor_zenith(j) 
     lamb%pe => pe(:,j)
     lamb%ze => he(:,j)
     lamb%te => te(:,j) 

     do i = 1, nch
        if (abs(albedo(j,i) - MISSING) < 0.1) then
        radiance_VL(j,i) = MISSING
        refl(j,i) = MISSING
        else
        lamb%wavelength = channels(i)
        lamb%albedo = albedo(j,i)
        lamb%tau => tau(:,i,j)
        lamb%ssa => ssa(:,i,j)
        lamb%g => g(:,i,j)
        lamb%pmom => pmom(:,i,j,:,:)
       
        call VLIDORT_LambRun (lamb, radiance_VL(j,i),reflectance_VL(j,i), AOT(j,i), rc, .false., .true.)
       
        if (abs(lamb%wavelength - 354.00) < 0.01) sfc_albedo_354nm = albedo(j,i)! surface albedo at 354 nm

        if (abs(lamb%wavelength - 388.00) < 0.01) then  ! calcul of S, T and R* at 388 nm
        call VLIDORT_LER(lamb, radiance_VL(j,i), refl(j,i), rc)
        sfc_albedo_388nm = albedo(j,i)
        reflectivity_388nm_adj = refl(j,i)-(sfc_albedo_388nm -sfc_albedo_354nm) 
        end if
        end if               
     end do
      
     do i = 1, nch  
        if (abs(albedo(j,i) - MISSING) < 0.1) then
        AI_VL(j) = MISSING
        else
        lamb%wavelength = channels(i)
        if (abs(lamb%wavelength - 354.00) < 0.01) then! calcul AI only at 354 nm using refl_388nm
        lamb%albedo = albedo(j,i)
        lamb%tau => tau(:,i,j)
        lamb%ssa => ssa(:,i,j)
        lamb%g => g(:,i,j) 
        lamb%pmom => pmom(:,i,j,:,:)
     
        call VLIDORT_AI(lamb, radiance_VL(j,i),reflectivity_388nm_adj, AI_VL(j), rc)
        end if 
        end if     
     end do 
  
    end if
  end do
  
end subroutine Vector
