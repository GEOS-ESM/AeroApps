!
!  Simple fortran wrapper for the Python interface to VLIDORT 
!  using a MODIS RTLS BRDF surface
!
!  Patricia Castellanos
!  May 2017
!.............................................................................

subroutine VECTOR_BRDF_MODIS(km, nch, nobs, channels, nMom,  &
                     nPol,tau, ssa, pmom, pe, he, te, kernel_wt, param, &
                     solar_zenith, relat_azymuth, sensor_zenith, &
                     MISSING,verbose, radiance_VL_SURF,reflectance_VL_SURF, ROT, BR, Q, U, rc, DO_2OS_CORRECTION)

    use VLIDORT_BRDF_MODIS  
    implicit None

  ! !INPUT PARAMETERS:

    integer,          intent(in)            :: km    ! number of vertical levels 
    integer,          intent(in)            :: nch   ! number of channels
    integer,          intent(in)            :: nobs  ! number of observations

    integer, target,  intent(in)            :: nMom  ! number of phase function moments 
    integer, target,  intent(in)            :: nPol  ! number of scattering matrix components                               
                    
    real*8, target,   intent(in)            :: channels(nch)    ! wavelengths [nm]

  !                                                   ! --- Aerosol Optical Properties ---
    real*8, target,   intent(in)            :: tau(km,nch,nobs) ! aerosol optical depth
    real*8, target,   intent(in)            :: ssa(km,nch,nobs) ! single scattering albedo    
    real*8, target,   intent(in)            :: pmom(km,nch,nobs,nMom,nPol) !components of the scat phase matrix

    real*8, target,   intent(in)            :: MISSING          ! MISSING VALUE
    real*8, target,   intent(in)            :: pe(km+1,nobs)    ! pressure at layer edges [Pa]
    real*8, target,   intent(in)            :: he(km+1,nobs)    ! height above sea-level  [m]
    real*8, target,   intent(in)            :: te(km+1,nobs)    ! temperature at layer edges [K]

    real*8, target,   intent(in)            :: kernel_wt(nkernel,nch,nobs)    ! kernel weights (/fiso,fgeo,fvol/)
    real*8, target,   intent(in)            :: param(nparam,nch,nobs)         ! Li-Sparse parameters 
                                                                              ! param1 = crown relative height (h/b)
                                                                              ! param2 = shape parameter (b/r)
                         
    real*8, target,   intent(in)            :: solar_zenith(nobs)  
    real*8, target,   intent(in)            :: relat_azymuth(nobs) 
    real*8, target,   intent(in)            :: sensor_zenith(nobs) 

    integer,          intent(in)            :: verbose

  ! !OUTPUT PARAMETERS:
    real*8,           intent(out)           :: radiance_VL_SURF(nobs,nch)     ! TOA normalized radiance from VLIDORT using surface module
    real*8,           intent(out)           :: reflectance_VL_SURF(nobs, nch) ! TOA reflectance from VLIDORT using surface module
    integer,          intent(out)           :: rc                             ! return code

    real*8,           intent(out)           :: ROT(km,nobs,nch)               ! rayleigh optical thickness
    real*8,           intent(out)           :: BR(nobs,nch)                   ! bidirectional reflectance 
    real*8,           intent(out)           :: Q(nobs, nch)                   ! Stokes parameter Q
    real*8,           intent(out)           :: U(nobs, nch)                   ! Stokes parameter U   


    call VLIDORT_Vector_LandMODIS (km, nch, nobs, channels, nMom, &
                                   nPol, tau, ssa, pmom, pe, ze, te, &
                                   kernel_wt, param, &
                                   solar_zenith, &
                                   relat_azymuth, &
                                   sensor_zenith, &
                                   MISSING,verbose,
                                   radiance_VL_int,
                                   reflectance_VL_int, 
                                   ROT, albedo, Q, U, ierr )  


end subroutine VECTOR_BRDF_MODIS