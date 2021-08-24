!
!  Simple fortran wrapper for the Python interface to VLIDORT 
!  
!
!  Patricia Castellanos
!  May 2017
!.............................................................................

subroutine VECTOR_BRDF_MODIS(km, nch, nobs, channels, nstreams, plane_parallel, nMom, nPol, nkernel,nparam, &
                     ROT, depol, &
                     tau, ssa, pmom, pe, he, te, kernel_wt, param, &
                     solar_zenith, relat_azymuth, sensor_zenith, &
                     MISSING,verbose, radiance_VL_SURF,reflectance_VL_SURF, BR, Q, U, BR_Q, BR_U, rc)

    use VLIDORT_BRDF_MODIS, only: VLIDORT_Vector_LandMODIS  
    implicit None

  ! !INPUT PARAMETERS:

    integer,          intent(in)            :: km    ! number of vertical levels 
    integer,          intent(in)            :: nch   ! number of channels
    integer,          intent(in)            :: nobs  ! number of observations

    integer,          intent(in)            :: nMom  ! number of phase function moments 
    integer,          intent(in)            :: nPol  ! number of scattering matrix components                               
    logical,          intent(in)            :: plane_parallel ! do plane parallel flag
    integer,          intent(in)            :: nstreams  ! number of half space streams

    integer,          intent(in)            :: nkernel ! number of kernels
    integer,          intent(in)            :: nparam  ! number of kernel parameters

    real*8,           intent(in)            :: channels(nch)    ! wavelengths [nm]

!                                                     ! --- Rayleigh Parameters ---
    real*8, target,   intent(in)            :: ROT(km,nobs,nch) ! rayleigh optical thickness
    real*8, target,   intent(in)            :: depol(nch)       ! rayleigh depolarization ratio used in phase matrix

  !                                                   ! --- Aerosol Optical Properties ---
    real*8,           intent(in)            :: tau(km,nch,nobs) ! aerosol optical depth
    real*8,           intent(in)            :: ssa(km,nch,nobs) ! single scattering albedo    
    real*8,           intent(in)            :: pmom(km,nch,nobs,nMom,nPol) !components of the scat phase matrix

    real*8,           intent(in)            :: MISSING          ! MISSING VALUE
    real*8,           intent(in)            :: pe(km+1,nobs)    ! pressure at layer edges [Pa]
    real*8,           intent(in)            :: he(km+1,nobs)    ! height above sea-level  [m]
    real*8,           intent(in)            :: te(km+1,nobs)    ! temperature at layer edges [K]

    real*8,           intent(in)            :: kernel_wt(nkernel,nch,nobs)    ! kernel weights (/fiso,fgeo,fvol/)
    real*8,           intent(in)            :: param(nparam,nch,nobs)         ! Li-Sparse parameters 
                                                                              ! param1 = crown relative height (h/b)
                                                                              ! param2 = shape parameter (b/r)
                         
    real*8,           intent(in)            :: solar_zenith(nobs)  
    real*8,           intent(in)            :: relat_azymuth(nobs) 
    real*8,           intent(in)            :: sensor_zenith(nobs) 

    integer,          intent(in)            :: verbose

  ! !OUTPUT PARAMETERS:
    real*8,           intent(out)           :: radiance_VL_SURF(nobs,nch)     ! TOA normalized radiance from VLIDORT using surface module
    real*8,           intent(out)           :: reflectance_VL_SURF(nobs, nch) ! TOA reflectance from VLIDORT using surface module
    integer,          intent(out)           :: rc                             ! return code

    real*8,           intent(out)           :: BR(nobs,nch)                   ! bidirectional reflectance 
    real*8,           intent(out)           :: Q(nobs, nch)                   ! Stokes parameter Q
    real*8,           intent(out)           :: U(nobs, nch)                   ! Stokes parameter U   
    real*8,           intent(out)           :: BR_Q(nobs, nch)                ! Stokes parameter Q
    real*8,           intent(out)           :: BR_U(nobs, nch)                ! Stokes parameter U

    logical, parameter                      :: DO_2OS_CORRECTION = .false.
    logical, parameter                      :: DO_BOA = .true.


    call VLIDORT_Vector_LandMODIS (km, nch, nobs, channels, nstreams, &
                                   plane_parallel, nMom, &
                                   nPol, ROT, depol, tau, ssa, pmom, pe, he, te, &
                                   kernel_wt, param, &
                                   solar_zenith, &
                                   relat_azymuth, &
                                   sensor_zenith, &
                                   MISSING,verbose, &
                                   radiance_VL_SURF, &
                                   reflectance_VL_SURF, &
                                   BR, Q, U, BR_Q, BR_U, rc, DO_2OS_CORRECTION, DO_BOA)  


end subroutine VECTOR_BRDF_MODIS


subroutine VECTOR_BRDF_MODIS_BPDF(km, nch, nobs, channels, nstreams, plane_parallel, nMom, nPol,nkernel,nparam,nparam_bpdf, &
                     ROT, depol, tau, ssa, pmom, pe, he, te, kernel_wt, RTLSparam, BPDFparam, &
                     solar_zenith, relat_azymuth, sensor_zenith, &
                     MISSING,verbose, radiance_VL_SURF,reflectance_VL_SURF, BR, Q, U, BR_Q, BR_U, rc)

    use VLIDORT_BRDF_MODIS_BPDF, only: VLIDORT_Vector_LandMODIS_BPDF 

    implicit None

  ! !INPUT PARAMETERS:

    integer,          intent(in)            :: km    ! number of vertical levels 
    integer,          intent(in)            :: nch   ! number of channels
    integer,          intent(in)            :: nobs  ! number of observations

    logical,          intent(in)            :: plane_parallel ! do plane parallel flag
    integer,          intent(in)            :: nstreams  ! number of half space streams

    integer,          intent(in)            :: nMom  ! number of phase function moments 
    integer,          intent(in)            :: nPol  ! number of scattering matrix components                               

    integer,          intent(in)            :: nkernel ! number of kernels
    integer,          intent(in)            :: nparam  ! number of kernel parameters
    integer,          intent(in)            :: nparam_bpdf  ! number of BPDF kernel parameters

    real*8,           intent(in)            :: channels(nch)    ! wavelengths [nm]

!                                                     ! --- Rayleigh Parameters ---
    real*8, target,   intent(in)            :: ROT(km,nobs,nch) ! rayleigh optical thickness
    real*8, target,   intent(in)            :: depol(nch)       ! rayleigh depolarization ratio used in phase matrix

  !                                                   ! --- Aerosol Optical Properties ---
    real*8,           intent(in)            :: tau(km,nch,nobs) ! aerosol optical depth
    real*8,           intent(in)            :: ssa(km,nch,nobs) ! single scattering albedo    
    real*8,           intent(in)            :: pmom(km,nch,nobs,nMom,nPol) !components of the scat phase matrix

    real*8,           intent(in)            :: MISSING          ! MISSING VALUE
    real*8,           intent(in)            :: pe(km+1,nobs)    ! pressure at layer edges [Pa]
    real*8,           intent(in)            :: he(km+1,nobs)    ! height above sea-level  [m]
    real*8,           intent(in)            :: te(km+1,nobs)    ! temperature at layer edges [K]

    real*8,           intent(in)            :: kernel_wt(nkernel,nch,nobs)    ! kernel weights (/fiso,fgeo,fvol/)
    real*8,           intent(in)            :: RTLSparam(nparam,nch,nobs)         ! Li-Sparse parameters 
                                                                              ! param1 = crown relative height (h/b)
                                                                              ! param2 = shape parameter (b/r)
    real*8,           intent(in)            :: BPDFparam(nparam_bpdf,nch,nobs)     ! BPRDF parameters 
                                                                              ! param1 = Refractive index of water (1.5)
                                                                              ! param2 = NDVI , must be <=1 or >=1
                                                                              ! param3 = Scaling factor (C, from Maignan 2009 RSE Table 1)                         
    real*8,           intent(in)            :: solar_zenith(nobs)  
    real*8,           intent(in)            :: relat_azymuth(nobs) 
    real*8,           intent(in)            :: sensor_zenith(nobs) 

    integer,          intent(in)            :: verbose

  ! !OUTPUT PARAMETERS:
    real*8,           intent(out)           :: radiance_VL_SURF(nobs,nch)     ! TOA normalized radiance from VLIDORT using surface module
    real*8,           intent(out)           :: reflectance_VL_SURF(nobs, nch) ! TOA reflectance from VLIDORT using surface module
    integer,          intent(out)           :: rc                             ! return code

    real*8,           intent(out)           :: BR(nobs,nch)                   ! bidirectional reflectance 
    real*8,           intent(out)           :: Q(nobs, nch)                   ! Stokes parameter Q
    real*8,           intent(out)           :: U(nobs, nch)                   ! Stokes parameter U   
    real*8,           intent(out)           :: BR_Q(nobs, nch)                ! Stokes parameter Q
    real*8,           intent(out)           :: BR_U(nobs, nch)                ! Stokes parameter U


    logical, parameter                      :: DO_BOA = .true.

    call VLIDORT_Vector_LandMODIS_BPDF (km, nch, nobs, channels, nstreams, plane_parallel, nMom, &
                                   nPol, ROT, depol, tau, ssa, pmom, pe, he, te, &
                                   kernel_wt, RTLSparam, BPDFparam, &
                                   solar_zenith, &
                                   relat_azymuth, &
                                   sensor_zenith, &
                                   MISSING,verbose, &
                                   radiance_VL_SURF, &
                                   reflectance_VL_SURF, &
                                   BR, Q, U, BR_Q, BR_U, rc, DO_BOA )  


end subroutine VECTOR_BRDF_MODIS_BPDF

subroutine VECTOR_LAMBERT(km, nch, nobs, channels, nstreams, plane_parallel, nMom, nPol, &
                     ROT, depol, alpha, tau, ssa, pmom, pe, he, te, albedo, &
                     solar_zenith, relat_azymuth, sensor_zenith, flux_factor, &
                     MISSING,verbose, radiance_VL_SURF,reflectance_VL_SURF, Q, U, rc)

    use VLIDORT_LAMBERT, only: VLIDORT_Vector_Lambert  
    implicit None

  ! !INPUT PARAMETERS:

    integer,          intent(in)            :: km    ! number of vertical levels 
    integer,          intent(in)            :: nch   ! number of channels
    integer,          intent(in)            :: nobs  ! number of observations

    logical,          intent(in)            :: plane_parallel ! do plane parallel flag

    integer,          intent(in)            :: nMom  ! number of phase function moments 
    integer,          intent(in)            :: nPol  ! number of scattering matrix components                               
    
    integer,          intent(in)            :: nstreams  ! number of half space streams

    real*8,           intent(in)            :: channels(nch)    ! wavelengths [nm]

!                                                   ! --- Rayleigh Parameters ---
    real*8,           intent(in)            :: ROT(km,nobs,nch) ! rayleigh optical thickness
    real*8,           intent(in)            :: depol(nch)       ! rayleigh depolarization ratio

    real*8,           intent(in)            :: alpha(km,nobs,nch)       ! trace gas absorption

  !                                                   ! --- Aerosol Optical Properties ---
    real*8,           intent(in)            :: tau(km,nch,nobs) ! aerosol optical depth
    real*8,           intent(in)            :: ssa(km,nch,nobs) ! single scattering albedo    
    real*8,           intent(in)            :: pmom(km,nch,nobs,nMom,nPol) !components of the scat phase matrix

    real*8,           intent(in)            :: MISSING          ! MISSING VALUE
    real*8,           intent(in)            :: pe(km+1,nobs)    ! pressure at layer edges [Pa]
    real*8,           intent(in)            :: he(km+1,nobs)    ! height above sea-level  [m]
    real*8,           intent(in)            :: te(km+1,nobs)    ! temperature at layer edges [K]

    real*8, target,   intent(in)            :: albedo(nobs,nch)       ! surface albedo
                         
    real*8,           intent(in)            :: solar_zenith(nobs)  
    real*8,           intent(in)            :: relat_azymuth(nobs) 
    real*8,           intent(in)            :: sensor_zenith(nobs) 

    real*8,           intent(in)  :: flux_factor(nch,nobs) ! solar flux (F0)

    integer,          intent(in)            :: verbose

  ! !OUTPUT PARAMETERS:
    real*8,           intent(out)           :: radiance_VL_SURF(nobs,nch)     ! TOA normalized radiance from VLIDORT using surface module
    real*8,           intent(out)           :: reflectance_VL_SURF(nobs, nch) ! TOA reflectance from VLIDORT using surface module
    integer,          intent(out)           :: rc                             ! return code

    real*8,           intent(out)           :: Q(nobs, nch)                   ! Stokes parameter Q
    real*8,           intent(out)           :: U(nobs, nch)                   ! Stokes parameter U   

    logical, parameter                      :: DO_2OS_CORRECTION = .false.
    logical, parameter                      :: DO_BOA = .true.


    call VLIDORT_Vector_Lambert (km, nch, nobs, channels,  nstreams, plane_parallel, nMom, &
                                   nPol, ROT, depol, alpha, tau, ssa, pmom, pe, he, te, &
                                   albedo, &
                                   solar_zenith, &
                                   relat_azymuth, &
                                   sensor_zenith, &
                                   flux_factor, &
                                   MISSING,verbose, &
                                   radiance_VL_SURF, &
                                   reflectance_VL_SURF, &
                                   Q, U, rc, DO_2OS_CORRECTION, DO_BOA )  


end subroutine VECTOR_LAMBERT
