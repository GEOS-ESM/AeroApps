!
!  Simple fortran wrapper for the Python interface to VLIDORT 
!  
!
!  Patricia Castellanos
!  May 2017
!.............................................................................

subroutine ROT_CALC(km, nch, nobs, channels, pe, he, te, &
                     MISSING,verbose, ROT, rc)

    use VLIDORT_ROT, only: VLIDORT_ROT_CALC
    implicit None

  ! !INPUT PARAMETERS:

    integer,          intent(in)            :: km    ! number of vertical levels 
    integer,          intent(in)            :: nch   ! number of channels
    integer,          intent(in)            :: nobs  ! number of observations
                    
    real*8,           intent(in)            :: channels(nch)    ! wavelengths [nm]

    real*8,           intent(in)            :: MISSING          ! MISSING VALUE
    real*8,           intent(in)            :: pe(km+1,nobs)    ! pressure at layer edges [Pa]
    real*8,           intent(in)            :: he(km+1,nobs)    ! height above sea-level  [m]
    real*8,           intent(in)            :: te(km+1,nobs)    ! temperature at layer edges [K]

    integer,          intent(in)            :: verbose

  ! !OUTPUT PARAMETERS:
    integer,          intent(out)           :: rc                             ! return code

    real*8,           intent(out)           :: ROT(km,nobs,nch)               ! rayleigh optical thickness

    call VLIDORT_ROT_CALC (km, nch, nobs, channels, pe, he, te, &
                                   MISSING,verbose, &
                                   ROT, rc )  


end subroutine ROT_CALC



subroutine VECTOR_BRDF_MODIS(km, nch, nobs, channels, nMom, nPol,nkernel,nparam, &
                     ROT, tau, ssa, pmom, pe, he, te, kernel_wt, param, &
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

    integer,          intent(in)            :: nkernel ! number of kernels
    integer,          intent(in)            :: nparam  ! number of kernel parameters
                    
    real*8,           intent(in)            :: channels(nch)    ! wavelengths [nm]

  !                                                   ! --- Rayleigh Parameters ---
    real*8, target,   intent(in)            :: ROT(km,nobs,nch) ! rayleigh optical thickness

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


    call VLIDORT_Vector_LandMODIS (km, nch, nobs, channels, nMom, &
                                   nPol, ROT, tau, ssa, pmom, pe, he, te, &
                                   kernel_wt, param, &
                                   solar_zenith, &
                                   relat_azymuth, &
                                   sensor_zenith, &
                                   MISSING,verbose, &
                                   radiance_VL_SURF, &
                                   reflectance_VL_SURF, &
                                   BR, Q, U, BR_Q, BR_U, rc )  


end subroutine VECTOR_BRDF_MODIS


subroutine VECTOR_BRDF_MODIS_BPDF(km, nch, nobs, channels, nMom, nPol,nkernel,nparam,nparam_bpdf, &
                     ROT, tau, ssa, pmom, pe, he, te, kernel_wt, RTLSparam, BPDFparam, &
                     solar_zenith, relat_azymuth, sensor_zenith, &
                     MISSING,verbose, radiance_VL_SURF,reflectance_VL_SURF, BR, Q, U, BR_Q, BR_U, rc)

    use VLIDORT_BRDF_MODIS_BPDF, only: VLIDORT_Vector_LandMODIS_BPDF 

    implicit None

  ! !INPUT PARAMETERS:

    integer,          intent(in)            :: km    ! number of vertical levels 
    integer,          intent(in)            :: nch   ! number of channels
    integer,          intent(in)            :: nobs  ! number of observations

    integer,          intent(in)            :: nMom  ! number of phase function moments 
    integer,          intent(in)            :: nPol  ! number of scattering matrix components                               

    integer,          intent(in)            :: nkernel ! number of kernels
    integer,          intent(in)            :: nparam  ! number of kernel parameters
    integer,          intent(in)            :: nparam_bpdf  ! number of BPDF kernel parameters
                    
    real*8,           intent(in)            :: channels(nch)    ! wavelengths [nm]

!                                                   ! --- Rayleigh Parameters ---
    real*8, target,   intent(in)            :: ROT(km,nobs,nch) ! rayleigh optical thickness

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


    call VLIDORT_Vector_LandMODIS_BPDF (km, nch, nobs, channels, nMom, &
                                   nPol, ROT, tau, ssa, pmom, pe, he, te, &
                                   kernel_wt, RTLSparam, BPDFparam, &
                                   solar_zenith, &
                                   relat_azymuth, &
                                   sensor_zenith, &
                                   MISSING,verbose, &
                                   radiance_VL_SURF, &
                                   reflectance_VL_SURF, &
                                   BR, Q, U, BR_Q, BR_U, rc )  


end subroutine VECTOR_BRDF_MODIS_BPDF

subroutine VECTOR_BPDF(km, nch, nobs, channels, nMom, nPol,nparam, &
                     ROT, tau, ssa, pmom, pe, he, te, BPDFparam, &
                     solar_zenith, relat_azymuth, sensor_zenith, &
                     MISSING,verbose, radiance_VL_SURF,reflectance_VL_SURF, BR, Q, U, BR_Q, BR_U, rc)

    use VLIDORT_BRDF_MODIS_BPDF, only: VLIDORT_Vector_BPDF 

    implicit None

  ! !INPUT PARAMETERS:

    integer,          intent(in)            :: km    ! number of vertical levels 
    integer,          intent(in)            :: nch   ! number of channels
    integer,          intent(in)            :: nobs  ! number of observations

    integer,          intent(in)            :: nMom  ! number of phase function moments 
    integer,          intent(in)            :: nPol  ! number of scattering matrix components                               

    integer,          intent(in)            :: nparam  ! number of kernel parameters
                    
    real*8,           intent(in)            :: channels(nch)    ! wavelengths [nm]

!                                                   ! --- Rayleigh Parameters ---
    real*8, target,   intent(in)            :: ROT(km,nobs,nch) ! rayleigh optical thickness

  !                                                   ! --- Aerosol Optical Properties ---
    real*8,           intent(in)            :: tau(km,nch,nobs) ! aerosol optical depth
    real*8,           intent(in)            :: ssa(km,nch,nobs) ! single scattering albedo    
    real*8,           intent(in)            :: pmom(km,nch,nobs,nMom,nPol) !components of the scat phase matrix

    real*8,           intent(in)            :: MISSING          ! MISSING VALUE
    real*8,           intent(in)            :: pe(km+1,nobs)    ! pressure at layer edges [Pa]
    real*8,           intent(in)            :: he(km+1,nobs)    ! height above sea-level  [m]
    real*8,           intent(in)            :: te(km+1,nobs)    ! temperature at layer edges [K]

    real*8,           intent(in)            :: BPDFparam(nparam,nch,nobs)     ! BPRDF parameters 
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


    call VLIDORT_Vector_BPDF (km, nch, nobs, channels, nMom, &
                                   nPol, ROT, tau, ssa, pmom, pe, he, te, &
                                   BPDFparam, &
                                   solar_zenith, &
                                   relat_azymuth, &
                                   sensor_zenith, &
                                   MISSING,verbose, &
                                   radiance_VL_SURF, &
                                   reflectance_VL_SURF, &
                                   BR, Q, U, BR_Q, BR_U, rc )  


end subroutine VECTOR_BPDF


subroutine VECTOR_LAMBERT(km, nch, nobs, channels, nMom, nPol, &
                     ROT, tau, ssa, pmom, pe, he, te, albedo, &
                     solar_zenith, relat_azymuth, sensor_zenith, &
                     MISSING,verbose, radiance_VL_SURF,reflectance_VL_SURF, Q, U, rc)

    use VLIDORT_LAMBERT, only: VLIDORT_Vector_Lambert  
    implicit None

  ! !INPUT PARAMETERS:

    integer,          intent(in)            :: km    ! number of vertical levels 
    integer,          intent(in)            :: nch   ! number of channels
    integer,          intent(in)            :: nobs  ! number of observations

    integer,          intent(in)            :: nMom  ! number of phase function moments 
    integer,          intent(in)            :: nPol  ! number of scattering matrix components                               
                    
    real*8,           intent(in)            :: channels(nch)    ! wavelengths [nm]

!                                                   ! --- Rayleigh Parameters ---
    real*8, target,   intent(in)            :: ROT(km,nobs,nch) ! rayleigh optical thickness

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

    integer,          intent(in)            :: verbose

  ! !OUTPUT PARAMETERS:
    real*8,           intent(out)           :: radiance_VL_SURF(nobs,nch)     ! TOA normalized radiance from VLIDORT using surface module
    real*8,           intent(out)           :: reflectance_VL_SURF(nobs, nch) ! TOA reflectance from VLIDORT using surface module
    integer,          intent(out)           :: rc                             ! return code

    real*8,           intent(out)           :: Q(nobs, nch)                   ! Stokes parameter Q
    real*8,           intent(out)           :: U(nobs, nch)                   ! Stokes parameter U   


    call VLIDORT_Vector_Lambert (km, nch, nobs, channels, nMom, &
                                   nPol, ROT, tau, ssa, pmom, pe, he, te, &
                                   albedo, &
                                   solar_zenith, &
                                   relat_azymuth, &
                                   sensor_zenith, &
                                   MISSING,verbose, &
                                   radiance_VL_SURF, &
                                   reflectance_VL_SURF, &
                                   Q, U, rc )  


end subroutine VECTOR_LAMBERT

subroutine VECTOR_LAMBERT_BPDF(km, nch, nobs, channels, nMom, nPol, nparam, &
                     ROT, tau, ssa, pmom, pe, he, te, albedo, BPDFparam, &
                     solar_zenith, relat_azymuth, sensor_zenith, &
                     MISSING,verbose, radiance_VL_SURF,reflectance_VL_SURF, BR, Q, U, BR_Q, BR_U, rc)

    use VLIDORT_LAMBERT_BPDF, only: VLIDORT_Vector_Lambert_BPDF  
    implicit None

  ! !INPUT PARAMETERS:

    integer,          intent(in)            :: km    ! number of vertical levels 
    integer,          intent(in)            :: nch   ! number of channels
    integer,          intent(in)            :: nobs  ! number of observations

    integer,          intent(in)            :: nMom  ! number of phase function moments 
    integer,          intent(in)            :: nPol  ! number of scattering matrix components                               

    integer,          intent(in)            :: nparam  ! number of kernel parameters
                     
    real*8,           intent(in)            :: channels(nch)    ! wavelengths [nm]

!                                                   ! --- Rayleigh Parameters ---
    real*8, target,   intent(in)            :: ROT(km,nobs,nch) ! rayleigh optical thickness

  !                                                   ! --- Aerosol Optical Properties ---
    real*8,           intent(in)            :: tau(km,nch,nobs) ! aerosol optical depth
    real*8,           intent(in)            :: ssa(km,nch,nobs) ! single scattering albedo    
    real*8,           intent(in)            :: pmom(km,nch,nobs,nMom,nPol) !components of the scat phase matrix

    real*8,           intent(in)            :: MISSING          ! MISSING VALUE
    real*8,           intent(in)            :: pe(km+1,nobs)    ! pressure at layer edges [Pa]
    real*8,           intent(in)            :: he(km+1,nobs)    ! height above sea-level  [m]
    real*8,           intent(in)            :: te(km+1,nobs)    ! temperature at layer edges [K]

    real*8, target,   intent(in)            :: albedo(nobs,nch)       ! surface albedo
    real*8,           intent(in)            :: BPDFparam(nparam,nch,nobs)     ! BPRDF parameters 
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


    call VLIDORT_Vector_Lambert_BPDF (km, nch, nobs, channels, nMom, &
                                   nPol, ROT, tau, ssa, pmom, pe, he, te, &
                                   albedo, BPDFparam, &
                                   solar_zenith, &
                                   relat_azymuth, &
                                   sensor_zenith, &
                                   MISSING,verbose, &
                                   radiance_VL_SURF, &
                                   reflectance_VL_SURF, &
                                   BR, Q, U, BR_Q, BR_U, rc )  


end subroutine VECTOR_LAMBERT_BPDF


subroutine SCALAR_LAMBERT(km, nch, nobs, channels, nMom, nPol, &
                     ROT, tau, ssa, g, pmom, pe, he, te, albedo, &
                     solar_zenith, relat_azymuth, sensor_zenith, &
                     MISSING,verbose, radiance_VL_SURF,reflectance_VL_SURF, rc)

    use VLIDORT_LAMBERT, only: VLIDORT_Scalar_Lambert  
    implicit None

  ! !INPUT PARAMETERS:

    integer,          intent(in)            :: km    ! number of vertical levels 
    integer,          intent(in)            :: nch   ! number of channels
    integer,          intent(in)            :: nobs  ! number of observations

    integer,          intent(in)            :: nMom  ! number of phase function moments 
    integer,          intent(in)            :: nPol  ! number of scattering matrix components                               
                    
    real*8,           intent(in)            :: channels(nch)    ! wavelengths [nm]

!                                                   ! --- Rayleigh Parameters ---
    real*8, target,   intent(in)            :: ROT(km,nobs,nch) ! rayleigh optical thickness

  !                                                   ! --- Aerosol Optical Properties ---
    real*8,           intent(in)            :: tau(km,nch,nobs) ! aerosol optical depth
    real*8,           intent(in)            :: ssa(km,nch,nobs) ! single scattering albedo  
    real*8, target,   intent(in)            :: g(km,nch,nobs)   ! asymmetry factor  
    real*8,           intent(in)            :: pmom(km,nch,nobs,nMom,nPol) !components of the scat phase matrix

    real*8,           intent(in)            :: MISSING          ! MISSING VALUE
    real*8,           intent(in)            :: pe(km+1,nobs)    ! pressure at layer edges [Pa]
    real*8,           intent(in)            :: he(km+1,nobs)    ! height above sea-level  [m]
    real*8,           intent(in)            :: te(km+1,nobs)    ! temperature at layer edges [K]

    real*8, target,   intent(in)            :: albedo(nobs,nch)       ! surface albedo
                         
    real*8,           intent(in)            :: solar_zenith(nobs)  
    real*8,           intent(in)            :: relat_azymuth(nobs) 
    real*8,           intent(in)            :: sensor_zenith(nobs) 

    integer,          intent(in)            :: verbose

  ! !OUTPUT PARAMETERS:
    real*8,           intent(out)           :: radiance_VL_SURF(nobs,nch)     ! TOA normalized radiance from VLIDORT using surface module
    real*8,           intent(out)           :: reflectance_VL_SURF(nobs, nch) ! TOA reflectance from VLIDORT using surface module
    integer,          intent(out)           :: rc                             ! return code



    call VLIDORT_Scalar_Lambert (km, nch, nobs, channels, nMom, &
                                   nPol, ROT, tau, ssa, g, pmom, pe, he, te, &
                                   albedo, &
                                   solar_zenith, &
                                   relat_azymuth, &
                                   sensor_zenith, &
                                   MISSING,verbose, &
                                   radiance_VL_SURF, &
                                   reflectance_VL_SURF, &
                                   rc )  


end subroutine SCALAR_LAMBERT

subroutine VECTOR_GissCX(km, nch, nobs, channels, nMom,  &
                     nPol, ROT, tau, ssa, pmom, pe, he, te, U10m, V10m, &
                     mr, solar_zenith, relat_azymuth, sensor_zenith, &
                     MISSING,verbose, radiance_VL_SURF,reflectance_VL_SURF, &
                     BR, Q, U, BR_Q, BR_U, rc)

    use VLIDORT_BRDF_CX, only: VLIDORT_Vector_GissCX  
    implicit None

  ! !INPUT PARAMETERS:

    integer,          intent(in)            :: km    ! number of vertical levels 
    integer,          intent(in)            :: nch   ! number of channels
    integer,          intent(in)            :: nobs  ! number of observations

    integer,          intent(in)            :: nMom  ! number of phase function moments 
    integer,          intent(in)            :: nPol  ! number of scattering matrix components                               
                    
    real*8,           intent(in)            :: channels(nch)    ! wavelengths [nm]

!                                                   ! --- Rayleigh Parameters ---
    real*8, target,   intent(in)            :: ROT(km,nobs,nch) ! rayleigh optical thickness

  !                                                   ! --- Aerosol Optical Properties ---
    real*8,           intent(in)            :: tau(km,nch,nobs) ! aerosol optical depth
    real*8,           intent(in)            :: ssa(km,nch,nobs) ! single scattering albedo    
    real*8,           intent(in)            :: pmom(km,nch,nobs,nMom,nPol) !components of the scat phase matrix

    real*8,           intent(in)            :: MISSING          ! MISSING VALUE
    real*8,           intent(in)            :: pe(km+1,nobs)    ! pressure at layer edges [Pa]
    real*8,           intent(in)            :: he(km+1,nobs)    ! height above sea-level  [m]
    real*8,           intent(in)            :: te(km+1,nobs)    ! temperature at layer edges [K]

    real*8,           intent(in)            :: U10m(nobs)   ! Wind speed components [m/s]
    real*8,           intent(in)            :: V10m(nobs)    
    real*8,           intent(in)            :: mr(nch)       ! refractive index
                         
    real*8,           intent(in)            :: solar_zenith(nobs)  
    real*8,           intent(in)            :: relat_azymuth(nobs) 
    real*8,           intent(in)            :: sensor_zenith(nobs) 

    integer,          intent(in)            :: verbose

  ! !OUTPUT PARAMETERS:
    real*8,           intent(out)           :: radiance_VL_SURF(nobs,nch)     ! TOA normalized radiance from VLIDORT using surface module
    real*8,           intent(out)           :: reflectance_VL_SURF(nobs, nch) ! TOA reflectance from VLIDORT using surface module
    integer,          intent(out)           :: rc                             ! return code

    real*8,           intent(out)           :: BR(nobs,nch)                   ! polarized bidirectional reflectance
    real*8,           intent(out)           :: BR_Q(nobs,nch)                 ! polarized bidirectional reflectance
    real*8,           intent(out)           :: BR_U(nobs,nch)                 ! polarized bidirectional reflectance    
    real*8,           intent(out)           :: Q(nobs, nch)                   ! Stokes parameter Q
    real*8,           intent(out)           :: U(nobs, nch)                   ! Stokes parameter U   


    call VLIDORT_Vector_GissCX (km, nch, nobs, channels, nMom, &
                                   nPol, ROT, tau, ssa, pmom, pe, he, te, &
                                   U10m, V10m, mr, &
                                   solar_zenith, &
                                   relat_azymuth, &
                                   sensor_zenith, &
                                   MISSING,verbose, &
                                   radiance_VL_SURF, &
                                   reflectance_VL_SURF, &
                                   BR, Q, U, BR_Q, BR_U, rc )  


end subroutine VECTOR_GissCX

subroutine VECTOR_CX(km, nch, nobs, channels, nMom,  &
                     nPol, ROT, tau, ssa, pmom, pe, he, te, U10m, V10m, &
                     mr, solar_zenith, relat_azymuth, sensor_zenith, &
                     MISSING,verbose, radiance_VL_SURF,reflectance_VL_SURF, &
                     BR, Q, U, BR_Q, BR_U, rc)

    use VLIDORT_BRDF_CX, only: VLIDORT_Vector_CX  
    implicit None

  ! !INPUT PARAMETERS:

    integer,          intent(in)            :: km    ! number of vertical levels 
    integer,          intent(in)            :: nch   ! number of channels
    integer,          intent(in)            :: nobs  ! number of observations

    integer,          intent(in)            :: nMom  ! number of phase function moments 
    integer,          intent(in)            :: nPol  ! number of scattering matrix components                               
                    
    real*8,           intent(in)            :: channels(nch)    ! wavelengths [nm]

!                                                   ! --- Rayleigh Parameters ---
    real*8, target,   intent(in)            :: ROT(km,nobs,nch) ! rayleigh optical thickness

  !                                                   ! --- Aerosol Optical Properties ---
    real*8,           intent(in)            :: tau(km,nch,nobs) ! aerosol optical depth
    real*8,           intent(in)            :: ssa(km,nch,nobs) ! single scattering albedo    
    real*8,           intent(in)            :: pmom(km,nch,nobs,nMom,nPol) !components of the scat phase matrix

    real*8,           intent(in)            :: MISSING          ! MISSING VALUE
    real*8,           intent(in)            :: pe(km+1,nobs)    ! pressure at layer edges [Pa]
    real*8,           intent(in)            :: he(km+1,nobs)    ! height above sea-level  [m]
    real*8,           intent(in)            :: te(km+1,nobs)    ! temperature at layer edges [K]

    real*8,           intent(in)            :: U10m(nobs)   ! Wind speed components [m/s]
    real*8,           intent(in)            :: V10m(nobs)    
    real*8,           intent(in)            :: mr(nch)       ! refractive index
                         
    real*8,           intent(in)            :: solar_zenith(nobs)  
    real*8,           intent(in)            :: relat_azymuth(nobs) 
    real*8,           intent(in)            :: sensor_zenith(nobs) 

    integer,          intent(in)            :: verbose

  ! !OUTPUT PARAMETERS:
    real*8,           intent(out)           :: radiance_VL_SURF(nobs,nch)     ! TOA normalized radiance from VLIDORT using surface module
    real*8,           intent(out)           :: reflectance_VL_SURF(nobs, nch) ! TOA reflectance from VLIDORT using surface module
    integer,          intent(out)           :: rc                             ! return code

    real*8,           intent(out)           :: BR(nobs,nch)                   ! polarized bidirectional reflectance
    real*8,           intent(out)           :: BR_Q(nobs,nch)                 ! polarized bidirectional reflectance
    real*8,           intent(out)           :: BR_U(nobs,nch)                 ! polarized bidirectional reflectance    
    real*8,           intent(out)           :: Q(nobs, nch)                   ! Stokes parameter Q
    real*8,           intent(out)           :: U(nobs, nch)                   ! Stokes parameter U   


    call VLIDORT_Vector_CX (km, nch, nobs, channels, nMom, &
                                   nPol, ROT, tau, ssa, pmom, pe, he, te, &
                                   U10m, V10m, mr, &
                                   solar_zenith, &
                                   relat_azymuth, &
                                   sensor_zenith, &
                                   MISSING,verbose, &
                                   radiance_VL_SURF, &
                                   reflectance_VL_SURF, &
                                   BR, Q, U, BR_Q, BR_U, rc )  


end subroutine VECTOR_CX
