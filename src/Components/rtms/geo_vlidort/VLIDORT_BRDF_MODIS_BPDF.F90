module VLIDORT_BRDF_MODIS_BPDF  
!
!  Simple f90 wrapper for the Python interface to VLIDORT for 
!  multiangle polarimeter simulator
!
!.............................................................................
  implicit NONE

  PUBLIC VLIDORT_Vector_LandMODIS_BPDF
  PUBLIC VLIDORT_Vector_LandMODIS_BPDF_Cloud
  PUBLIC VLIDORT_Vector_BPDF

  contains

  logical function IS_MISSING(x,MISSING)
    real*8, intent(in)     :: x
    real*8, intent(in)     :: MISSING

    IS_MISSING = abs(x/MISSING-1)<0.001
    return
  end function IS_MISSING

  subroutine VLIDORT_Vector_LandMODIS_BPDF (km, nch, nobs, channels, nstreams, plane_parallel, nMom,  &
                     nPol, ROT, depol, tau, ssa, pmom, pe, he, te, kernel_wt, RTLSparam, BPDFparam, &
                     solar_zenith, relat_azymuth, sensor_zenith, &
                     MISSING,verbose, radiance_VL_SURF,reflectance_VL_SURF, BR, Q, U, BR_Q, BR_U, rc, DO_BOA)
  !
  ! Place holder.
  !
     use VLIDORT_ScatMod
   
     implicit NONE
 
    logical                                 :: scalar
    integer, parameter                      :: nkernel = 3   ! number of kernels
    integer, parameter                      :: nparam  = 2   ! number of kernel parameters
    integer, parameter                      :: nparam_bpdf  = 3   ! number of BPDF kernel parameters

  ! !INPUT PARAMETERS:

    integer,          intent(in)            :: km    ! number of vertical levels 
    integer,          intent(in)            :: nch   ! number of channels
    integer,          intent(in)            :: nobs  ! number of observations

    logical,          intent(in)            :: plane_parallel ! do plane parallel flag

    integer, target,  intent(in)            :: nMom  ! number of phase function moments 
    integer, target,  intent(in)            :: nPol  ! number of scattering matrix components                               
    integer,          intent(in)            :: nstreams  ! number of half space streams
                    
    real*8, target,   intent(in)            :: channels(nch)    ! wavelengths [nm]

!                                                     ! --- Rayleigh Parameters ---
    real*8, target,   intent(in)            :: ROT(km,nobs,nch) ! rayleigh optical thickness
    real*8, target,   intent(in)            :: depol(nch)       ! rayleigh depolarization ratio used in phase matrix

  !                                                   ! --- Aerosol Optical Properties ---
    real*8, target,   intent(in)            :: tau(km,nch,nobs) ! aerosol optical depth
    real*8, target,   intent(in)            :: ssa(km,nch,nobs) ! single scattering albedo    
    real*8, target,   intent(in)            :: pmom(km,nch,nobs,nMom,nPol) !components of the scat phase matrix

    real*8, target,   intent(in)            :: MISSING          ! MISSING VALUE
    real*8, target,   intent(in)            :: pe(km+1,nobs)    ! pressure at layer edges [Pa]
    real*8, target,   intent(in)            :: he(km+1,nobs)    ! height above sea-level  [m]
    real*8, target,   intent(in)            :: te(km+1,nobs)    ! temperature at layer edges [K]

    real*8, target,   intent(in)            :: kernel_wt(nkernel,nch,nobs)    ! kernel weights (/fiso,fgeo,fvol/)
    real*8, target,   intent(in)            :: RTLSparam(nparam,nch,nobs)     ! Li-Sparse parameters 
                                                                              ! param1 = crown relative height (h/b)
                                                                              ! param2 = shape parameter (b/r)
    real*8, target,   intent(in)            :: BPDFparam(nparam_bpdf,nch,nobs)     ! BPRDF parameters 
                                                                              ! param1 = Refractive index of water (1.5)
                                                                              ! param2 = NDVI , must be <=1 or >=1
                                                                              ! param3 = Scaling factor (C, from Maignan 2009 RSE Table 1)
                         
    real*8, target,   intent(in)            :: solar_zenith(nobs)  
    real*8, target,   intent(in)            :: relat_azymuth(nobs) 
    real*8, target,   intent(in)            :: sensor_zenith(nobs) 

    integer,          intent(in)            :: verbose

  ! !OUTPUT PARAMETERS:
    real*8,           intent(out)           :: radiance_VL_SURF(nobs,nch)     ! TOA normalized radiance from VLIDORT using surface module
    real*8,           intent(out)           :: reflectance_VL_SURF(nobs, nch) ! TOA reflectance from VLIDORT using surface module
    integer,          intent(out)           :: rc                             ! return code

    real*8,           intent(out)           :: BR(nobs,nch)                   ! bidirectional reflectance 
    real*8,           intent(out)           :: BR_Q(nobs,nch)                   ! bidirectional reflectance 
    real*8,           intent(out)           :: BR_U(nobs,nch)                   ! bidirectional reflectance 
    real*8,           intent(out)           :: Q(nobs, nch)                   ! Stokes parameter Q
    real*8,           intent(out)           :: U(nobs, nch)                   ! Stokes parameter U   

  ! !OPTIONAL PARAMETERS
    logical, optional, intent(in) :: DO_BOA
  !                               ---
    
    integer             :: i,j,n,p,ier

    
    type(VLIDORT_scat) :: SCAT
    type(VLIDORT_output_vector)  :: output

    rc = 0
    ier = 0
    if (present(DO_BOA)) SCAT%DO_BOA = DO_BOA
    SCAT%Surface%Base%NSTREAMS = nstreams
    SCAT%Surface%Base%DO_PLANE_PARALLEL = plane_parallel
    SCAT%Surface%Base%NGREEK_MOMENTS_INPUT = nMom
    call VLIDORT_Init( SCAT%Surface%Base, km, rc, SCAT%DO_BOA)
    if ( rc /= 0 ) return

    SCAT%nMom    = nMom
    SCAT%nPol    = nPol
    SCAT%NSTOKES = 3

    if ( SCAT%NSTOKES  .GT. MAXSTOKES  )   return


    do j = 1, nobs
       
      ! Make sure angles are available
      ! -----------------------------------------
      if ( IS_MISSING(solar_zenith(j),MISSING)  .OR. & 
           IS_MISSING(sensor_zenith(j),MISSING) .OR. &
           IS_MISSING(relat_azymuth(j),MISSING)  )  then

        radiance_VL_SURF(j,:) = MISSING
        reflectance_VL_SURF(j,:) = MISSING
        
        cycle

      end if

      SCAT%pe => pe(:,j)
      SCAT%ze => he(:,j)
      SCAT%te => te(:,j) 

      do i = 1, nch
        ! Make sure kernel weights and parameters are defined
        
        do n = 1, nkernel
          
          if (IS_MISSING(kernel_wt(n,i,j),MISSING)) then
            radiance_VL_SURF(j,i) = MISSING
            reflectance_VL_SURF(j,i) = MISSING
            cycle
          end if
        end do
        
        do p = 1, nparam
          if (IS_MISSING(RTLSparam(p,i,j),MISSING)) then
            radiance_VL_SURF(j,i) = MISSING
            reflectance_VL_SURF(j,i) = MISSING
            cycle
          end if
        end do

        do p = 1, nparam_bpdf
          if (IS_MISSING(BPDFparam(p,i,j),MISSING)) then
            radiance_VL_SURF(j,i) = MISSING
            reflectance_VL_SURF(j,i) = MISSING
            cycle
          end if
        end do        

        if ( verbose > 0 ) then
          print*, 'DO MODIS BRDF + BPDF'
        end if

        scalar = .false.
        call VLIDORT_LANDMODIS_BPDF(SCAT%Surface,solar_zenith(j),&
                               sensor_zenith(j),relat_azymuth(j),&
                               kernel_wt(1,i,j),kernel_wt(2,i,j),kernel_wt(3,i,j),&
                               reshape(RTLSparam(:,i,j),(/nparam/)),reshape(BPDFparam(:,i,j),(/nparam_bpdf/)),&
                               scalar,rc)

        SCAT%wavelength = channels(i)
        SCAT%rot => ROT(:,j,i)  
        SCAT%depol_ratio => depol(i)      
        SCAT%tau => tau(:,i,j)
        SCAT%ssa => ssa(:,i,j)
        SCAT%pmom => pmom(:,i,j,:,:)

        BR(j,i) = SCAT%Surface%Base%VIO%VBRDF_Sup_Out%BS_DBOUNCE_BRDFUNC(1,1,1,1)
        BR_Q(j,i) = SCAT%Surface%Base%VIO%VBRDF_Sup_Out%BS_DBOUNCE_BRDFUNC(2,1,1,1) 
        BR_U(j,i) = SCAT%Surface%Base%VIO%VBRDF_Sup_Out%BS_DBOUNCE_BRDFUNC(3,1,1,1)     

        ! Check to make sure kernel combination is not unphysical
        if (BR(j,i) < 0) then
          radiance_VL_SURF(j,i)    = -500
          reflectance_VL_SURF(j,i) = -500
          BR(j,i)    = -500
          BR_Q(j,i)  = -500
          BR_U(j,i)  = -500
          Q(j,i)     = -500
          U(j,i)     = -500
          cycle
        end if

        call VLIDORT_Run_Vector (SCAT, output, ier)

        if (SCAT%DO_BOA) then
          radiance_VL_SURF(j,i)    = output%BOA_radiance
          reflectance_VL_SURF(j,i) = output%BOA_reflectance
          Q(j,i)                   = output%BOA_Q
          U(j,i)                   = output%BOA_U                          
        else
          radiance_VL_SURF(j,i)    = output%radiance
          reflectance_VL_SURF(j,i) = output%reflectance
          Q(j,i)                   = output%Q
          U(j,i)                   = output%U                
        end if
        
        if ( ier /= 0 ) then
          radiance_VL_SURF(j,i) = MISSING
          reflectance_VL_SURF(j,i) = MISSING               
          cycle
        end if

      end do ! end loop over channels
       
      if ( verbose > 0 ) then
        if ( mod(j-1,1000) == 0 ) then
          print *, '<> VLIDORT Vector: ', nint(j*100./nobs), '%'
        end if
      end if

    end do ! Loop over obs

  end subroutine VLIDORT_Vector_LandMODIS_BPDF

  subroutine VLIDORT_Vector_LandMODIS_BPDF_Cloud (km, nch, nobs, channels, nstreams, plane_parallel, nMom,  &
                     nPol, ROT, depol, alpha, tau, ssa, pmom, tauI, ssaI, pmomI, tauL, ssaL, pmomL, &
                     pe, he, te, kernel_wt, RTLSparam, BPDFparam, &
                     solar_zenith, relat_azymuth, sensor_zenith, flux_factor, &
                     MISSING,verbose, radiance_VL_SURF,reflectance_VL_SURF, BR, Q, U, BR_Q, BR_U, rc, DO_BOA)
  !
  ! Place holder.
  !
     use VLIDORT_ScatMod
   
     implicit NONE
 
    logical                                 :: scalar
    integer, parameter                      :: nkernel = 3   ! number of kernels
    integer, parameter                      :: nparam  = 2   ! number of kernel parameters
    integer, parameter                      :: nparam_bpdf  = 3   ! number of BPDF kernel parameters

  ! !INPUT PARAMETERS:

    integer,          intent(in)            :: km    ! number of vertical levels 
    integer,          intent(in)            :: nch   ! number of channels
    integer,          intent(in)            :: nobs  ! number of observations

    logical,          intent(in)            :: plane_parallel ! do plane parallel flag

    integer, target,  intent(in)            :: nMom  ! number of phase function moments 
    integer, target,  intent(in)            :: nPol  ! number of scattering matrix components                               
    integer,          intent(in)            :: nstreams  ! number of half space streams
                    
    real*8, target,   intent(in)            :: channels(nch)    ! wavelengths [nm]

!                                                     ! --- Rayleigh Parameters ---
    real*8, target,   intent(in)            :: ROT(km,nobs,nch) ! rayleigh optical thickness
    real*8, target,   intent(in)            :: depol(nch)       ! rayleigh depolarization ratio used in phase matrix

!                                                   ! --- Trace Gas Absorption ---
  real*8, target,   intent(in)  :: alpha(km,nobs,nch) ! trace gas absoprtion optical thickness

  !                                                   ! --- Aerosol Optical Properties ---
    real*8, target,   intent(in)            :: tau(km,nch,nobs) ! aerosol optical depth
    real*8, target,   intent(in)            :: ssa(km,nch,nobs) ! single scattering albedo    
    real*8, target,   intent(in)            :: pmom(km,nch,nobs,nMom,nPol) !components of the scat phase matrix

    real*8, target,   intent(in)            :: tauI(km,nch,nobs) ! ice cloud optical depth
    real*8, target,   intent(in)            :: ssaI(km,nch,nobs) ! ice cloud single scattering albedo   
    real*8, target,   intent(in)            :: pmomI(km,nch,nobs,nMom,nPol) !ice cloud components of the scat phase matrix

    real*8, target,   intent(in)            :: tauL(km,nch,nobs) ! liquid cloud optical depth
    real*8, target,   intent(in)            :: ssaL(km,nch,nobs) ! liquid cloud single scattering albedo   
    real*8, target,   intent(in)            :: pmomL(km,nch,nobs,nMom,nPol) !liquid cloud components of the scat phase matrix

    real*8, target,   intent(in)            :: MISSING          ! MISSING VALUE
    real*8, target,   intent(in)            :: pe(km+1,nobs)    ! pressure at layer edges [Pa]
    real*8, target,   intent(in)            :: he(km+1,nobs)    ! height above sea-level  [m]
    real*8, target,   intent(in)            :: te(km+1,nobs)    ! temperature at layer edges [K]

    real*8, target,   intent(in)            :: kernel_wt(nkernel,nch,nobs)    ! kernel weights (/fiso,fgeo,fvol/)
    real*8, target,   intent(in)            :: RTLSparam(nparam,nch,nobs)     ! Li-Sparse parameters 
                                                                              ! param1 = crown relative height (h/b)
                                                                              ! param2 = shape parameter (b/r)
    real*8, target,   intent(in)            :: BPDFparam(nparam_bpdf,nch,nobs)     ! BPRDF parameters 
                                                                              ! param1 = Refractive index of water (1.5)
                                                                              ! param2 = NDVI , must be <=1 or >=1
                                                                              ! param3 = Scaling factor (C, from Maignan 2009 RSE Table 1)
                         
    real*8, target,   intent(in)            :: solar_zenith(nobs)  
    real*8, target,   intent(in)            :: relat_azymuth(nobs) 
    real*8, target,   intent(in)            :: sensor_zenith(nobs) 

    real*8,           intent(in)  :: flux_factor(nch,nobs) ! solar flux (F0)
    
    integer,          intent(in)            :: verbose

  ! !OUTPUT PARAMETERS:
    real*8,           intent(out)           :: radiance_VL_SURF(nobs,nch)     ! TOA normalized radiance from VLIDORT using surface module
    real*8,           intent(out)           :: reflectance_VL_SURF(nobs, nch) ! TOA reflectance from VLIDORT using surface module
    integer,          intent(out)           :: rc                             ! return code

    real*8,           intent(out)           :: BR(nobs,nch)                   ! bidirectional reflectance
    real*8,           intent(out)           :: BR_Q(nobs,nch)                 ! bidirectional reflectance Q 
    real*8,           intent(out)           :: BR_U(nobs,nch)                 ! bidirectional reflectance U     
    real*8,           intent(out)           :: Q(nobs, nch)                   ! Stokes parameter Q
    real*8,           intent(out)           :: U(nobs, nch)                   ! Stokes parameter U   

  ! !OPTIONAL PARAMETERS
    logical, optional, intent(in) :: DO_BOA
  !                               ---
    
    integer             :: i,j,n,p,ier

    
    type(VLIDORT_scat) :: SCAT
    type(VLIDORT_output_vector)  :: output

    rc = 0
    ier = 0
    if (present(DO_BOA)) SCAT%DO_BOA = DO_BOA
    SCAT%Surface%Base%NSTREAMS = nstreams
    SCAT%Surface%Base%DO_PLANE_PARALLEL = plane_parallel
    SCAT%Surface%Base%NGREEK_MOMENTS_INPUT = nMom
    call VLIDORT_Init( SCAT%Surface%Base, km, rc, SCAT%DO_BOA)
    if ( rc /= 0 ) return

    SCAT%nMom    = nMom
    SCAT%nPol    = nPol
    SCAT%NSTOKES = 3

    if ( SCAT%NSTOKES  .GT. MAXSTOKES  )   return


    do j = 1, nobs
       
      ! Make sure angles are available
      ! -----------------------------------------
      if ( IS_MISSING(solar_zenith(j),MISSING)  .OR. & 
           IS_MISSING(sensor_zenith(j),MISSING) .OR. &
           IS_MISSING(relat_azymuth(j),MISSING)  )  then

        radiance_VL_SURF(j,:) = MISSING
        reflectance_VL_SURF(j,:) = MISSING
        
        cycle

      end if

      SCAT%pe => pe(:,j)
      SCAT%ze => he(:,j)
      SCAT%te => te(:,j) 

      do i = 1, nch
         ! set solar flux
          SCAT%Surface%Base%VIO%VLIDORT_FixIn%SunRays%TS_FLUX_FACTOR = flux_factor(i,j)

        ! Make sure kernel weights and parameters are defined
        
        do n = 1, nkernel
          
          if (IS_MISSING(kernel_wt(n,i,j),MISSING)) then
            radiance_VL_SURF(j,i) = MISSING
            reflectance_VL_SURF(j,i) = MISSING
            cycle
          end if
        end do
        
        do p = 1, nparam
          if (IS_MISSING(RTLSparam(p,i,j),MISSING)) then
            radiance_VL_SURF(j,i) = MISSING
            reflectance_VL_SURF(j,i) = MISSING
            cycle
          end if
        end do

        do p = 1, nparam_bpdf
          if (IS_MISSING(BPDFparam(p,i,j),MISSING)) then
            radiance_VL_SURF(j,i) = MISSING
            reflectance_VL_SURF(j,i) = MISSING
            cycle
          end if
        end do        

        if ( verbose > 0 ) then
          print*, 'DO MODIS BRDF + BPDF'
        end if

        scalar = .false.
        call VLIDORT_LANDMODIS_BPDF(SCAT%Surface,solar_zenith(j),&
                               sensor_zenith(j),relat_azymuth(j),&
                               kernel_wt(1,i,j),kernel_wt(2,i,j),kernel_wt(3,i,j),&
                               reshape(RTLSparam(:,i,j),(/nparam/)),reshape(BPDFparam(:,i,j),(/nparam_bpdf/)),&
                               scalar,rc)

        SCAT%wavelength = channels(i) 
        SCAT%rot => ROT(:,j,i)
        SCAT%depol_ratio => depol(i)      
        SCAT%alpha => alpha(:,j,i) 
        SCAT%tau => tau(:,i,j)
        SCAT%ssa => ssa(:,i,j)
        SCAT%pmom => pmom(:,i,j,:,:)

        BR(j,i) = SCAT%Surface%Base%VIO%VBRDF_Sup_Out%BS_DBOUNCE_BRDFUNC(1,1,1,1)
        BR_Q(j,i) = SCAT%Surface%Base%VIO%VBRDF_Sup_Out%BS_DBOUNCE_BRDFUNC(2,1,1,1) 
        BR_U(j,i) = SCAT%Surface%Base%VIO%VBRDF_Sup_Out%BS_DBOUNCE_BRDFUNC(3,1,1,1)     

        ! Check to make sure kernel combination is not unphysical
        if (BR(j,i) < 0) then
          radiance_VL_SURF(j,i)    = -500
          reflectance_VL_SURF(j,i) = -500
          BR(j,i)    = -500
          BR_Q(j,i)    = -500
          BR_U(j,i)    = -500
          Q(j,i)     = -500
          U(j,i)     = -500
          cycle
        end if

        call VLIDORT_Run_Vector (SCAT, output, ier)

        if (SCAT%DO_BOA) then
          radiance_VL_SURF(j,i)    = output%BOA_radiance
          reflectance_VL_SURF(j,i) = output%BOA_reflectance
          Q(j,i)                   = output%BOA_Q
          U(j,i)                   = output%BOA_U                          
        else
          radiance_VL_SURF(j,i)    = output%radiance
          reflectance_VL_SURF(j,i) = output%reflectance
          Q(j,i)                   = output%Q
          U(j,i)                   = output%U                
        end if
        
        if ( ier /= 0 ) then
          radiance_VL_SURF(j,i) = MISSING
          reflectance_VL_SURF(j,i) = MISSING               
          cycle
        end if

      end do ! end loop over channels
       
      if ( verbose > 0 ) then
        if ( mod(j-1,1000) == 0 ) then
          print *, '<> VLIDORT Vector: ', nint(j*100./nobs), '%'
        end if
      end if

    end do ! Loop over obs

  end subroutine VLIDORT_Vector_LandMODIS_BPDF_cloud

  subroutine VLIDORT_Vector_BPDF (km, nch, nobs, channels, nstreams, plane_parallel, nMom,  &
                     nPol, ROT, depol, tau, ssa, pmom, pe, he, te, BPDFparam, &
                     solar_zenith, relat_azymuth, sensor_zenith, &
                     MISSING,verbose, radiance_VL_SURF,reflectance_VL_SURF, BR, Q, U, BR_Q, BR_U, rc, DO_BOA)
  !
  ! Place holder.
  !
     use VLIDORT_ScatMod
   
     implicit NONE
 
    logical                                 :: scalar
    integer, parameter                      :: nparam  = 3   ! number of kernel parameters

  ! !INPUT PARAMETERS:

    integer,          intent(in)            :: km    ! number of vertical levels 
    integer,          intent(in)            :: nch   ! number of channels
    integer,          intent(in)            :: nobs  ! number of observations

    logical,          intent(in)            :: plane_parallel ! do plane parallel flag

    integer, target,  intent(in)            :: nMom  ! number of phase function moments 
    integer, target,  intent(in)            :: nPol  ! number of scattering matrix components                               
    integer,          intent(in)            :: nstreams  ! number of half space streams
                    
    real*8, target,   intent(in)            :: channels(nch)    ! wavelengths [nm]

!                                                     ! --- Rayleigh Parameters ---
    real*8, target,   intent(in)            :: ROT(km,nobs,nch) ! rayleigh optical thickness
    real*8, target,   intent(in)            :: depol(nch)       ! rayleigh depolarization ratio used in phase matrix

  !                                                   ! --- Aerosol Optical Properties ---
    real*8, target,   intent(in)            :: tau(km,nch,nobs) ! aerosol optical depth
    real*8, target,   intent(in)            :: ssa(km,nch,nobs) ! single scattering albedo    
    real*8, target,   intent(in)            :: pmom(km,nch,nobs,nMom,nPol) !components of the scat phase matrix

    real*8, target,   intent(in)            :: MISSING          ! MISSING VALUE
    real*8, target,   intent(in)            :: pe(km+1,nobs)    ! pressure at layer edges [Pa]
    real*8, target,   intent(in)            :: he(km+1,nobs)    ! height above sea-level  [m]
    real*8, target,   intent(in)            :: te(km+1,nobs)    ! temperature at layer edges [K]

    real*8, target,   intent(in)            :: BPDFparam(nparam,nch,nobs)     ! BPRDF parameters 
                                                                              ! param1 = Refractive index of water (1.5)
                                                                              ! param2 = NDVI , must be <=1 or >=1
                                                                              ! param3 = Scaling factor (C, from Maignan 2009 RSE Table 1)
                         
    real*8, target,   intent(in)            :: solar_zenith(nobs)  
    real*8, target,   intent(in)            :: relat_azymuth(nobs) 
    real*8, target,   intent(in)            :: sensor_zenith(nobs) 

    integer,          intent(in)            :: verbose

  ! !OUTPUT PARAMETERS:
    real*8,           intent(out)           :: radiance_VL_SURF(nobs,nch)     ! TOA normalized radiance from VLIDORT using surface module
    real*8,           intent(out)           :: reflectance_VL_SURF(nobs, nch) ! TOA reflectance from VLIDORT using surface module
    integer,          intent(out)           :: rc                             ! return code

    real*8,           intent(out)           :: BR(nobs,nch)                   ! bidirectional reflectance 
    real*8,           intent(out)           :: BR_Q(nobs,nch)                   ! bidirectional reflectance 
    real*8,           intent(out)           :: BR_U(nobs,nch)                   ! bidirectional reflectance 
    real*8,           intent(out)           :: Q(nobs, nch)                   ! Stokes parameter Q
    real*8,           intent(out)           :: U(nobs, nch)                   ! Stokes parameter U   

  ! !OPTIONAL PARAMETERS
    logical, optional, intent(in) :: DO_BOA
  !                               ---
    
    integer             :: i,j,n,p,ier

    
    type(VLIDORT_scat) :: SCAT
    type(VLIDORT_output_vector)  :: output

    rc = 0
    ier = 0
    if (present(DO_BOA)) SCAT%DO_BOA = DO_BOA
    SCAT%Surface%Base%NSTREAMS = nstreams
    SCAT%Surface%Base%DO_PLANE_PARALLEL = plane_parallel
    SCAT%Surface%Base%NGREEK_MOMENTS_INPUT = nMom
    call VLIDORT_Init( SCAT%Surface%Base, km, rc, SCAT%DO_BOA)
    if ( rc /= 0 ) return

    SCAT%nMom    = nMom
    SCAT%nPol    = nPol
    SCAT%NSTOKES = 3

    if ( SCAT%NSTOKES  .GT. MAXSTOKES  )   return


    do j = 1, nobs
       
      ! Make sure angles are available
      ! -----------------------------------------
      if ( IS_MISSING(solar_zenith(j),MISSING)  .OR. & 
           IS_MISSING(sensor_zenith(j),MISSING) .OR. &
           IS_MISSING(relat_azymuth(j),MISSING)  )  then

        radiance_VL_SURF(j,:) = MISSING
        reflectance_VL_SURF(j,:) = MISSING
        
        cycle

      end if

      SCAT%pe => pe(:,j)
      SCAT%ze => he(:,j)
      SCAT%te => te(:,j) 

      do i = 1, nch
        ! Make sure kernel weights and parameters are defined
                
        do p = 1, nparam
          if (IS_MISSING(BPDFparam(p,i,j),MISSING)) then
            radiance_VL_SURF(j,i) = MISSING
            reflectance_VL_SURF(j,i) = MISSING
            cycle
          end if

        end do

        if ( verbose > 0 ) then
          print*, 'DO BPDF'
        end if

        scalar = .false.
        call VLIDORT_BPDF(SCAT%Surface,solar_zenith(j),&
                               sensor_zenith(j),relat_azymuth(j),&
                               reshape(BPDFparam(:,i,j),(/nparam/)),&
                               scalar,rc)

        SCAT%wavelength = channels(i)   
        SCAT%rot => ROT(:,j,i)  
        SCAT%depol_ratio => depol(i)   
        SCAT%tau => tau(:,i,j)
        SCAT%ssa => ssa(:,i,j)
        SCAT%pmom => pmom(:,i,j,:,:)

        BR(j,i) = SCAT%Surface%Base%VIO%VBRDF_Sup_Out%BS_DBOUNCE_BRDFUNC(1,1,1,1)
        BR_Q(j,i) = SCAT%Surface%Base%VIO%VBRDF_Sup_Out%BS_DBOUNCE_BRDFUNC(2,1,1,1) 
        BR_U(j,i) = SCAT%Surface%Base%VIO%VBRDF_Sup_Out%BS_DBOUNCE_BRDFUNC(3,1,1,1)     

        ! Check to make sure kernel combination is not unphysical
        if (BR(j,i) < 0) then
          radiance_VL_SURF(j,i)    = -500
          reflectance_VL_SURF(j,i) = -500
          BR(j,i)    = -500
          BR_Q(j,i)  = -500
          BR_U(j,i)  = -500
          Q(j,i)     = -500
          U(j,i)     = -500
          cycle
        end if

        call VLIDORT_Run_Vector (SCAT, output, ier)

        if (SCAT%DO_BOA) then
          radiance_VL_SURF(j,i)    = output%BOA_radiance
          reflectance_VL_SURF(j,i) = output%BOA_reflectance
          Q(j,i)                   = output%BOA_Q
          U(j,i)                   = output%BOA_U                          
        else
          radiance_VL_SURF(j,i)    = output%radiance
          reflectance_VL_SURF(j,i) = output%reflectance
          Q(j,i)                   = output%Q
          U(j,i)                   = output%U                
        end if
        
        if ( ier /= 0 ) then
          radiance_VL_SURF(j,i) = MISSING
          reflectance_VL_SURF(j,i) = MISSING               
          cycle
        end if

      end do ! end loop over channels
       
      if ( verbose > 0 ) then
        if ( mod(j-1,1000) == 0 ) then
          print *, '<> VLIDORT Vector: ', nint(j*100./nobs), '%'
        end if
      end if

    end do ! Loop over obs

  end subroutine VLIDORT_Vector_BPDF

end module VLIDORT_BRDF_MODIS_BPDF
