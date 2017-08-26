module LIDORT_BRDF_MODIS   
!
!  Simple f77 wrapper for the Python interface to LIDORT for OMI aerosol
!  channels.
!
!.............................................................................
  implicit NONE

  PUBLIC LIDORT_Scalar_LandMODIS
  PUBLIC LIDORT_Scalar_LandMODIS_Cloud

  contains

  logical function IS_MISSING(x,MISSING)
    real*8, intent(in)     :: x
    real*8, intent(in)     :: MISSING

    IS_MISSING = abs(x/MISSING-1)<0.001
    return
  end function IS_MISSING

  subroutine LIDORT_Scalar_LandMODIS (km, nch, nobs,channels, nMom,       &
                     nPol, tau, ssa, g, pmom, pe, he, te, kernel_wt, param, &
                     solar_zenith, relat_azymuth, sensor_zenith, &
                     MISSING,verbose,radiance_L_SURF,reflectance_L_SURF, ROT, BR, rc )
  !
  ! Uses LIDORT in scalar mode to compute OMI aerosol TOA radiances.
  !
    use LIDORT_ScatMod


    implicit NONE

    logical,parameter             :: aerosol = .true.   ! whether or not simulation contains aerosol
    integer, parameter            :: nkernel = 3
    integer, parameter            :: nparam  = 2
  ! !INPUT PARAMETERS:
    integer,          intent(in)  :: km    ! number of vertical levels
    integer,          intent(in)  :: nch   ! number of channels
    integer,          intent(in)  :: nobs  ! number of observations
                                        
    integer, target,  intent(in)  :: nMom             ! number of phase function moments     
    integer, target,  intent(in)  :: nPol  ! number of scattering matrix components                               

    real*8, target,   intent(in)  :: channels(nch)    ! wavelengths [nm]

  !                                                   ! --- Aerosol Optical Properties ---
    real*8, target,   intent(in)  :: tau(km,nch,nobs) ! aerosol optical depth
    real*8, target,   intent(in)  :: ssa(km,nch,nobs) ! single scattering albedo
    real*8, target,   intent(in)  :: g(km,nch,nobs)   ! asymmetry factor
    real*8, target,   intent(in)  :: pmom(km,nch,nobs,nMom,nPol) !components of the scat phase matrix


    real*8, target,   intent(in)  :: pe(km+1,nobs)    ! pressure at layer edges [Pa]
    real*8, target,   intent(in)  :: he(km+1,nobs)    ! height above sea-level  [m]
    real*8, target,   intent(in)  :: te(km+1,nobs)    ! temperature at layer edges [K]

    real*8, target,   intent(in)  :: kernel_wt(nkernel,nch,nobs)   ! kernel weights (/fiso,fgeo,fvol/)
    real*8, target,   intent(in)  :: param(nparam,nch,nobs)        ! Li-Sparse parameters 
                                                                   ! param1 = crown relative height (h/b)
                                                                   ! param2 = shape parameter (b/r)
    
    real*8, target,   intent(in)  :: MISSING          ! MISSING VALUE                                      
    real*8, target,   intent(in)  :: solar_zenith(nobs)  
    real*8, target,   intent(in)  :: relat_azymuth(nobs) 
    real*8, target,   intent(in)  :: sensor_zenith(nobs) 
    
    integer,          intent(in)            :: verbose

  ! !OUTPUT PARAMETERS:

    real*8,           intent(out) :: radiance_L_SURF(nobs,nch)       ! TOA normalized radiance from LIDORT using surface module
    real*8,           intent(out) :: reflectance_L_SURF(nobs, nch)   ! TOA reflectance from LIDORT using surface module
    integer,          intent(out) :: rc                               ! return code
    real*8,           intent(out) :: ROT(km,nobs,nch)                 ! rayleigh optical thickness  
    real*8,           intent(out) :: BR(nobs,nch)                     ! bi-directional reflectance 

    integer             :: i,j,n,p,ier

  
    type(LIDORT_scat)           :: SCAT
    type(LIDORT_output_scalar)  :: output

    rc = 0
    ier = 0

    call LIDORT_Init( SCAT%Surface%Base, km, rc)
    if ( rc /= 0 ) return

    SCAT%nMom    = nMom

    do j = 1,nobs

       ! Make sure albedo and angles are available
       ! -----------------------------------------
      if ( IS_MISSING(solar_zenith(j),MISSING)  .OR. & 
           IS_MISSING(sensor_zenith(j),MISSING) .OR. &
           IS_MISSING(relat_azymuth(j),MISSING)  )  then

        radiance_L_SURF(j,:) = MISSING
        reflectance_L_SURF(j,:) = MISSING
        cycle
      end if
        
      SCAT%pe => pe(:,j)
      SCAT%ze => he(:,j)
      SCAT%te => te(:,j) 
   
      ! Loop over channels
      ! ------------------
      do i = 1, nch 
        ! Make sure kernel weights and parameters are defined
        do n = 1, nkernel
          if (IS_MISSING(kernel_wt(n,i,j),MISSING)) then
            radiance_L_SURF(j,i) = MISSING
            reflectance_L_SURF(j,i) = MISSING
            cycle
          end if
        end do
        do p = 1, nparam
          if (IS_MISSING(param(p,i,j),MISSING)) then
            radiance_L_SURF(j,i) = MISSING
            reflectance_L_SURF(j,i) = MISSING
            cycle
          end if
        end do

        if ( verbose > 0 ) then
          print*, 'DO MODIS BRDF'
        end if
        call LIDORT_LANDMODIS(SCAT%Surface,solar_zenith(j),&
                               sensor_zenith(j),relat_azymuth(j),&
                               kernel_wt(1,i,j),kernel_wt(2,i,j),kernel_wt(3,i,j),&
                               reshape(param(:,i,j),(/nparam/)),rc)

        if ( rc /= 0 ) return

        SCAT%wavelength = channels(i)
        SCAT%tau => tau(:,i,j)
        SCAT%ssa => ssa(:,i,j)
        SCAT%g => g(:,i,j)
        SCAT%pmom => pmom(:,i,j,:,:)
         
        call LIDORT_Run_Scalar (SCAT, output, ier)

        radiance_L_SURF(j,i)    = output%radiance
        reflectance_L_SURF(j,i) = output%reflectance
        BR(j,i) = SCAT%Surface%Base%VIO%BRDF_Sup_Out%BS_DBOUNCE_BRDFUNC(1,1,1)
        ROT(:,j,i) = SCAT%rot   

        if ( verbose > 0 ) then
          print *, 'My radiance land modis',radiance_L_SURF(j,i), reflectance_L_SURF(j,i) 
        end if

        if ( ier /= 0 ) then
          radiance_L_SURF(j,i) = MISSING
          reflectance_L_SURF(j,i) = MISSING
          cycle
        end if

        if (BR(j,i) < 0) then
          radiance_L_SURF(j,i)    = -500
          reflectance_L_SURF(j,i) = -500
          BR(j,i)    = -500
          ROT(:,j,i) = -500
          cycle
        end if        

      end do ! end loop over channels
     
      if ( verbose > 0 ) then
        if ( mod(j-1,1000) == 0 ) then
          print *, '<> LIDORT Scalar: ', nint(j*100./nobs), '%'
        end if
      end if

    end do ! Loop over obs

  end subroutine LIDORT_Scalar_LandMODIS

  !..........................................................................
  subroutine LIDORT_Scalar_LandMODIS_Cloud (km, nch, nobs,channels, nMom,       &
                     nPol, tau, ssa, g, pmom, tauI, ssaI, gI, pmomI, tauL, ssaL, gL, pmomL, &
                     pe, he, te, kernel_wt, param, &
                     solar_zenith, relat_azymuth, sensor_zenith, &
                     MISSING,verbose,radiance_L_SURF,reflectance_L_SURF, ROT, BR, rc )
  !
  ! Uses LIDORT in scalar mode to compute OMI aerosol TOA radiances.
  !
    use LIDORT_ScatMod


    implicit NONE

    logical,parameter             :: aerosol = .true.   ! whether or not simulation contains aerosol
    integer, parameter            :: nkernel = 3
    integer, parameter            :: nparam  = 2
  ! !INPUT PARAMETERS:
    integer,          intent(in)  :: km    ! number of vertical levels
    integer,          intent(in)  :: nch   ! number of channels
    integer,          intent(in)  :: nobs  ! number of observations
                                        
    integer, target,  intent(in)  :: nMom             ! number of phase function moments     
    integer, target,  intent(in)  :: nPol  ! number of scattering matrix components                               

    real*8, target,   intent(in)  :: channels(nch)    ! wavelengths [nm]

  !                                                   ! --- Aerosol Optical Properties ---
    real*8, target,   intent(in)  :: tau(km,nch,nobs) ! aerosol optical depth
    real*8, target,   intent(in)  :: ssa(km,nch,nobs) ! single scattering albedo
    real*8, target,   intent(in)  :: g(km,nch,nobs)   ! asymmetry factor
    real*8, target,   intent(in)  :: pmom(km,nch,nobs,nMom,nPol) !components of the scat phase matrix

    real*8, target,   intent(in)  :: tauI(km,nch,nobs) ! ice cloud optical depth
    real*8, target,   intent(in)  :: ssaI(km,nch,nobs) ! ice cloud single scattering albedo
    real*8, target,   intent(in)  :: gI(km,nch,nobs)   ! ice cloud asymmetry factor
    real*8, target,   intent(in)  :: pmomI(km,nch,nobs,nMom,nPol) !components of the scat phase matrix

    real*8, target,   intent(in)  :: tauL(km,nch,nobs) ! liquid cloud optical depth
    real*8, target,   intent(in)  :: ssaL(km,nch,nobs) ! ice cloud single scattering albedo
    real*8, target,   intent(in)  :: gL(km,nch,nobs)   ! liquid cloud asymmetry factor
    real*8, target,   intent(in)  :: pmomL(km,nch,nobs,nMom,nPol) !components of the scat phase matrix

    real*8, target,   intent(in)  :: pe(km+1,nobs)    ! pressure at layer edges [Pa]
    real*8, target,   intent(in)  :: he(km+1,nobs)    ! height above sea-level  [m]
    real*8, target,   intent(in)  :: te(km+1,nobs)    ! temperature at layer edges [K]

    real*8, target,   intent(in)  :: kernel_wt(nkernel,nch,nobs)   ! kernel weights (/fiso,fgeo,fvol/)
    real*8, target,   intent(in)  :: param(nparam,nch,nobs)        ! Li-Sparse parameters 
                                                                   ! param1 = crown relative height (h/b)
                                                                   ! param2 = shape parameter (b/r)
    
    real*8, target,   intent(in)  :: MISSING          ! MISSING VALUE                                      
    real*8, target,   intent(in)  :: solar_zenith(nobs)  
    real*8, target,   intent(in)  :: relat_azymuth(nobs) 
    real*8, target,   intent(in)  :: sensor_zenith(nobs) 
    
    integer,          intent(in)            :: verbose

  ! !OUTPUT PARAMETERS:

    real*8,           intent(out) :: radiance_L_SURF(nobs,nch)       ! TOA normalized radiance from LIDORT using surface module
    real*8,           intent(out) :: reflectance_L_SURF(nobs, nch)   ! TOA reflectance from LIDORT using surface module
    integer,          intent(out) :: rc                               ! return code
    real*8,           intent(out) :: ROT(km,nobs,nch)                 ! rayleigh optical thickness  
    real*8,           intent(out) :: BR(nobs,nch)                     ! bi-directional reflectance 

    integer             :: i,j,n,p,ier

  
    type(LIDORT_scat)           :: SCAT
    type(LIDORT_output_scalar)  :: output

    rc = 0
    ier = 0

    call LIDORT_Init( SCAT%Surface%Base, km, rc)
    if ( rc /= 0 ) return

    SCAT%nMom    = nMom

    do j = 1,nobs

       ! Make sure albedo and angles are available
       ! -----------------------------------------
      if ( IS_MISSING(solar_zenith(j),MISSING)  .OR. & 
           IS_MISSING(sensor_zenith(j),MISSING) .OR. &
           IS_MISSING(relat_azymuth(j),MISSING)  )  then

        radiance_L_SURF(j,:) = MISSING
        reflectance_L_SURF(j,:) = MISSING
        cycle
      end if
        
      SCAT%pe => pe(:,j)
      SCAT%ze => he(:,j)
      SCAT%te => te(:,j) 
   
      ! Loop over channels
      ! ------------------
      do i = 1, nch 
        ! Make sure kernel weights and parameters are defined
        do n = 1, nkernel
          if (IS_MISSING(kernel_wt(n,i,j),MISSING)) then
            radiance_L_SURF(j,i) = MISSING
            reflectance_L_SURF(j,i) = MISSING
            cycle
          end if
        end do
        do p = 1, nparam
          if (IS_MISSING(param(p,i,j),MISSING)) then
            radiance_L_SURF(j,i) = MISSING
            reflectance_L_SURF(j,i) = MISSING
            cycle
          end if
        end do

        if ( verbose > 0 ) then
          print*, 'DO MODIS BRDF'
        end if
        call LIDORT_LANDMODIS(SCAT%Surface,solar_zenith(j),&
                               sensor_zenith(j),relat_azymuth(j),&
                               kernel_wt(1,i,j),kernel_wt(2,i,j),kernel_wt(3,i,j),&
                               reshape(param(:,i,j),(/nparam/)),rc)

        if ( rc /= 0 ) return

        SCAT%wavelength = channels(i)
        SCAT%tau => tau(:,i,j)
        SCAT%ssa => ssa(:,i,j)
        SCAT%g => g(:,i,j)
        SCAT%pmom => pmom(:,i,j,:,:)
        SCAT%tauI => tauI(:,i,j)
        SCAT%ssaI => ssaI(:,i,j)
        SCAT%gI => gI(:,i,j)
        SCAT%pmomI => pmomI(:,i,j,:,:)
        SCAT%tauL => tauL(:,i,j)
        SCAT%ssaL => ssaL(:,i,j)
        SCAT%gL => gL(:,i,j)
        SCAT%pmomL => pmomL(:,i,j,:,:)
         
        call LIDORT_Run_Scalar (SCAT, output, ier)

        radiance_L_SURF(j,i)    = output%radiance
        reflectance_L_SURF(j,i) = output%reflectance
        BR(j,i) = SCAT%Surface%Base%VIO%BRDF_Sup_Out%BS_DBOUNCE_BRDFUNC(1,1,1)
        ROT(:,j,i) = SCAT%rot   

        if ( verbose > 0 ) then
          print *, 'My radiance land modis',radiance_L_SURF(j,i), reflectance_L_SURF(j,i) 
        end if

        if ( ier /= 0 ) then
          radiance_L_SURF(j,i) = MISSING
          reflectance_L_SURF(j,i) = MISSING
          cycle
        end if

        if (BR(j,i) < 0) then
          radiance_L_SURF(j,i)    = -500
          reflectance_L_SURF(j,i) = -500
          BR(j,i)    = -500
          ROT(:,j,i) = -500
          cycle
        end if        

      end do ! end loop over channels
     
      if ( verbose > 0 ) then
        if ( mod(j-1,1000) == 0 ) then
          print *, '<> LIDORT Scalar: ', nint(j*100./nobs), '%'
        end if
      end if

    end do ! Loop over obs

  end subroutine LIDORT_Scalar_LandMODIS_Cloud

end module LIDORT_BRDF_MODIS