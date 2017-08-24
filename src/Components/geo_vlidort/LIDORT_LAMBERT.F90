!
!  Simple f77 wrapper for the Python interface to LIDORT for OMI aerosol
!  channels.
!
!.............................................................................

subroutine LIDORT_Scalar_Lambert (km, nch, nobs,channels, nMom, &
                   nPol, tau, ssa, g, pmom, pe, he, te, albedo,            &
                   solar_zenith, relat_azymuth, sensor_zenith, &
                   MISSING,verbose,radiance_L,reflectance_L, ROT, rc)
!
! Uses LIDORT in scalar mode to compute OMI aerosol TOA radiances.
!
  use LIDORT_ScatMod

  implicit NONE

! !INPUT PARAMETERS:

  integer,          intent(in)  :: km    ! number of levels on file
  integer,          intent(in)  :: nch   ! number of channels
  integer,          intent(in)  :: nobs  ! number of observations

  integer, target,  intent(in)  :: nMom             ! number of phase function moments     
  integer, target,  intent(in)  :: nPol  ! number of scattering matrix components                               
                                                        
  real*8, target,   intent(in)  :: channels(nch)    ! wavelengths [nm]

!                                                   ! --- Mie Parameters ---
  real*8, target,   intent(in)  :: tau(km,nch,nobs) ! aerosol optical depth
  real*8, target,   intent(in)  :: ssa(km,nch,nobs) ! single scattering albedo
  real*8, target,   intent(in)  :: g(km,nch,nobs)   ! asymmetry factor
  real*8, target,   intent(in)  :: pmom(km,nch,nobs,nMom,nPol) !components of the scat phase matrix

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

  real*8,           intent(out) :: radiance_L(nobs,nch)       ! TOA normalized radiance from LIDORT
  integer,          intent(out) :: rc                          ! return code
  real*8,           intent(out) :: reflectance_L(nobs, nch)   ! TOA reflectance from LIDORT
  real*8,           intent(out) :: ROT(km,nobs,nch)                 ! rayleigh optical thickness 
 
!                               ---  
  integer             :: i,j, ier
 
  type(LIDORT_scat)            :: SCAT
  type(LIDORT_output_scalar)   :: output  

#define IS_MISSING(x) (abs(x/MISSING-1)<0.001)
  rc = 0
  ier = 0

  call LIDORT_Init( SCAT%Surface%Base, km, rc)
  if ( rc /= 0 ) then
    write(*,*) 'LIDORT_Init returning with error'
    return
  endif

  SCAT%nMom    = nMom

  do j = 1, nobs

     ! Make sure albedo and angles are available
     ! -----------------------------------------
     if ( IS_MISSING(solar_zenith(j))  .OR. & 
          IS_MISSING(sensor_zenith(j)) .OR. &
          IS_MISSING(relat_azymuth(j))  )  then

        radiance_L(j,:) = MISSING
        reflectance_L(j,:) = MISSING
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
              radiance_L(j,i) = MISSING
              reflectance_L(j,i) = MISSING
              cycle
           end if
           
           call LIDORT_SurfaceLamb(SCAT%Surface,albedo(j,i),solar_zenith (j),sensor_zenith(j),&
                               relat_azymuth(j))
           
           SCAT%wavelength = channels(i)
           SCAT%tau => tau(:,i,j)
           SCAT%ssa => ssa(:,i,j)
           SCAT%g => g(:,i,j)
           SCAT%pmom => pmom(:,i,j,:,:)
         
           call LIDORT_Run_Scalar (SCAT, output, ier)

           radiance_L(j,i)    = output%radiance
           reflectance_L(j,i) = output%reflectance
           ROT(:,j,i) = SCAT%rot              


           if ( ier /= 0 ) then
              radiance_L(j,i) = MISSING
              reflectance_L(j,i) = MISSING
              cycle
           end if

        end do ! end loop over channels
     
        if ( verbose > 0 ) then
           if ( mod(j-1,1000) == 0 ) then
              print *, '<> LIDORT Scalar: ', nint(j*100./nobs), '%'
           end if
        end if

  end do ! Loop over obs

end subroutine LIDORT_Scalar_Lambert

!.............................................................................

subroutine LIDORT_AI_Scalar (km, nch, nobs, channels,        &
                   tau, ssa, g, pe, he, te, albedo,               &
                   solar_zenith, relat_azymuth, sensor_zenith, &
                   MISSING,radiance, AI, verbose,radiance_L,          &
                   AI_VL, refl, rc)
!
! Uses LIDORT in scalar mode to compute OMI aerosol TOA radiances and AEROSOL INDEX.
! Look into LIDORT f90 code when the LIDORT RUN is called if scalar is set to .True.(not the default)
! Better to use LIDORT AI_Vector

  use LIDORT_ScatMod

  implicit NONE

  logical,parameter             :: aerosol = .true.   ! whether or not simulation contains aerosol
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
                       
  real*8, target,   intent(in)  :: solar_zenith(nobs)  
  real*8, target,   intent(in)  :: relat_azymuth(nobs) 
  real*8, target,   intent(in)  :: sensor_zenith(nobs) 
  real*8, target,   intent(in)  :: MISSING 

  real*8, target,   intent(in)  :: albedo(nobs,nch)       ! surface albedo
  real*8,           intent(in)  :: radiance(nobs,nch)     ! radiance observed OMI
!  real*8,           intent(in)  :: reflectivity(nobs,nch) ! reflectivity fichier OMI
 
  
  integer,          intent(in)  :: verbose
  real*8,           intent(in)  :: AI(nobs)                    ! Aerosol Index OMI file
! !OUTPUT PARAMETERS:

  real*8,           intent(out) :: radiance_L(nobs,nch)       ! TOA radiance from LIDORT
  integer,          intent(out) :: rc                          ! return code
  real*8,           intent(out) :: AI_VL(nobs)                 ! Aerosol Index model LIDORT
  real*8,           intent(out) :: refl(nobs)                  ! surface reflectivity R* at 388nm
  
!                               ---

  integer             :: n, i,j, k, ier
  real*8              :: spher_alb(nobs,nch)      ! spherical albedo
  real*8              :: trans(nobs, nch)         ! transmission
  real*8              :: reflectivity_388nm       ! Reflectivity lidort 388 nm
  real*8              :: reflectivity_354nm       ! Reflectivity lidort 354 nm
  real*8              :: reflectivity_388nm_OMI   ! Reflectivity OMI file
  real*8              :: reflectivity_354nm_OMI   ! Reflectivity OMI file
  real*8              :: reflectivity_388nm_adj   ! Reflectivity OMI file
  real*8              :: delta_R_model            ! Difference R
  real*8              :: delta_R_OMI              ! Difference R
 
  real*8              :: ROT(km,nobs,nch)            ! rayleigh optical thickness 
  real*8              :: g_tot(nobs,nch)          ! total asymetry factor (column - one profile) 
  real*8              :: sfc_albedo_354nm         ! surface albedo 354 nm OMI file (climato TOMS)
  real*8              :: sfc_albedo_388nm         ! surface albedo 388 nm OMI file (climato TOMS)
  real*8              :: raymodel(nobs)           ! radiance Ray 354 nm from LIDORT
  real*8              :: ray_OMI                  ! radiance Ray 354 nm OMI file
  real*8              :: reflectance_L(nobs,nch)    ! TOA reflectance from LIDORT
  
  type(LIDORT_scat) :: SCAT

#define IS_MISSING(x) (abs(x/MISSING-1)<0.001)

  sfc_albedo_354nm = MISSING
  ier = 0
  rc = 0
  
  call LIDORT_Init(  SCAT%Surface%Base, km, rc)
  if ( rc /= 0 ) return
  
  do j = 1, nobs
      if ( IS_MISSING(solar_zenith(j))  .OR. & 
          IS_MISSING(sensor_zenith(j)) .OR. &
          IS_MISSING(relat_azymuth(j))  )  then

        radiance_L(j,:) = MISSING
        refl(j) = MISSING
        AI_VL(j) = MISSING
        cycle

      end if
         
     SCAT%pe => pe(:,j)
     SCAT%ze => he(:,j)
     SCAT%te => te(:,j) 
 
     do i = 1, nch 
       
        ! Make sure albedo is defined
        ! ---------------------------
           if ( IS_MISSING(albedo(j,i)) ) then
              radiance_L(j,i) = MISSING
              refl(j) = MISSING
              AI_VL(j) = MISSING
              cycle
           end if
           
        call LIDORT_SurfaceLamb(SCAT%Surface,albedo(j,i),solar_zenith (j),sensor_zenith(j),&
                               relat_azymuth(j))

        SCAT%wavelength = channels(i)        
        SCAT%tau => tau(:,i,j)
        SCAT%ssa => ssa(:,i,j)
        SCAT%g => g(:,i,j)

        call LIDORT_Run (SCAT, radiance_L(j,i),reflectance_L(j,i), &
                          ROT(:,j,i), aerosol, ier)

        if ( ier /= 0 ) then
                 radiance_L(j,i) = MISSING
                 AI_VL(j) = MISSING
                 refl(j) = MISSING           
                 cycle
        end if 
               
        if (abs(SCAT%wavelength - 354.00) < 0.01)  sfc_albedo_354nm = albedo(j,i) ! surface albedo  at 354 nm
        
        if (abs(SCAT%wavelength - 388.00) < 0.01) then           
            if ((.not. IS_MISSING(sfc_albedo_354nm)) .and.&
               (.not. IS_MISSING(radiance_L(j,i)) ))  then ! calcul of S, T and R* at 388 nm
                
                call LIDORT_LER(SCAT, radiance_L(j,i), refl(j), rc)
                sfc_albedo_388nm = albedo(j,i)
                reflectivity_388nm_adj = refl(j)-(sfc_albedo_388nm -sfc_albedo_354nm) 
            else              
                reflectivity_388nm_adj = MISSING
            end if       
        end if  ! end calcul of reflectivity 388nm
        
         
     end do ! end loop over nch
      
!     do calcul AI only at 354 nm using refl_388nm (if not = MISSING) -> i = 1  
        
        if (IS_MISSING(reflectivity_388nm_adj)) then
            AI_VL(j) = MISSING        
        else      
        SCAT%wavelength = channels(1)                  
        SCAT%tau => tau(:,1,j)
        SCAT%ssa => ssa(:,1,j)
        SCAT%g => g(:,1,j) 
 
        call LIDORT_AI(SCAT, radiance_L(j,1),reflectivity_388nm_adj, AI_VL(j), rc)    
        
        end if 
  end do

end subroutine LIDORT_AI_Scalar

!.............................................................................
subroutine LIDORT_Scalar_Lambert_Cloud (km, nch, nobs,channels, nMom, &
                   nPol, tau, ssa, g, pmom, tauI, ssaI, gI, pmomI, tauL, ssaL, gL, pmomL, &
                   pe, he, te, albedo,            &
                   solar_zenith, relat_azymuth, sensor_zenith, &
                   MISSING,verbose,radiance_L,reflectance_L, ROT, rc)
!
! Uses LIDORT in scalar mode to compute OMI aerosol TOA radiances.
!
  use LIDORT_ScatMod

  implicit NONE

! !INPUT PARAMETERS:

  integer,          intent(in)  :: km    ! number of levels on file
  integer,          intent(in)  :: nch   ! number of channels
  integer,          intent(in)  :: nobs  ! number of observations

  integer, target,  intent(in)  :: nMom             ! number of phase function moments     
  integer, target,  intent(in)  :: nPol  ! number of scattering matrix components                               
                                                        
  real*8, target,   intent(in)  :: channels(nch)    ! wavelengths [nm]

!                                                   ! --- Mie Parameters ---
  real*8, target,   intent(in)  :: tau(km,nch,nobs) ! aerosol optical depth
  real*8, target,   intent(in)  :: ssa(km,nch,nobs) ! single scattering albedo
  real*8, target,   intent(in)  :: g(km,nch,nobs)   ! asymmetry factor
  real*8, target,   intent(in)  :: pmom(km,nch,nobs,nMom,nPol) !components of the scat phase matrix

  real*8, target,   intent(in)  :: tauI(km,nch,nobs) ! ice cloud optical depth
  real*8, target,   intent(in)  :: ssaI(km,nch,nobs) ! ice cloud single scattering albedo
  real*8, target,   intent(in)  :: gI(km,nch,nobs)   ! ice cloud asymmetry factor
  real*8, target,   intent(in)  :: pmomI(km,nch,nobs,nMom,nPol) !components of the scat phase matrix

  real*8, target,   intent(in)  :: tauL(km,nch,nobs) ! liquid cloud optical depth
  real*8, target,   intent(in)  :: ssaL(km,nch,nobs) ! liquid cloud single scattering albedo
  real*8, target,   intent(in)  :: gL(km,nch,nobs)   ! liquid cloud asymmetry factor
  real*8, target,   intent(in)  :: pmomL(km,nch,nobs,nMom,nPol) !components of the scat phase matrix

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

  real*8,           intent(out) :: radiance_L(nobs,nch)       ! TOA normalized radiance from LIDORT
  integer,          intent(out) :: rc                          ! return code
  real*8,           intent(out) :: reflectance_L(nobs, nch)   ! TOA reflectance from LIDORT
  real*8,           intent(out) :: ROT(km,nobs,nch)                 ! rayleigh optical thickness 
 
!                               ---  
  integer             :: i,j, ier
 
  type(LIDORT_scat)            :: SCAT
  type(LIDORT_output_scalar)   :: output  

#define IS_MISSING(x) (abs(x/MISSING-1)<0.001)
  rc = 0
  ier = 0

  call LIDORT_Init( SCAT%Surface%Base, km, rc)
  if ( rc /= 0 ) then
    write(*,*) 'LIDORT_Init returning with error'
    return
  endif

  SCAT%nMom    = nMom

  do j = 1, nobs

     ! Make sure albedo and angles are available
     ! -----------------------------------------
     if ( IS_MISSING(solar_zenith(j))  .OR. & 
          IS_MISSING(sensor_zenith(j)) .OR. &
          IS_MISSING(relat_azymuth(j))  )  then

        radiance_L(j,:) = MISSING
        reflectance_L(j,:) = MISSING
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
              radiance_L(j,i) = MISSING
              reflectance_L(j,i) = MISSING
              cycle
           end if
           
           call LIDORT_SurfaceLamb(SCAT%Surface,albedo(j,i),solar_zenith (j),sensor_zenith(j),&
                               relat_azymuth(j))
           
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

           radiance_L(j,i)    = output%radiance
           reflectance_L(j,i) = output%reflectance
           ROT(:,j,i) = SCAT%rot              


           if ( ier /= 0 ) then
              radiance_L(j,i) = MISSING
              reflectance_L(j,i) = MISSING
              cycle
           end if

        end do ! end loop over channels
     
        if ( verbose > 0 ) then
           if ( mod(j-1,1000) == 0 ) then
              print *, '<> LIDORT Scalar: ', nint(j*100./nobs), '%'
           end if
        end if

  end do ! Loop over obs

end subroutine LIDORT_Scalar_Lambert_Cloud

!.............................................................................