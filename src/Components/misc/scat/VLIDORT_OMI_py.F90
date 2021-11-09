!
!  Simple f77 wrapper for the Python interface to VLIDORT for OMI aerosol
!  channels.
!
!.............................................................................

subroutine Scalar (km, nch, nobs,channels,        &
                   tau, ssa, g, pe, he, te, albedo,            &
                   solar_zenith, relat_azymuth, sensor_zenith, &
                   MISSING,verbose,radiance_VL,reflectance_VL, rc)
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
  real*8              :: SSA_VL(nobs,nch)   
  real*8              :: Q(nobs,nch),U(nobs,nch)    
  type(VLIDORT_scat) :: SCAT

#define IS_MISSING(x) (abs(x/MISSING-1)<0.001)
  rc = 0
  ier = 0

  call VLIDORT_Init( SCAT%Surface%Base, km, rc)
  if ( rc /= 0 ) return

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
         
           call VLIDORT_Run (SCAT, radiance_VL(j,i), reflectance_VL(j,i), &
                                 AOT(j,i),SSA_VL(j,i),ier,Q(j,i),U(j,i), .true., .true.)

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

!..........................................................................

subroutine Vector (km, nch, nobs, channels, nMom,  &
                   nPol,tau, ssa, g, pmom, pe, he, te, albedo, &
                   solar_zenith, relat_azymuth, sensor_zenith, &
                   MISSING,verbose, radiance_VL, reflectance_VL,Q,U, rc)
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

                       
  real*8, target,   intent(in)  :: solar_zenith(nobs)  
  real*8, target,   intent(in)  :: relat_azymuth(nobs) 
  real*8, target,   intent(in)  :: sensor_zenith(nobs) 

  real*8, target,   intent(in)  :: albedo(nobs,nch)       ! surface albedo

  integer,          intent(in)  :: verbose

! !OUTPUT PARAMETERS:

  real*8,           intent(out) :: radiance_VL(nobs,nch)       ! TOA radiance from VLIDORT
  integer,          intent(out) :: rc                          ! return code
  real*8,           intent(out) :: reflectance_VL(nobs, nch)   ! TOA reflectance from VLIDORT
  real*8,           intent(out) :: Q(nobs, nch)   ! Q Stokes component
  real*8,           intent(out) :: U(nobs, nch)   ! U Stokes component
!                               ---
  
  integer             :: i,j, ier
  real*8              :: AOT(nobs,nch)            ! total aerosol optical thickness 
  real*8              :: SSA_VL(nobs,nch)   
  
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
        Q(j,i) = MISSING
        U(j,i) = MISSING
        cycle

      end if
     
     SCAT%pe => pe(:,j)
     SCAT%ze => he(:,j)
     SCAT%te => te(:,j) 

     do i = 1, nch
       if ( IS_MISSING(albedo(j,i)) ) then
              radiance_VL(j,i) = MISSING
              reflectance_VL(j,i) = MISSING
              Q(j,i) = MISSING
              U(j,i) = MISSING
              cycle
       end if

       call VLIDORT_SurfaceLamb(SCAT%Surface,albedo(j,i),solar_zenith (j),sensor_zenith(j),&
                               relat_azymuth(j),.false.)

        SCAT%wavelength = channels(i)        
        SCAT%tau => tau(:,i,j)
        SCAT%ssa => ssa(:,i,j)
        SCAT%g => g(:,i,j)
        SCAT%pmom => pmom(:,:,:,i,j)

        call VLIDORT_Run (SCAT, radiance_VL(j,i),reflectance_VL(j,i), AOT(j,i),SSA_VL(j,i),&
                           ier, Q(j,i),U(j,i),.false., .true.)
        if ( ier /= 0 ) then

              if ( solar_zenith (j) > 89.9 ) then ! if SZA > 90
                   radiance_VL(j,i) = 0.
                   reflectance_VL(j,i) = 0.
                   Q(j,i) = 0.
                   U(j,i) = 0.
              else
                   radiance_VL(j,i) = MISSING
                   reflectance_VL(j,i) = MISSING
                   Q(j,i) = MISSING
                   U(j,i) = MISSING
              endif
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
!.............................................................................

subroutine AI_Scalar (km, nch, nobs, channels,        &
                   tau, ssa, g, pe, he, te, albedo,               &
                   solar_zenith, relat_azymuth, sensor_zenith, &
                   MISSING,radiance, AI, verbose,radiance_VL,          &
                   AI_VL, refl, rc)
!
! Uses VLIDORT in scalar mode to compute OMI aerosol TOA radiances and AEROSOL INDEX.
! Look into VLIDORT f90 code when the VLIDORT RUN is called if scalar is set to .True.(not the default)
! Better to use VLIDORT AI_Vector

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

  real*8,           intent(out) :: radiance_VL(nobs,nch)       ! TOA radiance from VLIDORT
  integer,          intent(out) :: rc                          ! return code
  real*8,           intent(out) :: AI_VL(nobs)                 ! Aerosol Index model VLIDORT
  real*8,           intent(out) :: refl(nobs)                  ! surface reflectivity R* at 388nm
  
!                               ---
  
  integer             :: n, i,j, k, ier
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
  real*8              :: Q(nobs,nch), U(nobs,nch)
  real*8              :: SSA_VL(nobs,nch) 
  type(VLIDORT_scat) :: SCAT

#define IS_MISSING(x) (abs(x/MISSING-1)<0.001)

  sfc_albedo_354nm = MISSING
  ier = 0
  rc = 0
  
  call VLIDORT_Init(  SCAT%Surface%Base, km, rc)
  if ( rc /= 0 ) return
  
  do j = 1, nobs
      if ( IS_MISSING(solar_zenith(j))  .OR. & 
          IS_MISSING(sensor_zenith(j)) .OR. &
          IS_MISSING(relat_azymuth(j))  )  then

        radiance_VL(j,:) = MISSING
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
              radiance_VL(j,i) = MISSING
              refl(j) = MISSING
              AI_VL(j) = MISSING
              cycle
           end if
           
        call VLIDORT_SurfaceLamb(SCAT%Surface,albedo(j,i),solar_zenith (j),sensor_zenith(j),&
                               relat_azymuth(j),.true.)

        SCAT%wavelength = channels(i)        
        SCAT%tau => tau(:,i,j)
        SCAT%ssa => ssa(:,i,j)
        SCAT%g => g(:,i,j)
                 
        call VLIDORT_Run (SCAT, radiance_VL(j,i),reflectance_VL(j,i), AOT(j,i),SSA_VL(j,i),&
                          ier, Q(j,i),U(j,i),.true., .true.)

        if ( ier /= 0 ) then
                 radiance_VL(j,i) = MISSING
                 AI_VL(j) = MISSING
                 refl(j) = MISSING           
                 cycle
        end if 
               
        if (abs(SCAT%wavelength - 354.00) < 0.01)  sfc_albedo_354nm = albedo(j,i) ! surface albedo  at 354 nm
        
        if (abs(SCAT%wavelength - 388.00) < 0.01) then           
            if ((.not. IS_MISSING(sfc_albedo_354nm) )   .and. &
                (.not. IS_MISSING(radiance_VL(j,i)) ))  then ! calcul of S, T and R* at 388 nm
                
                call VLIDORT_LER(SCAT, radiance_VL(j,i), refl(j), rc)
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
 
        call VLIDORT_AI(SCAT, radiance_VL(j,1),reflectivity_388nm_adj, AI_VL(j), rc)    
        
        end if 
  end do

end subroutine AI_Scalar
!..........................................................................

subroutine AI_Vector (km, nch, nobs,  channels, nMom,  &
                   nPol,tau, ssa, g, pmom, pe, he, te,oceanler,albedo,     &
                   solar_zenith,relat_azymuth, sensor_zenith,    &
                   qa_grd, &
                   MISSING,verbose, radiance_VL,  &
                   AI_VL,SSA_VL, AOT,refl, rc, &
                   icalc388, trans388, spher388, residue, BRDF)
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


  real*8, target,   intent(in)  :: pe(km+1,nobs)    ! pressure at layer edges [Pa]
  real*8, target,   intent(in)  :: he(km+1,nobs)    ! height above sea-level  [m]
  real*8, target,   intent(in)  :: te(km+1,nobs)    ! temperature at layer edges [K]
                       
  real*8, target,   intent(in)  :: solar_zenith(nobs)  
  real*8, target,   intent(in)  :: relat_azymuth(nobs) 
  real*8, target,   intent(in)  :: sensor_zenith(nobs) 
  integer, target,  intent(in)  :: qa_grd(nobs) 
  real*8, target,   intent(in)  :: MISSING 

  real*8,           intent(in)  :: oceanler(2,16,16,16)
  real*8, target,   intent(in)  :: albedo(nobs,nch)       ! surface albedo
!  real*8,           intent(in)  :: radiance(nobs,nch)     ! radiance observed OMI
!  real*8,           intent(in)  :: reflectivity(nobs,nch)! reflectivity fichier OMI
!  real*8,           intent(in)  :: AI(nobs)               ! Aerosol Index OMI file
  
  integer,          intent(in)  :: verbose

! !OUTPUT PARAMETERS:

  real*8,           intent(out) :: radiance_VL(nobs,nch)       ! TOA radiance from VLIDORT
  integer,          intent(out) :: rc                          ! return code
  real*8,           intent(out) :: AI_VL(nobs)                 ! Aerosol Index model VLIDORT
  real*8,           intent(out) :: refl(nobs)                  ! surface reflectivity R* at 388nm
  real*8,           intent(out)  :: SSA_VL(nobs,nch)  
  real*8,           intent(out)  :: AOT(nobs,nch)            ! total aerosol optical thickness                           

  real*8,           intent(out)  :: icalc388(nobs), trans388(nobs), spher388(nobs), &
                                    residue(nobs), BRDF(nobs,nch)
  
  integer             :: n, i,j, k, imom, ipol, i_354,ier
  real*8              :: spher_alb(nobs,nch)      ! spherical albedo
  real*8              :: trans(nobs, nch)         ! transmission
  real*8              :: reflectivity_388nm       ! Reflectivity Vlidort 388 nm
  real*8              :: reflectivity_354nm       ! Reflectivity Vlidort 354 nm
  real*8              :: reflectivity_388nm_OMI   ! Reflectivity OMI file
  real*8              :: reflectivity_354nm_OMI   ! Reflectivity OMI file
  real*8              :: reflectivity_388nm_adj   ! Reflectivity OMI file
  real*8              :: delta_R_model            ! Difference R
  real*8              :: delta_R_OMI              ! Difference R
 
  real*8              :: g_tot(nobs,nch)            ! total asymetry factor 
  real*8              :: sfc_albedo_354nm         ! surface albedo 354 nm OMI file (climato TOMS)
  real*8              :: sfc_albedo_388nm         ! surface albedo 388 nm OMI file (climato TOMS)
  real*8              :: raymodel(nobs)           ! radiance Ray 354 nm from VLIDORT
  real*8              :: ray_OMI                  ! radiance Ray 354 nm OMI file
  real*8              :: reflectance_VL(nobs,nch)    ! TOA reflectance from VLIDORT
  real*8              :: Q(nobs,nch), U(nobs,nch) 
  real*8              :: icalc354(nobs), trans354(nobs), spher354(nobs)
  real*8              :: cf
  real*8, dimension(2):: oceanlerout
! spectral real part of refractive index at 354, 388, and 470 nm from hitran
  real*8, dimension(3):: mr = (/1.34338,   1.34014,   1.33582/)
  real*8              :: u10m, v10m
  integer             :: qa
 
  type(VLIDORT_scat) :: SCAT

#define IS_MISSING(x) (abs(x/MISSING-1)<0.001)  
  
  rc = 0
  ier = 0

  call VLIDORT_Init(  SCAT%Surface%Base, km, rc)

  if ( rc /= 0 ) return

  SCAT%nMom = nMom
  SCAT%nPol = nPol

  do j = 1, nobs

!    Initialize values to a missing value
     qa = ibits(qa_grd(j),0,4)

     sfc_albedo_354nm = MISSING
     sfc_albedo_388nm = MISSING
     reflectivity_388nm_adj = MISSING
     radiance_VL(j,:) = MISSING
     refl(j) = MISSING
     AI_VL(j) = MISSING
     icalc388(j) = MISSING
     trans388(j) = MISSING
     spher388(j) = MISSING
     icalc354(j) = MISSING
     trans354(j) = MISSING
     spher354(j) = MISSING
     residue(j) = MISSING

!    If OMI not providing needed geometry then cycle to next obs
     if ( IS_MISSING(solar_zenith(j))  .OR. & 
          IS_MISSING(sensor_zenith(j)) .OR. &
          IS_MISSING(relat_azymuth(j))  )  then
        cycle
     end if
          
     SCAT%pe => pe(:,j)
     SCAT%ze => he(:,j)
     SCAT%te => te(:,j) 

     do i = 1, nch

       ! Make sure albedo is defined
           ! ---------------------------
           if ( IS_MISSING(albedo(j,i)) ) then
              radiance_VL(j,i) = MISSING
              cycle
           end if
          
! Introduce Fresnel reflectance over the ocean
          if(qa .eq. 0 .or. qa .ge. 6) then
            u10m = 6.  ! these should come from the model
            v10m = 0.
            call VLIDORT_GissCoxMunk(SCAT%Surface,U10m,V10m,mr(i),solar_zenith(j),&
                                     sensor_zenith(j),relat_azymuth(j),.false.,rc)
            BRDF(j,i) = SCAT%Surface%Base%VIO%VBRDF_Sup_Out%BS_DBOUNCE_BRDFUNC(1,1,1,1)
          else
            call VLIDORT_SurfaceLamb(SCAT%Surface,albedo(j,i),solar_zenith (j),sensor_zenith(j),&
                                     relat_azymuth(j),.false.)
          endif

          SCAT%wavelength = channels(i)          
          SCAT%tau => tau(:,i,j)
          SCAT%ssa => ssa(:,i,j)
          SCAT%g => g(:,i,j)
          SCAT%pmom => pmom(:,:,:,i,j)

          call VLIDORT_Run (SCAT, radiance_VL(j,i),reflectance_VL(j,i), AOT(j,i),SSA_VL(j,i),&
                            ier,Q(j,i),U(j,i), .false., .true.)

          if ( ier /= 0 ) then
             radiance_VL(j,i) = MISSING
             cycle
          end if 

! Calculations for 354 nm
          if (abs(SCAT%wavelength - 354.00) < 0.01) then
               sfc_albedo_354nm = albedo(j,i)  ! surface albedo at 354 nm
!              PRC: For 1.62.1 algorithm we want results of LER at 354 nm
               call VLIDORT_LER(SCAT, radiance_VL(j,i), refl(j), rc, &
                                icalc = icalc354(j), &
                                transmissivity = trans354(j), &
                                sphericalalbedo = spher354(j))
!write(*,'(1x,i3,4x,3(2x,f8.5))') 354, icalc354(j), trans354(j), spher354(j)
          endif

! Calculations for 388 nm
          if (abs(SCAT%wavelength - 388.00) < 0.01)  then  
                 
              sfc_albedo_388nm = albedo(j,i)
              if ((.not. IS_MISSING(sfc_albedo_354nm)).and. &  ! calcul of S, T and R* at 388 nm
                  (.not. IS_MISSING(sfc_albedo_388nm)).and. &
                  (.not. IS_MISSING(radiance_VL(j,i)) ))  then

                call VLIDORT_LER(SCAT, radiance_VL(j,i), refl(j), rc, &
                                 icalc = icalc388(j), &
                                 transmissivity = trans388(j), &
                                 sphericalalbedo = spher388(j))

!write(*,'(1x,i3,4x,3(2x,f8.5))') 388, icalc388(j), trans388(j), spher388(j)
                
!               PRC: For 1.62.1 algorithm we make a "correction factor" based on LER388
                if(refl(j) .lt. 0.15) then
                 cf = 1.0
                else if(refl(j) .gt. 0.8) then
                 cf = 0.0
                else
                 cf = 1.0-( (refl(j)-0.15) / (0.8-0.15) )
                endif

!               PRC: For 1.62.1 algorithm we make a surface albedo correction for over ocean
                if(qa .eq. 0 .or. qa .ge. 6) then
                  call AI_oceanler(oceanler,solar_zenith(j), sensor_zenith(j), relat_azymuth(j), oceanlerout )
                  reflectivity_388nm_adj = refl(j)-(  (sfc_albedo_388nm+oceanlerout(2)) &
                                                     -(sfc_albedo_354nm+oceanlerout(1))   )*cf
                  call AI_Local(radiance_VL(j,1), icalc354(j), reflectivity_388nm_adj, trans354(j), spher354(j), AI_VL(j))
!                  write(*,'(1x,i3,4x,5(2x,f8.5))') qa, refl(j), sfc_albedo_388nm-sfc_albedo_354nm, &
!                                                                       oceanlerout(2)-oceanlerout(1), &
!                                                   reflectivity_388nm_adj, AI_VL(j)
                else
                 reflectivity_388nm_adj = refl(j)-(sfc_albedo_388nm -sfc_albedo_354nm)*cf
                endif
!   print *, 'pete0', j, refl(j), reflectivity_388nm_adj
              else              
                reflectivity_388nm_adj = MISSING
              end if   
                
          end if  ! end calculation of reflectivity 388nm
                        
     end do ! end loop over nch
      
!     do calcul AI only at 354 nm using refl_388nm (if not = MISSING) -> i = 1  
        
        if (IS_MISSING(radiance_VL(j,1)) .or. IS_MISSING(reflectivity_388nm_adj)) then
            AI_VL(j) = MISSING        
            refl(j)  = MISSING
        else      
            SCAT%wavelength = channels(1)                  
            SCAT%tau => tau(:,1,j)
            SCAT%ssa => ssa(:,1,j)
            SCAT%g => g(:,1,j) 
            SCAT%pmom => pmom(:,:,:,1,j)        
!           This is the so-called Corrected AI
            AI_VL(j) = MISSING
!           PRC: For 1.62.1 algorithm we calculate AI using expression from expression based on icalc354
!            call VLIDORT_AI(SCAT, radiance_VL(j,1),reflectivity_388nm_adj, AI_VL(j), rc)
            call AI_Local(radiance_VL(j,1), icalc354(j), reflectivity_388nm_adj, trans354(j), spher354(j), AI_VL(j))
!           This is the so-called Uncorrected AI - aka, the residue
            residue(j) = MISSING
            if(.not. IS_MISSING(refl(j))) &
!             call VLIDORT_AI(SCAT, radiance_VL(j,1),refl(j), residue(j), rc)    
             call AI_Local(radiance_VL(j,1), icalc354(j), refl(j), trans354(j), spher354(j), residue(j))
!   print *, 'pete1', j, refl(j), reflectivity_388nm_adj
        end if     

!        write(*,'(1x,i6,4x,5(2x,f8.5))') j, radiance_VL(j,1), reflectance_VL(j,1),  AOT(j,1),SSA_VL(j,1), AI_VL(j)

  end do ! end loop over nobs
  
end subroutine AI_Vector

subroutine AI_Local (obs, mod, ref, trans, spher, ai)
 real*8  :: obs     ! "observed" reflectance
 real*8  :: mod     ! "modeled" reflectance for pure Rayleigh atmosphere
 real*8  :: ref     ! adjusted surface reflectivity to use in calculation
 real*8  :: trans   ! atmospheric transmission
 real*8  :: spher   ! atmospheric spherical albedo
 real*8  :: ai      ! output aerosol index
 real*8  :: calrad  ! calculated radiance based on inputs

 calrad = mod + (ref*trans/(1.0-ref*spher))
 ai     = -100. * log10(obs/calrad)
! write(*,'(7(f8.5,2x))')  obs, mod, ref, trans, spher, calrad, ai

end subroutine AI_Local

!=======  Find ocean LER with given geometry =============================
! PRC: This is a correction factor as a function of viewing geometry to
!      be applied to the surface albedos over ocean points for Fresnel
!      reflection (routine provided by S. Gasso)
  subroutine AI_oceanler(oceanler,sza_in, vza_in, raa_in, out_ler )

    IMPLICIT NONE
!
    real*8,                INTENT(IN)  :: sza_in, vza_in, raa_in
    real*8                             :: szain, vzain, raain
    real*8,                INTENT(IN)  :: oceanler(2,16,16,16)
    real*8, DIMENSION(2), INTENT(OUT)  :: out_ler
    real*8                             :: frac_sza, frac_vza, frac_raa
    real*8                             :: slope, bias
    real*8, DIMENSION(2)               :: interp_oceanler
    integer                            :: it, i,j,k, wave_in
    real*8, DIMENSION(16) :: sza_list=(/1.5, 6.0, 12.0,18.0,&
                                        24.0, 30.0, 36.0, 42.0,48.0, 54.0, &
                                        60.0, 66.0, 72.0, 78.0, 84.0, 88.5/)
    real*8, DIMENSION(16) :: vza_list=(/1.5, 6.0, 12.0,18.0,&
                                        24.0, 30.0, 36.0, 42.0,48.0, 54.0, &
                                        60.0, 66.0, 72.0, 78.0, 84.0, 88.5/)
    real*8, DIMENSION(16) :: raa_list=(/0.0, 12.0, 24.0, &
                                        36.0, 48.0, 60.0, 72.0, 84.0, 96.0, 108.0,&
                                        120.0, 132.0, 144.0, 156.0, 168.0, 180.0/)
      !-- Input angles -----------
      szain = sza_in
      vzain = vza_in
      raain = raa_in
      if (sza_in <= 1.50) szain = 1.51
      if (sza_in >= 88.50) szain = 88.49
      if (vza_in <= 1.50) vzain = 1.51
      if (vza_in >= 88.50) vzain = 88.49
      if (raa_in <= 0.0) raain = 0.01
      if (raa_in >= 180.0) raain = 179.9
!
      DO it = 1, 2
        wave_in = it
        !-- Determine interpolation indices ----
        i = 1
        DO
          if (sza_list(i) > szain) EXIT
          i = i +1
        END DO
        frac_sza = (sza_list(i) - szain)/( sza_list(i) - sza_list(i-1) )
!
        j = 1
        DO
          if (vza_list(j) > vzain) EXIT
          j = j +1
        END DO
        frac_vza = (vza_list(j) - vzain)/( vza_list(j) - vza_list(j-1) )
!
        k = 1
        DO
          if (raa_list(k) > raain) EXIT
          k = k +1
        END DO
        frac_raa = (raa_list(k) - raain)/( raa_list(k) - raa_list(k-1) )
!
!write(*,*) ' size(oceanler) ', size(oceanler)
!write(*,*) 'sza_in,vza_in,raa_in', sza_in,vza_in,raa_in
!write(*,*)  'i, j, k ', i, j, k

!
        interp_oceanler(it) =  (1.0 - frac_sza)*(1.0 - frac_vza)* &
                               (1.0 - frac_raa)*oceanler(wave_in, i,j,k) &
        + (frac_sza)*(1.0 - frac_vza)*(1.0 - frac_raa)*oceanler(wave_in, i-1,j,k) &
        + (frac_sza)*(frac_vza)*(1.0 - frac_raa)*oceanler(wave_in, i-1,j-1,k) &
        + (frac_sza)*(frac_vza)*(frac_raa)*oceanler(wave_in, i-1,j-1,k-1) &
        + (frac_sza)*(1.0 - frac_vza)*(frac_raa)*oceanler(wave_in, i-1,j,k-1) &
        + (1.0 - frac_sza)*(frac_vza)*(1.0 - frac_raa)*oceanler(wave_in, i,j-1,k) &
        + (1.0 - frac_sza)*(frac_vza)*(frac_raa)*oceanler(wave_in, i,j-1,k-1) &
        + (1.0 - frac_sza)*(1.0 - frac_vza)*(frac_raa)*oceanler(wave_in, i,j,k-1)

!write(*,*) 'Interpolated ocean LER : ',  interp_oceanler
!write(*,*) 'frac_sza,frac_vza,frac_raa ', frac_sza,frac_vza,frac_raa

      END DO
     !--- Compute LER at 354.0 nm (interpolated) and 388.0 nm (extrapolated) ----
      slope = (interp_oceanler(1) -interp_oceanler(2))/(380.0 - 340.0)
      bias = interp_oceanler(2) - slope*340.0
!
      out_ler(1) = slope*354.0 + bias  ! 354 nm LER
      out_ler(2) = slope*388.0 + bias  ! 388 nm LER
!
      IF (out_ler(1) < 0.0 .OR. out_ler(1) > 0.5) out_ler(1) = -999.0
      IF (out_ler(2) < 0.0 .OR. out_ler(2) > 0.5) out_ler(2) = -999.0
!
!     write(*,*) 'out_ler', out_ler
!
END subroutine AI_oceanler
