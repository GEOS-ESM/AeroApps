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
  use VLIDORT_LambMod

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
 
  type(VLIDORT_Lamb) :: lamb

#define IS_MISSING(x) (abs(x/MISSING-1)<0.001)
  rc = 0
  ier = 0

  call VLIDORT_LambInit( lamb, km, rc, .true. )
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

      lamb%solar_zenith = solar_zenith (j)
      lamb%relative_azymuth = relat_azymuth(j)
      lamb%sensor_zenith = sensor_zenith(j) 
      lamb%pe => pe(:,j)
      lamb%ze => he(:,j)
      lamb%te => te(:,j) 
 
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

           lamb%wavelength = channels(i)
           lamb%albedo = albedo(j,i)
           lamb%tau => tau(:,i,j)
           lamb%ssa => ssa(:,i,j)
           lamb%g => g(:,i,j)
         
           call VLIDORT_LambRun (lamb, radiance_VL(j,i), reflectance_VL(j,i), &
                                 AOT(j,i), ier, .true., .true.)
           
           ! add if SZA > 90 deg --> rad =0 for lidar calc

           if ( ier /= 0 ) then
              if (lamb%solar_zenith > 89.9) then
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
                   MISSING,verbose, radiance_VL, reflectance_VL, rc)
!
! Place holder.
!
   use VLIDORT_LambMod
 
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

!                               ---
  
  integer             :: i,j,ier
  real*8              :: AOT(nobs,nch)            ! total aerosol optical thickness 

  
  type(VLIDORT_Lamb) :: lamb

#define IS_MISSING(x) (abs(x/MISSING-1)<0.001)
  rc = 0
  ier = 0

  call VLIDORT_LambInit( lamb, km, rc, .false. )
  if ( rc /= 0 ) return

  lamb%nmom = nMom
 
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

      lamb%solar_zenith = solar_zenith (j)
      lamb%relative_azymuth = relat_azymuth(j)
      lamb%sensor_zenith = sensor_zenith(j) 
      lamb%pe => pe(:,j)
      lamb%ze => he(:,j)
      lamb%te => te(:,j) 
 
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

           lamb%wavelength = channels(i)
           lamb%albedo = albedo(j,i)
           lamb%tau => tau(:,i,j)
           lamb%ssa => ssa(:,i,j)
           lamb%g => g(:,i,j)
           lamb%pmom => pmom(:,:,:,i,j)

        call VLIDORT_LambRun (lamb, radiance_VL(j,i),reflectance_VL(j,i),&
                              AOT(j,i), ier, .false., .true.)
        ! add if SZA > 90 deg --> rad =0 for lidar calc
        if ( ier /= 0 ) then
              if (lamb%solar_zenith > 89.9) then
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
!(Check in VLIDORT f77 if the subroutines have true for scalar option)
! It is better to use AI_vector for accuracy

  use VLIDORT_LambMod

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
  real*8,           intent(out) :: refl(nobs, nch)             ! surface reflectivity R*
  
!                               ---
  
  integer             :: n, i,j, k,ier
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

#define IS_MISSING(x) (abs(x/MISSING-1)<0.001)

  sfc_albedo_354nm = MISSING
  rc = 0
  ier = 0

  call VLIDORT_LambInit( lamb, km, rc, .true. )
  if ( rc /= 0 ) return
 
  do j = 1, nobs

     ! Make sure albedo and angles are available
     ! -----------------------------------------
     if ( IS_MISSING(solar_zenith(j))  .OR. & 
          IS_MISSING(sensor_zenith(j)) .OR. &
          IS_MISSING(relat_azymuth(j))  )  then

         radiance_VL(j,i) = MISSING
         AI_VL(j) = MISSING
         refl(j,i) = MISSING
        cycle

     end if

      lamb%solar_zenith = solar_zenith (j)
      lamb%relative_azymuth = relat_azymuth(j)
      lamb%sensor_zenith = sensor_zenith(j) 
      lamb%pe => pe(:,j)
      lamb%ze => he(:,j)
      lamb%te => te(:,j) 
 
      ! Loop over channels
      ! ------------------
      do i = 1, nch 
          
       ! Mare sure albedo is defined
       ! ---------------------------
        if ( IS_MISSING(albedo(j,i)) ) then
        radiance_VL(j,i) = MISSING
        AI_VL(j) = MISSING
        refl(j,i) = MISSING
        else
        
        lamb%wavelength = channels(i)
        lamb%albedo = albedo(j,i)
        lamb%tau => tau(:,i,j)
        lamb%ssa => ssa(:,i,j)
        lamb%g => g(:,i,j)
        
        call VLIDORT_LambRun (lamb, radiance_VL(j,i),reflectance_VL(j,i),&
                              AOT(j,i), ier, .true., .true.)
             if ( ier /= 0 ) then
                 radiance_VL(j,i) = MISSING
                 AI_VL(j) = MISSING
                 refl(j,i) = MISSING           
                 cycle
             end if        

        if (abs(lamb%wavelength - 354.00) < 0.01)  sfc_albedo_354nm = albedo(j,i) ! surface albedo  at 354 nm
        
        if (abs(lamb%wavelength - 388.00) < 0.01) then           
            if (.not. IS_MISSING(sfc_albedo_354nm))  then  ! calcul of S, T and R* at 388 nm                
                call VLIDORT_LER(lamb, radiance_VL(j,i), refl(j,i), rc)
                sfc_albedo_388nm = albedo(j,i)
                reflectivity_388nm_adj = refl(j,i)-(sfc_albedo_388nm -sfc_albedo_354nm) 
            else              
                reflectivity_388nm_adj = MISSING
            end if       
        end if  ! end calcul of reflectivity 388nm
        
        end if  ! end if albedo is MISSING            
     end do ! end loop over nch
      
!     do calcul AI only at 354 nm using refl_388nm (if not = MISSING) -> i = 1  
        
        if (IS_MISSING(reflectivity_388nm_adj)) then
            AI_VL(j) = MISSING   
 
        else      
        lamb%wavelength = channels(1)          
        lamb%albedo = albedo(j,1)
        lamb%tau => tau(:,1,j)
        lamb%ssa => ssa(:,1,j)
        lamb%g => g(:,1,j) 
        
        call VLIDORT_AI(lamb, radiance_VL(j,1),reflectivity_388nm_adj, AI_VL(j), rc)
       
        end if     
!     end do 
  
    
  end do

end subroutine AI_Scalar
!..........................................................................

subroutine AI_Vector (km, nch, nobs,  channels, nMom,  &
                   nPol,tau, ssa, g, pmom, pe, he, te,albedo,     &
                   solar_zenith, relat_azymuth, sensor_zenith, &
                   MISSING,radiance,AI, verbose, radiance_VL,  &
                   AI_VL, refl, rc)
!
! Place holder.
!
   use VLIDORT90_LambMod
 
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
  real*8, target,   intent(in)  :: MISSING 


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
  real*8              :: AOT(nobs,nch)            ! total aerosol optical thickness 
  real*8              :: g_tot(nobs,nch)            ! total asymetry factor 
  real*8              :: sfc_albedo_354nm         ! surface albedo 354 nm OMI file (climato TOMS)
  real*8              :: sfc_albedo_388nm         ! surface albedo 388 nm OMI file (climato TOMS)
  real*8              :: raymodel(nobs)           ! radiance Ray 354 nm from VLIDORT
  real*8              :: ray_OMI                  ! radiance Ray 354 nm OMI file
  real*8              :: reflectance_VL(nobs,nch)    ! TOA reflectance from VLIDORT

#define IS_MISSING(x) (abs(x/MISSING-1)<0.001)  

  type(VLIDORT_Lamb) :: lamb

  sfc_albedo_354nm = MISSING
  rc = 0
  ier = 0
 
  call VLIDORT_LambInit( lamb, km, rc, .false. )
  if ( rc /= 0 ) return

  lamb%nmom = nMom
 
  do j = 1, nobs

     ! Make sure albedo and angles are available
     ! -----------------------------------------
     if ( IS_MISSING(solar_zenith(j))  .OR. & 
          IS_MISSING(sensor_zenith(j)) .OR. &
          IS_MISSING(relat_azymuth(j))  )  then

         radiance_VL(j,i) = MISSING
         AI_VL(j) = MISSING
         refl(j,i) = MISSING
        cycle

      end if

      lamb%solar_zenith = solar_zenith (j)
      lamb%relative_azymuth = relat_azymuth(j)
      lamb%sensor_zenith = sensor_zenith(j) 
      lamb%pe => pe(:,j)
      lamb%ze => he(:,j)
      lamb%te => te(:,j)

     do i = 1, nch
       if (abs((albedo(j,i)/MISSING)-1.) < 0.001) then
          radiance_VL(j,i) = MISSING
          AI_VL(j) = MISSING
          refl(j,i) = MISSING
        else
          lamb%wavelength = channels(i)
          lamb%albedo = albedo(j,i)
          lamb%tau => tau(:,i,j)
          lamb%ssa => ssa(:,i,j)
          lamb%g => g(:,i,j)
          lamb%pmom => pmom(:,:,:,i,j)
        
          call VLIDORT_LambRun (lamb, radiance_VL(j,i),reflectance_VL(j,i),&
                                AOT(j,i), ier, .false., .true.)
       
               if ( ier /= 0 ) then
                 radiance_VL(j,i) = MISSING
                 AI_VL(j) = MISSING
                 refl(j,i) = MISSING           
                 cycle
             end if  


          if (abs(lamb%wavelength - 354.00) < 0.01)  sfc_albedo_354nm = albedo(j,i)  ! surface albedo at 354 nm
       
          if (abs(lamb%wavelength - 388.00) < 0.01) then           
              if (.not. IS_MISSING(sfc_albedo_354nm)) then  ! calcul of S, T and R* at 388 nm
                
                call VLIDORT_LER(lamb, radiance_VL(j,i), refl(j,i), rc)
                sfc_albedo_388nm = albedo(j,i)
                reflectivity_388nm_adj = refl(j,i)-(sfc_albedo_388nm -sfc_albedo_354nm) 
              else              
                reflectivity_388nm_adj = MISSING
              end if       
          end if  ! end calcul of reflectivity 388nm
        
       end if  ! end if albedo is MISSING            
     end do ! end loop over nch
      
!     do calcul AI only at 354 nm using refl_388nm (if not = MISSING) -> i = 1  
        
        if (IS_MISSING(reflectivity_388nm_adj)) then
            AI_VL(j) = MISSING   
        else      
        lamb%wavelength = channels(1)          
        lamb%albedo = albedo(j,1)
        lamb%tau => tau(:,1,j)
        lamb%ssa => ssa(:,1,j)
        lamb%g => g(:,1,j) 
        lamb%pmom => pmom(:,:,:,1,j)
        
        call VLIDORT_AI(lamb, radiance_VL(j,1),reflectivity_388nm_adj, AI_VL(j), rc)
       
        end if     
!     end do 
  
  end do
  
end subroutine AI_Vector
