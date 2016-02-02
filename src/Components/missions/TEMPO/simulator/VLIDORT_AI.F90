!.............................................................................

subroutine VLIDORT_AI_Scalar (km, nch, nobs, channels,        &
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

  logical, parameter            :: aerosol = .true.
  logical, parameter            :: scalar = .true.
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
 
  real*8              :: ROT(km,nobs,nch)            ! rayleigh optical thickness 
  real*8              :: g_tot(nobs,nch)          ! total asymetry factor (column - one profile) 
  real*8              :: sfc_albedo_354nm         ! surface albedo 354 nm OMI file (climato TOMS)
  real*8              :: sfc_albedo_388nm         ! surface albedo 388 nm OMI file (climato TOMS)
  real*8              :: raymodel(nobs)           ! radiance Ray 354 nm from VLIDORT
  real*8              :: ray_OMI                  ! radiance Ray 354 nm OMI file
  real*8              :: reflectance_VL(nobs,nch)    ! TOA reflectance from VLIDORT
  real*8              :: Q(nobs,nch), U(nobs,nch)
  
  type(VLIDORT_scat) :: SCAT

! #define IS_MISSING(x) (abs(x/MISSING-1)<0.001)

 ! sfc_albedo_354nm = MISSING
 ! ier = 0
  rc = 0
  
!  call VLIDORT_Init(  SCAT%Surface%Base, km, rc)
!  if ( rc /= 0 ) return
  
  radiance_VL = 0
  AI_VL = 0
  refl = 0
  ! do j = 1, nobs
  !     if ( IS_MISSING(solar_zenith(j))  .OR. & 
  !         IS_MISSING(sensor_zenith(j)) .OR. &
  !         IS_MISSING(relat_azymuth(j))  )  then

  !       radiance_VL(j,:) = MISSING
  !       refl(j) = MISSING
  !       AI_VL(j) = MISSING
  !       cycle

  !     end if
         
  !    SCAT%pe => pe(:,j)
  !    SCAT%ze => he(:,j)
  !    SCAT%te => te(:,j) 
 
!      do i = 1, nch 
       
!         ! Make sure albedo is defined
!         ! ---------------------------
!            if ( IS_MISSING(albedo(j,i)) ) then
!               radiance_VL(j,i) = MISSING
!               refl(j) = MISSING
!               AI_VL(j) = MISSING
!               cycle
!            end if
           
!         call VLIDORT_SurfaceLamb(SCAT%Surface,albedo(j,i),solar_zenith (j),sensor_zenith(j),&
!                                relat_azymuth(j),scalar)

!         SCAT%wavelength = channels(i)        
!         SCAT%tau => tau(:,i,j)
!         SCAT%ssa => ssa(:,i,j)
!         SCAT%g => g(:,i,j)
                 
!         call VLIDORT_Run (SCAT, radiance_VL(j,i),reflectance_VL(j,i), &
!                           ROT(:,j,i), Q(j,i),U(j,i),scalar, aerosol, ier)

!         if ( ier /= 0 ) then
!                  radiance_VL(j,i) = MISSING
!                  AI_VL(j) = MISSING
!                  refl(j) = MISSING           
!                  cycle
!         end if 
               
!         if (abs(SCAT%wavelength - 354.00) < 0.01)  sfc_albedo_354nm = albedo(j,i) ! surface albedo  at 354 nm
        
!         if (abs(SCAT%wavelength - 388.00) < 0.01) then           
!             if ((.not. IS_MISSING(sfc_albedo_354nm)) .and.&
!                (.not. IS_MISSING(radiance_VL(j,i)) ))  then ! calcul of S, T and R* at 388 nm
                
!                 call VLIDORT_LER(SCAT, radiance_VL(j,i), refl(j), rc)
!                 sfc_albedo_388nm = albedo(j,i)
!                 reflectivity_388nm_adj = refl(j)-(sfc_albedo_388nm -sfc_albedo_354nm) 
!             else              
!                 reflectivity_388nm_adj = MISSING
!             end if       
!         end if  ! end calcul of reflectivity 388nm
        
         
!      end do ! end loop over nch
      
! !     do calcul AI only at 354 nm using refl_388nm (if not = MISSING) -> i = 1  
        
!         if (IS_MISSING(reflectivity_388nm_adj)) then
!             AI_VL(j) = MISSING        
!         else      
!         SCAT%wavelength = channels(1)                  
!         SCAT%tau => tau(:,1,j)
!         SCAT%ssa => ssa(:,1,j)
!         SCAT%g => g(:,1,j) 
 
!         call VLIDORT_AI(SCAT, radiance_VL(j,1),reflectivity_388nm_adj, AI_VL(j), rc)    
        
!         end if 
!   end do

end subroutine VLIDORT_AI_Scalar
!..........................................................................

! subroutine VLIDORT_AI_Vector (km, nch, nobs,  channels, nMom,  &
!                    nPol,tau, ssa, g, pmom, pe, he, te,albedo,     &
!                    solar_zenith,relat_azymuth, sensor_zenith,    &
!                    MISSING,verbose, radiance_VL,  &
!                    AI_VL, refl, rc)
! !
! ! Place holder.
! !
!    use VLIDORT_ScatMod
 
!    implicit NONE

!   logical, parameter            :: aerosol = .true.
!   logical, parameter            :: scalar = .false.
! ! !INPUT PARAMETERS:

!   integer,          intent(in)  :: km    ! number of levels on file
!   integer,          intent(in)  :: nch   ! number of channels
!   integer,          intent(in)  :: nobs  ! number of observations

!   integer, target,  intent(in)  :: nMom  ! number of moments 
!   integer, target,  intent(in)  :: nPol  ! number of components                               
                  
!   real*8, target,   intent(in)  :: channels(nch)    ! wavelengths [nm]

! !                                                   ! --- Mie Parameters ---
!   real*8, target,   intent(in)  :: tau(km,nch,nobs) ! aerosol optical depth
!   real*8, target,   intent(in)  :: ssa(km,nch,nobs) ! single scattering albedo
!   real*8, target,   intent(in)  :: g(km,nch,nobs)   ! asymmetry factor
  
!   real*8, target,   intent(in)  :: pmom(km,nMom,nPol,nch,nobs) !components of the scat phase matrix


!   real*8, target,   intent(in)  :: pe(km+1,nobs)    ! pressure at layer edges [Pa]
!   real*8, target,   intent(in)  :: he(km+1,nobs)    ! height above sea-level  [m]
!   real*8, target,   intent(in)  :: te(km+1,nobs)    ! temperature at layer edges [K]
                       
!   real*8, target,   intent(in)  :: solar_zenith(nobs)  
!   real*8, target,   intent(in)  :: relat_azymuth(nobs) 
!   real*8, target,   intent(in)  :: sensor_zenith(nobs) 
!   real*8, target,   intent(in)  :: MISSING 


!   real*8, target,   intent(in)  :: albedo(nobs,nch)       ! surface albedo
! !  real*8,           intent(in)  :: radiance(nobs,nch)     ! radiance observed OMI
! !  real*8,           intent(in)  :: reflectivity(nobs,nch)! reflectivity fichier OMI
! !  real*8,           intent(in)  :: AI(nobs)               ! Aerosol Index OMI file
  
!   integer,          intent(in)  :: verbose

! ! !OUTPUT PARAMETERS:

!   real*8,           intent(out) :: radiance_VL(nobs,nch)       ! TOA radiance from VLIDORT
!   integer,          intent(out) :: rc                          ! return code
!   real*8,           intent(out) :: AI_VL(nobs)                 ! Aerosol Index model VLIDORT
!   real*8,           intent(out) :: refl(nobs)                  ! surface reflectivity R* at 388nm
 
!   real*8                        :: ROT(km,nobs,nch)            ! rayleigh optical thickness  
  
!   integer             :: n, i,j, k, imom, ipol, i_354,ier
!   real*8              :: spher_alb(nobs,nch)      ! spherical albedo
!   real*8              :: trans(nobs, nch)         ! transmission
!   real*8              :: reflectivity_388nm       ! Reflectivity Vlidort 388 nm
!   real*8              :: reflectivity_354nm       ! Reflectivity Vlidort 354 nm
!   real*8              :: reflectivity_388nm_OMI   ! Reflectivity OMI file
!   real*8              :: reflectivity_354nm_OMI   ! Reflectivity OMI file
!   real*8              :: reflectivity_388nm_adj   ! Reflectivity OMI file
!   real*8              :: delta_R_model            ! Difference R
!   real*8              :: delta_R_OMI              ! Difference R
 
!   real*8              :: g_tot(nobs,nch)            ! total asymetry factor 
!   real*8              :: sfc_albedo_354nm         ! surface albedo 354 nm OMI file (climato TOMS)
!   real*8              :: sfc_albedo_388nm         ! surface albedo 388 nm OMI file (climato TOMS)
!   real*8              :: raymodel(nobs)           ! radiance Ray 354 nm from VLIDORT
!   real*8              :: ray_OMI                  ! radiance Ray 354 nm OMI file
!   real*8              :: reflectance_VL(nobs,nch)    ! TOA reflectance from VLIDORT
!   real*8              :: Q(nobs,nch), U(nobs,nch)
  
!   type(VLIDORT_scat) :: SCAT

! #define IS_MISSING(x) (abs(x/MISSING-1)<0.001)  
  
!   sfc_albedo_354nm = MISSING
!   reflectivity_388nm_adj = MISSING
!   rc = 0
!   ier = 0
 
!   call VLIDORT_Init(  SCAT%Surface%Base, km, rc)
!   if ( rc /= 0 ) return

!   SCAT%nMom = nMom
!   SCAT%nPol = nPol

!   do j = 1, nobs
!      if ( IS_MISSING(solar_zenith(j))  .OR. & 
!           IS_MISSING(sensor_zenith(j)) .OR. &
!           IS_MISSING(relat_azymuth(j))  )  then

!         radiance_VL(j,:) = MISSING
!         refl(j) = MISSING
!         AI_VL(j) = MISSING
!         cycle

!      end if
!      print*, 'observation', j     
!      SCAT%pe => pe(:,j)
!      SCAT%ze => he(:,j)
!      SCAT%te => te(:,j) 

!      do i = 1, nch
!        print*, 'nchannel number', i
!        ! Mare sure albedo is defined
!            ! ---------------------------
!            if ( IS_MISSING(albedo(j,i)) ) then
!               radiance_VL(j,i) = MISSING
!               refl(j) = MISSING
!               AI_VL(j) = MISSING
!               cycle
!            end if
          
!           call VLIDORT_SurfaceLamb(SCAT%Surface,albedo(j,i),solar_zenith (j),sensor_zenith(j),&
!                                relat_azymuth(j),scalar)
          
!           SCAT%wavelength = channels(i)          
!           SCAT%tau => tau(:,i,j)
!           SCAT%ssa => ssa(:,i,j)
!           SCAT%g => g(:,i,j)
!           SCAT%pmom => pmom(:,:,:,i,j)
!           print*, 'everything ok' 
!           call VLIDORT_Run (SCAT, radiance_VL(j,i),reflectance_VL(j,i), &
!                             ROT(:,j,i),Q(j,i),U(j,i), scalar, aerosol, ier)
!                if ( ier /= 0 ) then
!                     radiance_VL(j,i) = MISSING
!                     AI_VL(j) = MISSING
!                     refl(j) = MISSING      
!                     cycle
!                end if 

!           if (abs(SCAT%wavelength - 354.00) < 0.01) &
!                sfc_albedo_354nm = albedo(j,i)  ! surface albedo at 354 nm
              
!           if (abs(SCAT%wavelength - 388.00) < 0.01)  then  
                 
!               if ((.not. IS_MISSING(sfc_albedo_354nm)).and. &  ! calcul of S, T and R* at 388 nm
!                  (.not. IS_MISSING(radiance_VL(j,i)) ))  then

!                 sfc_albedo_388nm = albedo(j,i)
!                 call VLIDORT_LER(SCAT, radiance_VL(j,i), refl(j), rc)
                
!                 reflectivity_388nm_adj = refl(j)-(sfc_albedo_388nm -sfc_albedo_354nm) 
!               else              
!                 reflectivity_388nm_adj = MISSING
!               end if   
                
!           end if  ! end calcul of reflectivity 388nm
                        
!      end do ! end loop over nch
!       print*, 'refl', reflectivity_388nm_adj
! !     do calcul AI only at 354 nm using refl_388nm (if not = MISSING) -> i = 1  
        
!         if (IS_MISSING(reflectivity_388nm_adj)) then
!             AI_VL(j) = MISSING        
!         else      
!         SCAT%wavelength = channels(1)                  
!         SCAT%tau => tau(:,1,j)
!         SCAT%ssa => ssa(:,1,j)
!         SCAT%g => g(:,1,j) 
!         SCAT%pmom => pmom(:,:,:,1,j)        
!         call VLIDORT_AI(SCAT, radiance_VL(j,1),reflectivity_388nm_adj, AI_VL(j), rc)    
!         end if     

!   end do ! end loop over nobs
  
! end subroutine VLIDORT_AI_Vector
