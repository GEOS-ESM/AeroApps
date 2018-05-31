module VLIDORT_BRDF_CX
!
!  Simple wrapper for the Python interface to VLIDORT for ocean surfaces
!
!.............................................................................

implicit NONE

PUBLIC VLIDORT_Scalar_CX
PUBLIC VLIDORT_Vector_CX
contains



logical function IS_MISSING(x,MISSING)
  real*8, intent(in)     :: x
  real*8, intent(in)     :: MISSING

  IS_MISSING = abs(x/MISSING-1)<0.001
  return
end function IS_MISSING



subroutine VLIDORT_Scalar_CX (km, nch, nobs,channels, nMom, &
                   nPol, tau, ssa, g, pmom, pe, he, te, U10m, V10m, &
                   mr, solar_zenith, relat_azymuth, sensor_zenith, &
                   MISSING,verbose,radiance_VL_SURF,reflectance_Vl_SURF, ROT, BRDF,rc)
!
! Uses VLIDORT in scalar mode to compute OMI aerosol TOA radiances.
!
  use VLIDORT_ScatMod

  implicit NONE

  logical, parameter            :: scalar = .true.

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

  real*8, target,   intent(in)  :: U10m(nobs)   ! Wind speed components [m/s]
  real*8, target,   intent(in)  :: V10m(nobs)    
  real*8, target,   intent(in)  :: mr(nch)       ! refractive index
  
  real*8, target,   intent(in)  :: MISSING          ! MISSING VALUE                                      
  real*8, target,   intent(in)  :: solar_zenith(nobs)  
  real*8, target,   intent(in)  :: relat_azymuth(nobs) 
  real*8, target,   intent(in)  :: sensor_zenith(nobs) 
  
  integer,          intent(in)  :: verbose

! !OUTPUT PARAMETERS:

  real*8,           intent(out) :: radiance_VL_SURF(nobs,nch)       ! TOA normalized radiance from VLIDORT
  integer,          intent(out) :: rc                          ! return code
  real*8,           intent(out) :: reflectance_VL_SURF(nobs, nch)   ! TOA reflectance from VLIDORT
  real*8,           intent(out) :: ROT(km,nobs,nch)                 ! rayleigh optical thickness  
  real*8,           intent(out) :: BRDF(1,nobs, nch)  
!                         ---  
  integer             :: i,j,n,p,ier

  type(VLIDORT_scat) :: SCAT
  type(VLIDORT_output_scalar)  :: output  

  rc = 0
  ier = 0

  call VLIDORT_Init( SCAT%Surface%Base, km, rc)
  if ( rc /= 0 ) return

  do j = 1,nobs

     ! Make sure albedo and angles are available
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
 
      ! Loop over channels
      ! ------------------
      do i = 1, nch 
       
          ! Mare sure winds and mr are defines
          ! --------------------------------------
          if ( IS_MISSING(U10m(j),MISSING)  .OR. & 
                IS_MISSING(V10m(j),MISSING) .OR. &
                IS_MISSING(mr(i),MISSING)  )  then

              radiance_VL_SURF(j,i) = MISSING
              reflectance_VL_SURF(j,i) = MISSING
              cycle
          end if

          if ( verbose > 0 ) then
            print*, 'DO COX MUNK'
          end if
          call VLIDORT_GissCoxMunk(SCAT%Surface,U10m(j),V10m(j),mr(i),solar_zenith (j),&
                                    sensor_zenith(j),relat_azymuth(j),scalar,rc)
          if ( rc /= 0 ) return

          SCAT%wavelength = channels(i)
          SCAT%tau => tau(:,i,j)
          SCAT%ssa => ssa(:,i,j)
          SCAT%g => g(:,i,j)
          SCAT%pmom => pmom(:,i,j,:,:)

          BRDF(1,j,i) = SCAT%Surface%Base%VIO%VBRDF_Sup_Out%BS_DBOUNCE_BRDFUNC(1,1,1,1)     

          call VLIDORT_Run_Scalar (SCAT, output, ier)

          radiance_VL_SURF(j,i)    = output%radiance
          reflectance_VL_SURF(j,i) = output%reflectance
          ROT(:,j,i) = SCAT%rot      

           if ( ier /= 0 ) then
              radiance_VL_SURF(j,i) = MISSING
              reflectance_VL_SURF(j,i) = MISSING
              cycle
           end if


        end do ! end loop over channels
     
        if ( verbose > 0 ) then
           if ( mod(j-1,1000) == 0 ) then
              print *, '<> VLIDORT Scalar: ', nint(j*100./nobs), '%'
           end if
        end if

  end do ! Loop over obs

end subroutine VLIDORT_Scalar_CX

!..........................................................................

subroutine VLIDORT_Vector_CX (km, nch, nobs,channels, nMom, &
                   nPol, tau, ssa, pmom, pe, he, te, U10m, V10m, &
                   mr, solar_zenith, relat_azymuth, sensor_zenith, &
                   MISSING,verbose,radiance_VL_SURF,reflectance_VL_SURF, ROT, Q, U, BRDF,rc)
!
! Place holder.
!
   use VLIDORT_ScatMod
 
   implicit NONE

  logical                                 :: scalar

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
  real*8, target,   intent(in)  :: pmom(km,nMom,nPol,nch,nobs) !components of the scat phase matrix

  real*8, target,   intent(in)  :: MISSING          ! MISSING VALUE
  real*8, target,   intent(in)  :: pe(km+1,nobs)    ! pressure at layer edges [Pa]
  real*8, target,   intent(in)  :: he(km+1,nobs)    ! height above sea-level  [m]
  real*8, target,   intent(in)  :: te(km+1,nobs)    ! temperature at layer edges [K]

  real*8, target,   intent(in)  :: U10m(nobs)   ! Wind speed components [m/s]
  real*8, target,   intent(in)  :: V10m(nobs)    
  real*8, target,   intent(in)  :: mr(nch)       ! refractive index
                       
  real*8, target,   intent(in)  :: solar_zenith(nobs)  
  real*8, target,   intent(in)  :: relat_azymuth(nobs) 
  real*8, target,   intent(in)  :: sensor_zenith(nobs) 

  integer,          intent(in)  :: verbose

! !OUTPUT PARAMETERS:

  real*8,           intent(out) :: radiance_VL_SURF(nobs,nch)       ! TOA normalized radiance from VLIDORT
  integer,          intent(out) :: rc                          ! return code
  real*8,           intent(out) :: reflectance_VL_SURF(nobs, nch)   ! TOA reflectance from VLIDORT
  real*8,           intent(out) :: ROT(km,nobs,nch)               ! rayleigh optical thickness  
  real*8,           intent(out) :: BRDF(3,nobs, nch)  
  real*8,           intent(out) :: Q(nobs, nch)   ! Stokes parameter Q
  real*8,           intent(out) :: U(nobs, nch)   ! Stokes parameter U
!                               ---
  
  integer             :: i,j,n,p,ier 
  
  type(VLIDORT_scat) :: SCAT
  type(VLIDORT_output_vector)  :: output  

  rc = 0
  ier = 0
 
  call VLIDORT_Init( SCAT%Surface%Base, km, rc)
  if ( rc /= 0 ) return

  SCAT%nMom = nMom
  SCAT%nPol = nPol
  SCAT%NSTOKES = 3
  if ( SCAT%NSTOKES  .GT. MAXSTOKES  )   return

  do j = 1, nobs
     
     ! Make sure albedo and angles are available
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
          ! Mare sure winds and mr are defines
          ! --------------------------------------
          if ( IS_MISSING(U10m(j),MISSING)  .OR. & 
                IS_MISSING(V10m(j),MISSING) .OR. &
                IS_MISSING(mr(i),MISSING)  )  then

              radiance_VL_SURF(j,i) = MISSING
              reflectance_VL_SURF(j,i) = MISSING
              cycle
          end if

          if ( verbose > 0 ) then
            print*, 'DO COX MUNK'
          end if
          scalar = .false.
          call VLIDORT_GissCoxMunk(SCAT%Surface,U10m(j),V10m(j),mr(i),solar_zenith (j),&
                                    sensor_zenith(j),relat_azymuth(j),scalar,rc)
           if ( rc /= 0 ) return

           SCAT%wavelength = channels(i)
           SCAT%tau => tau(:,i,j)
           SCAT%ssa => ssa(:,i,j)
           SCAT%pmom => pmom(:,i,j,:,:)

          BRDF(1,j,i) = SCAT%Surface%Base%VIO%VBRDF_Sup_Out%BS_DBOUNCE_BRDFUNC(1,1,1,1) 
          BRDF(2,j,i) = SCAT%Surface%Base%VIO%VBRDF_Sup_Out%BS_DBOUNCE_BRDFUNC(2,1,1,1) 
          BRDF(3,j,i) = SCAT%Surface%Base%VIO%VBRDF_Sup_Out%BS_DBOUNCE_BRDFUNC(3,1,1,1)     

         
          call VLIDORT_Run_Vector (SCAT, output, ier)


          ROT(:,j,i) = SCAT%rot
          radiance_VL_SURF(j,i)    = output%radiance
          reflectance_VL_SURF(j,i) = output%reflectance
          Q(j,i)                   = output%Q
          U(j,i)                   = output%U                

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

end subroutine VLIDORT_Vector_Cx

end module VLIDORT_BRDF_CX