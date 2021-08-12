module VLIDORT_BRDF_CX
!
!  Simple wrapper for the Python interface to VLIDORT for ocean surfaces
!
!.............................................................................

implicit NONE

PUBLIC VLIDORT_Scalar_GissCX
PUBLIC VLIDORT_Vector_GissCX
PUBLIC VLIDORT_Vector_GissCx_Cloud
PUBLIC VLIDORT_Scalar_GissCx_NOBM_Cloud
PUBLIC VLIDORT_Vector_GissCx_NOBM_Cloud

PUBLIC VLIDORT_Scalar_CX
PUBLIC VLIDORT_Vector_CX
PUBLIC VLIDORT_Vector_CX_NOBM

PUBLIC VLIDORT_Vector_OCICX
PUBLIC VLIDORT_Vector_OCIGissCX
PUBLIC VLIDORT_Scalar_OCIGissCx_Cloud
PUBLIC VLIDORT_Vector_OCIGissCx_Cloud
PUBLIC VLIDORT_Scalar_OCIGissCx_NOBM_Cloud
PUBLIC VLIDORT_Vector_OCIGissCx_NOBM_Cloud


interface VLIDORT_Vector_CX
  module procedure VLIDORT_Vector_CX_SingleGeom
  module procedure VLIDORT_Vector_CX_MultiGeom
end interface 


interface VLIDORT_Vector_CX_Cloud
  module procedure VLIDORT_Vector_CX_Cloud_SingleGeom
  module procedure VLIDORT_Vector_CX_Cloud_MultiGeom
end interface


contains



logical function IS_MISSING(x,MISSING)
  real*8, intent(in)     :: x
  real*8, intent(in)     :: MISSING

  IS_MISSING = abs(x/MISSING-1)<0.001
  return
end function IS_MISSING


!..........................................................................

subroutine VLIDORT_Vector_OCIGissCX (km, nch, nobs,channels, nstreams, plane_parallel, nMom, &
                   nPol, ROT, depol, alpha, tau, ssa, pmom, pe, he, te, U10m, V10m, &
                   mr, solar_zenith, relat_azymuth, sensor_zenith, flux_factor, &
                   MISSING,verbose,radiance_VL_SURF,reflectance_VL_SURF, &
                   Q, U, BRDF, BRDF_Q, BRDF_U, rc)
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

  logical,          intent(in)  :: plane_parallel

  integer, target,  intent(in)  :: nMom  ! number of moments 
  integer, target,  intent(in)  :: nPol  ! number of components                               
  integer,          intent(in)  :: nstreams  ! number of half space streams
                  
  real*8, target,   intent(in)  :: channels(nch)    ! wavelengths [nm]

!                                                   ! --- Rayleigh Parameters ---
  real*8, target,   intent(in)  :: ROT(km,nobs,nch) ! rayleigh optical thickness
  real*8, target,   intent(in)  :: depol(nch)       ! rayleigh depolarization ratio used in phase matrix

!                                                   ! --- Trace Gas Absorption---
  real*8, target,   intent(in)  :: alpha(km,nobs,nch) ! trace gas absoprtion optical thickness

!                                                   ! --- Mie Parameters ---
  real*8, target,   intent(in)  :: tau(km,nch,nobs) ! aerosol optical depth
  real*8, target,   intent(in)  :: ssa(km,nch,nobs) ! single scattering albedo
  real*8, target,   intent(in)  :: pmom(km,nch,nobs,nMom,nPol) !components of the scat phase matrix

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

  real*8,           intent(in)  :: flux_factor(nch,nobs) ! solar flux (F0)
  integer,          intent(in)  :: verbose

! !OUTPUT PARAMETERS:

  real*8,           intent(out) :: radiance_VL_SURF(nobs,nch)       ! TOA normalized radiance from VLIDORT
  integer,          intent(out) :: rc                          ! return code
  real*8,           intent(out) :: reflectance_VL_SURF(nobs, nch)   ! TOA reflectance from VLIDORT
  real*8,           intent(out) :: BRDF(nobs, nch)
  real*8,           intent(out) :: BRDF_Q(nobs, nch)
  real*8,           intent(out) :: BRDF_U(nobs, nch)  
  real*8,           intent(out) :: Q(nobs, nch)   ! Stokes parameter Q
  real*8,           intent(out) :: U(nobs, nch)   ! Stokes parameter U
!                               ---
  
  integer             :: i,j,n,p,ier 
  
  type(VLIDORT_scat) :: SCAT
  type(VLIDORT_output_vector)  :: output  

  rc = 0
  ier = 0
 
  ! set streams here
  SCAT%Surface%BASE%NSTREAMS = nstreams
  SCAT%Surface%BASE%DO_PLANE_PARALLEL = plane_parallel
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
          ! set solar flux
          SCAT%Surface%Base%VIO%VLIDORT_FixIn%SunRays%TS_FLUX_FACTOR = flux_factor(i,j)

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
            print*,U10m(j),V10m(j),mr(i),solar_zenith (j),&
                                    sensor_zenith(j),relat_azymuth(j)
          end if
          scalar = .false.
          call VLIDORT_OCIGissCoxMunk(SCAT%Surface,U10m(j),V10m(j),mr(i),solar_zenith (j),&
                                    sensor_zenith(j),relat_azymuth(j),scalar,rc)

          if ( verbose > 0 ) then
            print*, 'FINISHED COX MUNK'
          end if


          if ( rc /= 0 ) return

          SCAT%wavelength = channels(i)
          SCAT%rot => ROT(:,j,i)
          SCAT%depol_ratio => depol(i)
          SCAT%alpha => alpha(:,j,i)
          SCAT%tau => tau(:,i,j)
          SCAT%ssa => ssa(:,i,j)
          SCAT%pmom => pmom(:,i,j,:,:)

          if ( verbose > 0 ) then
            print *,'brdf output1', SCAT%Surface%Base%VIO%VBRDF_Sup_Out%BS_DBOUNCE_BRDFUNC(1,1,1,1) 
            print *,'brdf output2', SCAT%Surface%Base%VIO%VBRDF_Sup_Out%BS_DBOUNCE_BRDFUNC(2,1,1,1) 
            print *,'brdf output3', SCAT%Surface%Base%VIO%VBRDF_Sup_Out%BS_DBOUNCE_BRDFUNC(3,1,1,1) 
          end if

          BRDF(j,i) = SCAT%Surface%Base%VIO%VBRDF_Sup_Out%BS_DBOUNCE_BRDFUNC(1,1,1,1) 
          BRDF_Q(j,i) = SCAT%Surface%Base%VIO%VBRDF_Sup_Out%BS_DBOUNCE_BRDFUNC(2,1,1,1) 
          BRDF_U(j,i) = SCAT%Surface%Base%VIO%VBRDF_Sup_Out%BS_DBOUNCE_BRDFUNC(3,1,1,1)     

         
          call VLIDORT_Run_Vector (SCAT, output, ier)


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

end subroutine VLIDORT_Vector_OCIGissCx


!..........................................................................

subroutine VLIDORT_Vector_OCICX (km, nch, nobs,channels, nstreams, plane_parallel, nMom, &
                   nPol, ROT, depol, tau, ssa, pmom, pe, he, te, U10m, V10m, &
                   mr, solar_zenith, relat_azymuth, sensor_zenith, &
                   MISSING,verbose,radiance_VL_SURF,reflectance_VL_SURF, &
                   Q, U, BRDF, BRDF_Q, BRDF_U, rc)
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

  logical,          intent(in)  :: plane_parallel ! do plane parallel flag

  integer, target,  intent(in)  :: nMom  ! number of moments 
  integer, target,  intent(in)  :: nPol  ! number of components    
  integer,          intent(in)  :: nstreams  ! number of half space streams

                  
  real*8, target,   intent(in)  :: channels(nch)    ! wavelengths [nm]

!                                                   ! --- Rayleigh Parameters ---
  real*8, target,   intent(in)  :: ROT(km,nobs,nch) ! rayleigh optical thickness
  real*8, target,   intent(in)  :: depol(nch)       ! rayleigh depolarization ratio used in phase matrix

!                                                   ! --- Mie Parameters ---
  real*8, target,   intent(in)  :: tau(km,nch,nobs) ! aerosol optical depth
  real*8, target,   intent(in)  :: ssa(km,nch,nobs) ! single scattering albedo
  real*8, target,   intent(in)  :: pmom(km,nch,nobs,nMom,nPol) !components of the scat phase matrix

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
  real*8,           intent(out) :: BRDF(nobs, nch)
  real*8,           intent(out) :: BRDF_Q(nobs, nch)
  real*8,           intent(out) :: BRDF_U(nobs, nch)  
  real*8,           intent(out) :: Q(nobs, nch)   ! Stokes parameter Q
  real*8,           intent(out) :: U(nobs, nch)   ! Stokes parameter U
!                               ---
  
  integer             :: i,j,n,p,ier 
  
  type(VLIDORT_scat) :: SCAT
  type(VLIDORT_output_vector)  :: output  

  rc = 0
  ier = 0

  ! set streams here
  SCAT%Surface%Base%NSTREAMS = nstreams 
  SCAT%Surface%Base%DO_PLANE_PARALLEL = plane_parallel
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
            print*,U10m(j),V10m(j),mr(i),solar_zenith (j),&
                                    sensor_zenith(j),relat_azymuth(j)
          end if
          scalar = .false.
          call VLIDORT_OCICoxMunk(SCAT%Surface,U10m(j),V10m(j),mr(i),solar_zenith (j),&
                                    sensor_zenith(j),relat_azymuth(j),scalar,rc)

          if ( verbose > 0 ) then
            print*, 'FINISHED COX MUNK'
          end if


          if ( rc /= 0 ) return

          SCAT%wavelength = channels(i)
          SCAT%rot => ROT(:,j,i)
          SCAT%depol_ratio => depol(i)
          SCAT%tau => tau(:,i,j)
          SCAT%ssa => ssa(:,i,j)
          SCAT%pmom => pmom(:,i,j,:,:)

          if ( verbose > 0 ) then
            print *,'brdf output1', SCAT%Surface%Base%VIO%VBRDF_Sup_Out%BS_DBOUNCE_BRDFUNC(1,1,1,1) 
            print *,'brdf output2', SCAT%Surface%Base%VIO%VBRDF_Sup_Out%BS_DBOUNCE_BRDFUNC(2,1,1,1) 
            print *,'brdf output3', SCAT%Surface%Base%VIO%VBRDF_Sup_Out%BS_DBOUNCE_BRDFUNC(3,1,1,1) 
          end if

          BRDF(j,i) = SCAT%Surface%Base%VIO%VBRDF_Sup_Out%BS_DBOUNCE_BRDFUNC(1,1,1,1) 
          BRDF_Q(j,i) = SCAT%Surface%Base%VIO%VBRDF_Sup_Out%BS_DBOUNCE_BRDFUNC(2,1,1,1) 
          BRDF_U(j,i) = SCAT%Surface%Base%VIO%VBRDF_Sup_Out%BS_DBOUNCE_BRDFUNC(3,1,1,1)     

         
          call VLIDORT_Run_Vector (SCAT, output, ier)


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

end subroutine VLIDORT_Vector_OCICX

!..........................................................................

subroutine VLIDORT_Scalar_CX (km, nch, nobs,channels, nstreams, plane_parallel, nMom, &
                   nPol, ROT, depol, tau, ssa, g, pmom, pe, he, te, U10m, V10m, &
                   mr, solar_zenith, relat_azymuth, sensor_zenith, &
                   MISSING,verbose,radiance_VL_SURF,reflectance_Vl_SURF, BRDF,rc)
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

  logical,          intent(in)  :: plane_parallel ! do plane parallel flag

  integer, target,  intent(in)  :: nMom             ! number of phase function moments    
  integer, target,  intent(in)  :: nPol  ! number of scattering matrix components   
  integer,          intent(in)  :: nstreams  ! number of half space streams                                
                                      
                  
  real*8, target,   intent(in)  :: channels(nch)    ! wavelengths [nm]

!                                                   ! --- Rayleigh Parameters ---
  real*8, target,   intent(in)  :: ROT(km,nobs,nch) ! rayleigh optical thickness
  real*8, target,   intent(in)  :: depol(nch)       ! rayleigh depolarization ratio used in phase matrix

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
  real*8,           intent(out) :: BRDF(1,nobs, nch)  
!                         ---  
  integer             :: i,j,n,p,ier

  type(VLIDORT_scat) :: SCAT
  type(VLIDORT_output_scalar)  :: output  

  rc = 0
  ier = 0

  ! set streams here
  SCAT%Surface%Base%NSTREAMS = nstreams
  SCAT%Surface%Base%DO_PLANE_PARALLEL = plane_parallel
  call VLIDORT_Init( SCAT%Surface%Base, km, rc)
  if ( rc /= 0 ) return

  SCAT%nMom    = nMom
  SCAT%NSTOKES = 1    


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
            print*,U10m(j),V10m(j),mr(i),solar_zenith (j),&
                                    sensor_zenith(j),relat_azymuth(j),scalar

          end if
          call VLIDORT_CoxMunk(SCAT%Surface,U10m(j),V10m(j),mr(i),solar_zenith (j),&
                                    sensor_zenith(j),relat_azymuth(j),scalar,rc)
          if ( rc /= 0 ) return

          SCAT%wavelength = channels(i)
          SCAT%rot => ROT(:,j,i)
          SCAT%depol_ratio => depol(i)
          SCAT%tau => tau(:,i,j)
          SCAT%ssa => ssa(:,i,j)
          SCAT%g => g(:,i,j)
          SCAT%pmom => pmom(:,i,j,:,:)

          if ( verbose > 0 ) then
            print *,'brdf output1', SCAT%Surface%Base%VIO%VBRDF_Sup_Out%BS_DBOUNCE_BRDFUNC(1,1,1,1) 
          end if

          BRDF(1,j,i) = SCAT%Surface%Base%VIO%VBRDF_Sup_Out%BS_DBOUNCE_BRDFUNC(1,1,1,1) 

          if ( verbose > 0 ) then
            print *, 'VLIDORT inputs'
            print*,channels(i),tau(1,i,j),ssa(1,i,j),g(1,i,j),pmom(1,i,j,1,1)
          end if

          call VLIDORT_Run_Scalar (SCAT, output, ier)

          radiance_VL_SURF(j,i)    = output%radiance
          reflectance_VL_SURF(j,i) = output%reflectance

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

subroutine VLIDORT_Vector_CX_SingleGeom (km, nch, nobs,channels, nstreams, plane_parallel, nMom, &
                   nPol, ROT, depol, tau, ssa, pmom, pe, he, te, U10m, V10m, &
                   mr, solar_zenith, relat_azymuth, sensor_zenith, &
                   MISSING,verbose,radiance_VL_SURF,reflectance_VL_SURF, &
                   Q, U, BRDF, BRDF_Q, BRDF_U, rc)
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

  logical,          intent(in)  :: plane_parallel ! do plane parallel flag

  integer, target,  intent(in)  :: nMom  ! number of moments 
  integer, target,  intent(in)  :: nPol  ! number of components    
  integer,          intent(in)  :: nstreams  ! number of half space streams

                  
  real*8, target,   intent(in)  :: channels(nch)    ! wavelengths [nm]

!                                                   ! --- Rayleigh Parameters ---
  real*8, target,   intent(in)  :: ROT(km,nobs,nch) ! rayleigh optical thickness
  real*8, target,   intent(in)  :: depol(nch)       ! rayleigh depolarization ratio used in phase matrix

!                                                   ! --- Mie Parameters ---
  real*8, target,   intent(in)  :: tau(km,nch,nobs) ! aerosol optical depth
  real*8, target,   intent(in)  :: ssa(km,nch,nobs) ! single scattering albedo
  real*8, target,   intent(in)  :: pmom(km,nch,nobs,nMom,nPol) !components of the scat phase matrix

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
  real*8,           intent(out) :: BRDF(nobs, nch)
  real*8,           intent(out) :: BRDF_Q(nobs, nch)
  real*8,           intent(out) :: BRDF_U(nobs, nch)  
  real*8,           intent(out) :: Q(nobs, nch)   ! Stokes parameter Q
  real*8,           intent(out) :: U(nobs, nch)   ! Stokes parameter U
!                               ---
  
  integer             :: i,j,n,p,ier 
  
  type(VLIDORT_scat) :: SCAT
  type(VLIDORT_output_vector)  :: output  

  rc = 0
  ier = 0

  ! set streams here
  SCAT%Surface%Base%NSTREAMS = nstreams 
  SCAT%Surface%Base%DO_PLANE_PARALLEL = plane_parallel
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
            print*,U10m(j),V10m(j),mr(i),solar_zenith (j),&
                                    sensor_zenith(j),relat_azymuth(j)
          end if
          scalar = .false.
          call VLIDORT_CoxMunk(SCAT%Surface,U10m(j),V10m(j),mr(i),solar_zenith (j),&
                                    sensor_zenith(j),relat_azymuth(j),scalar,rc)

          if ( verbose > 0 ) then
            print*, 'FINISHED COX MUNK'
          end if


          if ( rc /= 0 ) return

          SCAT%wavelength = channels(i)
          SCAT%rot => ROT(:,j,i)
          SCAT%depol_ratio => depol(i)
          SCAT%tau => tau(:,i,j)
          SCAT%ssa => ssa(:,i,j)
          SCAT%pmom => pmom(:,i,j,:,:)

          if ( verbose > 0 ) then
            print *,'brdf output1', SCAT%Surface%Base%VIO%VBRDF_Sup_Out%BS_DBOUNCE_BRDFUNC(1,1,1,1) 
            print *,'brdf output2', SCAT%Surface%Base%VIO%VBRDF_Sup_Out%BS_DBOUNCE_BRDFUNC(2,1,1,1) 
            print *,'brdf output3', SCAT%Surface%Base%VIO%VBRDF_Sup_Out%BS_DBOUNCE_BRDFUNC(3,1,1,1) 
          end if

          BRDF(j,i) = SCAT%Surface%Base%VIO%VBRDF_Sup_Out%BS_DBOUNCE_BRDFUNC(1,1,1,1) 
          BRDF_Q(j,i) = 0 
          BRDF_U(j,i) = 0     

         
          call VLIDORT_Run_Vector (SCAT, output, ier)


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

end subroutine VLIDORT_Vector_CX_SingleGeom

!..........................................................................

subroutine VLIDORT_Vector_CX_MultiGeom (km, nch, nobs, ngeom, channels, nstreams, plane_parallel, nMom, &
                   nPol, ROT, depol, tau, ssa, pmom, pe, he, te, U10m, V10m, &
                   mr, solar_zenith, relat_azymuth, sensor_zenith, &
                   MISSING,verbose,radiance_VL_SURF,reflectance_VL_SURF, &
                   Q, U, BRDF, BRDF_Q, BRDF_U, rc)
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
  integer,          intent(in)  :: ngeom  ! number of geometries

  logical,          intent(in)  :: plane_parallel ! do plane parallel flag

  integer, target,  intent(in)  :: nMom  ! number of moments 
  integer, target,  intent(in)  :: nPol  ! number of components    
  integer,          intent(in)  :: nstreams  ! number of half space streams
                  
  real*8, target,   intent(in)  :: channels(nch)    ! wavelengths [nm]

!                                                   ! --- Rayleigh Parameters ---
  real*8, target,   intent(in)  :: ROT(km,nobs,nch) ! rayleigh optical thickness
  real*8, target,   intent(in)  :: depol(nch)       ! rayleigh depolarization ratio used in phase matrix

!                                                   ! --- Mie Parameters ---
  real*8, target,   intent(in)  :: tau(km,nch,nobs) ! aerosol optical depth
  real*8, target,   intent(in)  :: ssa(km,nch,nobs) ! single scattering albedo
  real*8, target,   intent(in)  :: pmom(km,nch,nobs,nMom,nPol) !components of the scat phase matrix

  real*8, target,   intent(in)  :: MISSING          ! MISSING VALUE
  real*8, target,   intent(in)  :: pe(km+1,nobs)    ! pressure at layer edges [Pa]
  real*8, target,   intent(in)  :: he(km+1,nobs)    ! height above sea-level  [m]
  real*8, target,   intent(in)  :: te(km+1,nobs)    ! temperature at layer edges [K]

  real*8, target,   intent(in)  :: U10m(nobs)   ! Wind speed components [m/s]
  real*8, target,   intent(in)  :: V10m(nobs)    
  real*8, target,   intent(in)  :: mr(nch)       ! refractive index
                       
  real*8, target,   intent(in)  :: solar_zenith(nobs,ngeom)  
  real*8, target,   intent(in)  :: relat_azymuth(nobs,ngeom) 
  real*8, target,   intent(in)  :: sensor_zenith(nobs,ngeom) 

  integer,          intent(in)  :: verbose

! !OUTPUT PARAMETERS:

  real*8,           intent(out) :: radiance_VL_SURF(nobs,nch,ngeom)       ! TOA normalized radiance from VLIDORT
  integer,          intent(out) :: rc                          ! return code
  real*8,           intent(out) :: reflectance_VL_SURF(nobs, nch,ngeom)   ! TOA reflectance from VLIDORT
  real*8,           intent(out) :: BRDF(nobs, nch,ngeom)
  real*8,           intent(out) :: BRDF_Q(nobs, nch,ngeom)
  real*8,           intent(out) :: BRDF_U(nobs, nch,ngeom)  
  real*8,           intent(out) :: Q(nobs, nch,ngeom)   ! Stokes parameter Q
  real*8,           intent(out) :: U(nobs, nch,ngeom)   ! Stokes parameter U
!                               ---
  
  integer             :: i,j,n,p,ier,g 
  
  type(VLIDORT_scat_multigeom) :: SCAT
  type(VLIDORT_output_vector_multigeom)  :: output  

  rc = 0
  ier = 0

  ! set streams here
  SCAT%Surface%Base%NSTREAMS = nstreams 
  SCAT%Surface%Base%DO_PLANE_PARALLEL = plane_parallel
  SCAT%Surface%Base%N_USER_OBSGEOMS   = ngeom
  SCAT%Surface%Base%NBEAMS            = ngeom
  SCAT%Surface%Base%N_USER_STREAMS    = ngeom
  SCAT%Surface%Base%N_USER_RELAZMS    = ngeom  
  call VLIDORT_Init( SCAT%Surface%Base, km, rc)
  if ( rc /= 0 ) return

  SCAT%nMom = nMom
  SCAT%nPol = nPol
  SCAT%NSTOKES = 3
  if ( SCAT%NSTOKES  .GT. MAXSTOKES  )   return

  do j = 1, nobs
     
     ! Make sure albedo and angles are available
     ! -----------------------------------------
     if ( IS_MISSING(solar_zenith(j,1),MISSING)  .OR. & 
          IS_MISSING(sensor_zenith(j,1),MISSING) .OR. &
          IS_MISSING(relat_azymuth(j,1),MISSING)  )  then

        radiance_VL_SURF(j,:,:) = MISSING
        reflectance_VL_SURF(j,:,:) = MISSING
      
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

              radiance_VL_SURF(j,i,:) = MISSING
              reflectance_VL_SURF(j,i,:) = MISSING
              cycle
          end if

          if ( verbose > 0 ) then
            print*, 'DO COX MUNK'
            print*,U10m(j),V10m(j),mr(i),solar_zenith (j,1),&
                                    sensor_zenith(j,1),relat_azymuth(j,1)
          end if
          scalar = .false.
          call VLIDORT_CoxMunk(SCAT%Surface,U10m(j),V10m(j),mr(i),solar_zenith (j,:),&
                                    sensor_zenith(j,:),relat_azymuth(j,:),scalar,rc)

          if ( verbose > 0 ) then
            print*, 'FINISHED COX MUNK'
          end if


          if ( rc /= 0 ) return

          SCAT%wavelength = channels(i)
          SCAT%rot => ROT(:,j,i)
          SCAT%depol_ratio => depol(i)
          SCAT%tau => tau(:,i,j)
          SCAT%ssa => ssa(:,i,j)
          SCAT%pmom => pmom(:,i,j,:,:)

          if ( verbose > 0 ) then
            print *,'brdf output1', SCAT%Surface%Base%VIO%VBRDF_Sup_Out%BS_DBOUNCE_BRDFUNC(1,1,1,1) 
            print *,'brdf output2', SCAT%Surface%Base%VIO%VBRDF_Sup_Out%BS_DBOUNCE_BRDFUNC(2,1,1,1) 
            print *,'brdf output3', SCAT%Surface%Base%VIO%VBRDF_Sup_Out%BS_DBOUNCE_BRDFUNC(3,1,1,1) 
          end if
   
          do g=1,ngeom
            BRDF(j,i,g) = SCAT%Surface%Base%VIO%VBRDF_Sup_Out%BS_DBOUNCE_BRDFUNC(1,1,1,g) 
            BRDF_Q(j,i,g) = 0 
            BRDF_U(j,i,g) = 0     
          end do
         
          call VLIDORT_Run_Vector (SCAT, output, ier)


          radiance_VL_SURF(j,i,:)    = output%radiance
          reflectance_VL_SURF(j,i,:) = output%reflectance
          Q(j,i,:)                   = output%Q
          U(j,i,:)                   = output%U                

          if ( ier /= 0 ) then
              radiance_VL_SURF(j,i,:) = MISSING
              reflectance_VL_SURF(j,i,:) = MISSING                          
              cycle
          end if


        end do ! end loop over channels
     
        if ( verbose > 0 ) then
           if ( mod(j-1,1000) == 0 ) then
              print *, '<> VLIDORT Vector: ', nint(j*100./nobs), '%'
           end if
        end if

  end do ! Loop over obs

end subroutine VLIDORT_Vector_Cx_MultiGeom


!..........................................................................

subroutine VLIDORT_Vector_CX_Cloud_SingleGeom (km, nch, nobs,channels, nstreams, plane_parallel, nMom, &
                   nPol, ROT, depol, tau, ssa, pmom, tauI, ssaI, pmomI, tauL, ssaL, pmomL, &
                   pe, he, te, U10m, V10m, &
                   mr, solar_zenith, relat_azymuth, sensor_zenith, &
                   MISSING,verbose,radiance_VL_SURF,reflectance_VL_SURF, &
                   Q, U, BRDF, BRDF_Q, BRDF_U, rc)
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

  logical,          intent(in)  :: plane_parallel ! do plane parallel flag

  integer, target,  intent(in)  :: nMom  ! number of moments 
  integer, target,  intent(in)  :: nPol  ! number of components    
  integer,          intent(in)  :: nstreams  ! number of half space streams

                  
  real*8, target,   intent(in)  :: channels(nch)    ! wavelengths [nm]

!                                                   ! --- Rayleigh Parameters ---
  real*8, target,   intent(in)  :: ROT(km,nobs,nch) ! rayleigh optical thickness
  real*8, target,   intent(in)  :: depol(nch)       ! rayleigh depolarization ratio used in phase matrix

!                                                   ! --- Mie Parameters ---
  real*8, target,   intent(in)  :: tau(km,nch,nobs) ! aerosol optical depth
  real*8, target,   intent(in)  :: ssa(km,nch,nobs) ! single scattering albedo
  real*8, target,   intent(in)  :: pmom(km,nch,nobs,nMom,nPol) !components of the scat phase matrix

  real*8, target,   intent(in)  :: tauI(km,nch,nobs) ! ice cloud optical depth
  real*8, target,   intent(in)  :: ssaI(km,nch,nobs) ! ice cloud single scattering albedo
  real*8, target,   intent(in)  :: pmomI(km,nch,nobs,nMom,nPol) !components of the scat phase matrix


  real*8, target,   intent(in)  :: tauL(km,nch,nobs) ! liquid cloud optical depth
  real*8, target,   intent(in)  :: ssaL(km,nch,nobs) ! liquid cloud single scattering albedo
  real*8, target,   intent(in)  :: pmomL(km,nch,nobs,nMom,nPol) !components of the scat phase matrix

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
  real*8,           intent(out) :: BRDF(nobs, nch)
  real*8,           intent(out) :: BRDF_Q(nobs, nch)
  real*8,           intent(out) :: BRDF_U(nobs, nch)  
  real*8,           intent(out) :: Q(nobs, nch)   ! Stokes parameter Q
  real*8,           intent(out) :: U(nobs, nch)   ! Stokes parameter U
!                               ---
  
  integer             :: i,j,n,p,ier 
  
  type(VLIDORT_scat) :: SCAT
  type(VLIDORT_output_vector)  :: output  

  rc = 0
  ier = 0

  ! set streams here
  SCAT%Surface%Base%NSTREAMS = nstreams 
  SCAT%Surface%Base%DO_PLANE_PARALLEL = plane_parallel
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
            print*,U10m(j),V10m(j),mr(i),solar_zenith (j),&
                                    sensor_zenith(j),relat_azymuth(j)
          end if
          scalar = .false.
          call VLIDORT_CoxMunk(SCAT%Surface,U10m(j),V10m(j),mr(i),solar_zenith (j),&
                                    sensor_zenith(j),relat_azymuth(j),scalar,rc)

          if ( verbose > 0 ) then
            print*, 'FINISHED COX MUNK'
          end if


          if ( rc /= 0 ) return

          SCAT%wavelength = channels(i)
          SCAT%rot => ROT(:,j,i)
          SCAT%depol_ratio => depol(i)
          SCAT%tau => tau(:,i,j)
          SCAT%ssa => ssa(:,i,j)
          SCAT%pmom => pmom(:,i,j,:,:)
          SCAT%tauI => tauI(:,i,j)
          SCAT%ssaI => ssaI(:,i,j)
          SCAT%pmomI => pmomI(:,i,j,:,:)
          SCAT%tauL => tauL(:,i,j)
          SCAT%ssaL => ssaL(:,i,j)
          SCAT%pmomL => pmomL(:,i,j,:,:)


          if ( verbose > 0 ) then
            print *,'brdf output1', SCAT%Surface%Base%VIO%VBRDF_Sup_Out%BS_DBOUNCE_BRDFUNC(1,1,1,1) 
            print *,'brdf output2', SCAT%Surface%Base%VIO%VBRDF_Sup_Out%BS_DBOUNCE_BRDFUNC(2,1,1,1) 
            print *,'brdf output3', SCAT%Surface%Base%VIO%VBRDF_Sup_Out%BS_DBOUNCE_BRDFUNC(3,1,1,1) 
          end if

          BRDF(j,i) = SCAT%Surface%Base%VIO%VBRDF_Sup_Out%BS_DBOUNCE_BRDFUNC(1,1,1,1) 
          BRDF_Q(j,i) = 0 
          BRDF_U(j,i) = 0     

         
          call VLIDORT_Run_Vector (SCAT, output, ier)


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

end subroutine VLIDORT_Vector_CX_Cloud_SingleGeom

!..........................................................................

subroutine VLIDORT_Vector_CX_Cloud_MultiGeom (km, nch, nobs, ngeom, channels, nstreams, plane_parallel, nMom, &
                   nPol, ROT, depol, tau, ssa, pmom, tauI, ssaI, pmomI, tauL, ssaL, pmomL, &
                   pe, he, te, U10m, V10m, &
                   mr, solar_zenith, relat_azymuth, sensor_zenith, &
                   MISSING,verbose,radiance_VL_SURF,reflectance_VL_SURF, &
                   Q, U, BRDF, BRDF_Q, BRDF_U, rc)
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
  integer,          intent(in)  :: ngeom  ! number of geometries

  logical,          intent(in)  :: plane_parallel ! do plane parallel flag

  integer, target,  intent(in)  :: nMom  ! number of moments 
  integer, target,  intent(in)  :: nPol  ! number of components    
  integer,          intent(in)  :: nstreams  ! number of half space streams

                  
  real*8, target,   intent(in)  :: channels(nch)    ! wavelengths [nm]

!                                                   ! --- Rayleigh Parameters ---
  real*8, target,   intent(in)  :: ROT(km,nobs,nch) ! rayleigh optical thickness
  real*8, target,   intent(in)  :: depol(nch)       ! rayleigh depolarization ratio used in phase matrix

!                                                   ! --- Mie Parameters ---
  real*8, target,   intent(in)  :: tau(km,nch,nobs) ! aerosol optical depth
  real*8, target,   intent(in)  :: ssa(km,nch,nobs) ! single scattering albedo
  real*8, target,   intent(in)  :: pmom(km,nch,nobs,nMom,nPol) !components of the scat phase matrix

  real*8, target,   intent(in)  :: tauI(km,nch,nobs) ! ice cloud optical depth
  real*8, target,   intent(in)  :: ssaI(km,nch,nobs) ! ice cloud single scattering albedo
  real*8, target,   intent(in)  :: pmomI(km,nch,nobs,nMom,nPol) !components of the scat phase matrix


  real*8, target,   intent(in)  :: tauL(km,nch,nobs) ! liquid cloud optical depth
  real*8, target,   intent(in)  :: ssaL(km,nch,nobs) ! liquid cloud single scattering albedo
  real*8, target,   intent(in)  :: pmomL(km,nch,nobs,nMom,nPol) !components of the scat phase matrix

  real*8, target,   intent(in)  :: MISSING          ! MISSING VALUE
  real*8, target,   intent(in)  :: pe(km+1,nobs)    ! pressure at layer edges [Pa]
  real*8, target,   intent(in)  :: he(km+1,nobs)    ! height above sea-level  [m]
  real*8, target,   intent(in)  :: te(km+1,nobs)    ! temperature at layer edges [K]

  real*8, target,   intent(in)  :: U10m(nobs)   ! Wind speed components [m/s]
  real*8, target,   intent(in)  :: V10m(nobs)    
  real*8, target,   intent(in)  :: mr(nch)       ! refractive index
                       
  real*8, target,   intent(in)  :: solar_zenith(nobs,ngeom)  
  real*8, target,   intent(in)  :: relat_azymuth(nobs,ngeom) 
  real*8, target,   intent(in)  :: sensor_zenith(nobs,ngeom) 

  integer,          intent(in)  :: verbose

! !OUTPUT PARAMETERS:

  real*8,           intent(out) :: radiance_VL_SURF(nobs,nch,ngeom)       ! TOA normalized radiance from VLIDORT
  integer,          intent(out) :: rc                          ! return code
  real*8,           intent(out) :: reflectance_VL_SURF(nobs, nch,ngeom)   ! TOA reflectance from VLIDORT
  real*8,           intent(out) :: BRDF(nobs, nch,ngeom)
  real*8,           intent(out) :: BRDF_Q(nobs, nch,ngeom)
  real*8,           intent(out) :: BRDF_U(nobs, nch,ngeom)  
  real*8,           intent(out) :: Q(nobs, nch,ngeom)   ! Stokes parameter Q
  real*8,           intent(out) :: U(nobs, nch,ngeom)   ! Stokes parameter U
!                               ---
  
  integer             :: i,j,n,p,ier,g 
  
  type(VLIDORT_scat_multigeom) :: SCAT
  type(VLIDORT_output_vector_multigeom)  :: output  

  rc = 0
  ier = 0

  ! set streams here
  SCAT%Surface%Base%NSTREAMS = nstreams 
  SCAT%Surface%Base%DO_PLANE_PARALLEL = plane_parallel
  SCAT%Surface%Base%N_USER_OBSGEOMS   = ngeom
  SCAT%Surface%Base%NBEAMS            = ngeom
  SCAT%Surface%Base%N_USER_STREAMS    = ngeom
  SCAT%Surface%Base%N_USER_RELAZMS    = ngeom
  call VLIDORT_Init( SCAT%Surface%Base, km, rc)
  if ( rc /= 0 ) return

  SCAT%nMom = nMom
  SCAT%nPol = nPol
  SCAT%NSTOKES = 3
  if ( SCAT%NSTOKES  .GT. MAXSTOKES  )   return

  do j = 1, nobs
     
     ! Make sure albedo and angles are available
     ! -----------------------------------------
     if ( IS_MISSING(solar_zenith(j,1),MISSING)  .OR. & 
          IS_MISSING(sensor_zenith(j,1),MISSING) .OR. &
          IS_MISSING(relat_azymuth(j,1),MISSING)  )  then

        radiance_VL_SURF(j,:,:) = MISSING
        reflectance_VL_SURF(j,:,:) = MISSING
      
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

              radiance_VL_SURF(j,i,:) = MISSING
              reflectance_VL_SURF(j,i,:) = MISSING
              cycle
          end if

          if ( verbose > 0 ) then
            print*, 'DO COX MUNK'
            print*,U10m(j),V10m(j),mr(i),solar_zenith (j,1),&
                                    sensor_zenith(j,1),relat_azymuth(j,1)
          end if
          scalar = .false.
          call VLIDORT_CoxMunk(SCAT%Surface,U10m(j),V10m(j),mr(i),solar_zenith (j,:),&
                                    sensor_zenith(j,:),relat_azymuth(j,:),scalar,rc)

          if ( verbose > 0 ) then
            print*, 'FINISHED COX MUNK'
          end if


          if ( rc /= 0 ) return

          SCAT%wavelength = channels(i)
          SCAT%rot => ROT(:,j,i)
          SCAT%depol_ratio => depol(i)
          SCAT%tau => tau(:,i,j)
          SCAT%ssa => ssa(:,i,j)
          SCAT%pmom => pmom(:,i,j,:,:)
          SCAT%tauI => tauI(:,i,j)
          SCAT%ssaI => ssaI(:,i,j)
          SCAT%pmomI => pmomI(:,i,j,:,:)
          SCAT%tauL => tauL(:,i,j)
          SCAT%ssaL => ssaL(:,i,j)
          SCAT%pmomL => pmomL(:,i,j,:,:)


          if ( verbose > 0 ) then
            print *,'brdf output1', SCAT%Surface%Base%VIO%VBRDF_Sup_Out%BS_DBOUNCE_BRDFUNC(1,1,1,1) 
            print *,'brdf output2', SCAT%Surface%Base%VIO%VBRDF_Sup_Out%BS_DBOUNCE_BRDFUNC(2,1,1,1) 
            print *,'brdf output3', SCAT%Surface%Base%VIO%VBRDF_Sup_Out%BS_DBOUNCE_BRDFUNC(3,1,1,1) 
          end if

          do g=1,ngeom
            BRDF(j,i,g) = SCAT%Surface%Base%VIO%VBRDF_Sup_Out%BS_DBOUNCE_BRDFUNC(1,1,1,g) 
            BRDF_Q(j,i,g) = 0 
            BRDF_U(j,i,g) = 0     
          end do
         
          call VLIDORT_Run_Vector (SCAT, output, ier)


          radiance_VL_SURF(j,i,:)    = output%radiance
          reflectance_VL_SURF(j,i,:) = output%reflectance
          Q(j,i,:)                   = output%Q
          U(j,i,:)                   = output%U                

          if ( ier /= 0 ) then
              radiance_VL_SURF(j,i,:) = MISSING
              reflectance_VL_SURF(j,i,:) = MISSING                          
              cycle
          end if


        end do ! end loop over channels
     
        if ( verbose > 0 ) then
           if ( mod(j-1,1000) == 0 ) then
              print *, '<> VLIDORT Vector: ', nint(j*100./nobs), '%'
           end if
        end if

  end do ! Loop over obs

end subroutine VLIDORT_Vector_CX_Cloud_MultiGeom

!..........................................................................



subroutine VLIDORT_Vector_CX_NOBM (km, nch, nobs,channels, nstreams, plane_parallel, nMom, &
                   nPol, ROT, depol, tau, ssa, pmom, pe, he, te, U10m, V10m, &
                   mr, sleave, sleave_iso, sleave_adjust, & 
                   solar_zenith, relat_azymuth, sensor_zenith, &
                   MISSING,verbose,radiance_VL_SURF,reflectance_VL_SURF, &
                   Q, U, BRDF, BRDF_Q, BRDF_U, rc, ADJUSTED_SLEAVE)
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

  logical,          intent(in)  :: plane_parallel ! do plane parallel flag

  integer, target,  intent(in)  :: nMom  ! number of moments 
  integer, target,  intent(in)  :: nPol  ! number of components                               
  integer,          intent(in)  :: nstreams  ! number of half space streams
                 
  real*8, target,   intent(in)  :: channels(nch)    ! wavelengths [nm]

!                                                   ! --- Rayleigh Parameters ---
  real*8, target,   intent(in)  :: ROT(km,nobs,nch) ! rayleigh optical thickness
  real*8, target,   intent(in)  :: depol(nch)       ! rayleigh depolarization ratio used in phase matrix

!                                                   ! --- Mie Parameters ---
  real*8, target,   intent(in)  :: tau(km,nch,nobs) ! aerosol optical depth
  real*8, target,   intent(in)  :: ssa(km,nch,nobs) ! single scattering albedo
  real*8, target,   intent(in)  :: pmom(km,nch,nobs,nMom,nPol) !components of the scat phase matrix

  real*8, target,   intent(in)  :: MISSING          ! MISSING VALUE
  real*8, target,   intent(in)  :: pe(km+1,nobs)    ! pressure at layer edges [Pa]
  real*8, target,   intent(in)  :: he(km+1,nobs)    ! height above sea-level  [m]
  real*8, target,   intent(in)  :: te(km+1,nobs)    ! temperature at layer edges [K]

  real*8, target,   intent(in)  :: U10m(nobs)   ! Wind speed components [m/s]
  real*8, target,   intent(in)  :: V10m(nobs)    
  real*8, target,   intent(in)  :: mr(nch)       ! refractive index
  real*8, target,   intent(in)  :: sleave(nch,nobs)       ! sun normalized water leaving radiance from NOBM
  logical,          intent(in)  :: sleave_iso     ! do isotropic water leaving radiance 
  logical, target,  intent(in)  :: sleave_adjust  ! do water leaving adjustment for transmittance.
                                                          ! described and used in Gregg & Rousseaux 2017   
                       
  real*8, target,   intent(in)  :: solar_zenith(nobs)  
  real*8, target,   intent(in)  :: relat_azymuth(nobs) 
  real*8, target,   intent(in)  :: sensor_zenith(nobs) 

  integer,          intent(in)  :: verbose

! !OUTPUT PARAMETERS:

  real*8,           intent(out) :: radiance_VL_SURF(nobs,nch)       ! TOA normalized radiance from VLIDORT
  integer,          intent(out) :: rc                          ! return code
  real*8,           intent(out) :: reflectance_VL_SURF(nobs, nch)   ! TOA reflectance from VLIDORT
  real*8,           intent(out) :: BRDF(nobs, nch)
  real*8,           intent(out) :: BRDF_Q(nobs, nch)
  real*8,           intent(out) :: BRDF_U(nobs, nch)  
  real*8,           intent(out) :: Q(nobs, nch)   ! Stokes parameter Q
  real*8,           intent(out) :: U(nobs, nch)   ! Stokes parameter U
  real*8, optional, intent(out) :: ADJUSTED_SLEAVE(nobs,nch)
!                               ---
  
  integer             :: i,j,n,p,ier 
  
  type(VLIDORT_scat) :: SCAT
  type(VLIDORT_output_vector)  :: output  

  rc = 0
  ier = 0
 
  ! set streams here
  SCAT%Surface%Base%NSTREAMS = nstreams
  SCAT%Surface%Base%DO_PLANE_PARALLEL = plane_parallel
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
            print*,U10m(j),V10m(j),mr(i),solar_zenith (j),&
                                    sensor_zenith(j),relat_azymuth(j)
          end if
          scalar = .false.
          call VLIDORT_CoxMunk(SCAT%Surface,U10m(j),V10m(j),mr(i),solar_zenith (j),&
                                    sensor_zenith(j),relat_azymuth(j),scalar,rc)

          if ( rc /= 0 ) return
          ! must be called second
          if (sleave_iso) then
             call VLIDORT_NOBM_ISO(SCAT%Surface,sleave(i,j),solar_zenith(j),&
                                    sensor_zenith(j),relat_azymuth(j),&
                                    scalar,rc,sleave_adjust)
          else
             call VLIDORT_NOBM(SCAT%Surface,sleave(i,j),solar_zenith(j),&
                                    sensor_zenith(j),relat_azymuth(j),&
                                    scalar,rc,sleave_adjust)
          end if

          if ( verbose > 0 ) then
            print*, 'FINISHED COX MUNK'
          end if


          if ( rc /= 0 ) return

          SCAT%wavelength = channels(i)
          SCAT%rot => ROT(:,j,i)
          SCAT%depol_ratio => depol(i)
          SCAT%tau => tau(:,i,j)
          SCAT%ssa => ssa(:,i,j)
          SCAT%pmom => pmom(:,i,j,:,:)

          if ( verbose > 0 ) then
            print *,'brdf output1', SCAT%Surface%Base%VIO%VBRDF_Sup_Out%BS_DBOUNCE_BRDFUNC(1,1,1,1) 
            print *,'brdf output2', SCAT%Surface%Base%VIO%VBRDF_Sup_Out%BS_DBOUNCE_BRDFUNC(2,1,1,1) 
            print *,'brdf output3', SCAT%Surface%Base%VIO%VBRDF_Sup_Out%BS_DBOUNCE_BRDFUNC(3,1,1,1) 
          end if

          BRDF(j,i) = SCAT%Surface%Base%VIO%VBRDF_Sup_Out%BS_DBOUNCE_BRDFUNC(1,1,1,1) 
          BRDF_Q(j,i) = SCAT%Surface%Base%VIO%VBRDF_Sup_Out%BS_DBOUNCE_BRDFUNC(2,1,1,1) 
          BRDF_U(j,i) = SCAT%Surface%Base%VIO%VBRDF_Sup_Out%BS_DBOUNCE_BRDFUNC(3,1,1,1)     

         
          call VLIDORT_Run_Vector (SCAT, output, ier)


          radiance_VL_SURF(j,i)    = output%radiance
          reflectance_VL_SURF(j,i) = output%reflectance
          Q(j,i)                   = output%Q
          U(j,i)                   = output%U                

          if (sleave_adjust) then
              ADJUSTED_SLEAVE(j,i) = output%ADJUSTED_SLEAVE
          else
              if (present(ADJUSTED_SLEAVE)) then
                ADJUSTED_SLEAVE(j,i) = 0.0
              end if 
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

end subroutine VLIDORT_Vector_Cx_NOBM


subroutine VLIDORT_Scalar_GissCX (km, nch, nobs,channels, nstreams, plane_parallel, nMom, &
                   nPol, ROT, depol, tau, ssa, g, pmom, pe, he, te, U10m, V10m, &
                   mr, solar_zenith, relat_azymuth, sensor_zenith, &
                   MISSING,verbose,radiance_VL_SURF,reflectance_Vl_SURF, BRDF,rc)
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

  logical,          intent(in)  :: plane_parallel ! do plane parallel flag

  integer, target,  intent(in)  :: nMom             ! number of phase function moments    
  integer, target,  intent(in)  :: nPol  ! number of scattering matrix components                               
  integer,          intent(in)  :: nstreams  ! number of half space streams                                      
                  
  real*8, target,   intent(in)  :: channels(nch)    ! wavelengths [nm]

!                                                   ! --- Rayleigh Parameters ---
  real*8, target,   intent(in)  :: ROT(km,nobs,nch) ! rayleigh optical thickness
  real*8, target,   intent(in)  :: depol(nch)       ! rayleigh depolarization ratio used in phase matrix

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
  real*8,           intent(out) :: BRDF(1,nobs, nch)  
!                         ---  
  integer             :: i,j,n,p,ier

  type(VLIDORT_scat) :: SCAT
  type(VLIDORT_output_scalar)  :: output  

  rc = 0
  ier = 0

  ! Set streams here
  SCAT%Surface%Base%NSTREAMS = nstreams
  SCAT%Surface%Base%DO_PLANE_PARALLEL = plane_parallel
  call VLIDORT_Init( SCAT%Surface%Base, km, rc)
  if ( rc /= 0 ) return

  SCAT%nMom    = nMom
  SCAT%NSTOKES = 1    


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
            print*,U10m(j),V10m(j),mr(i),solar_zenith (j),&
                                    sensor_zenith(j),relat_azymuth(j),scalar

          end if
          call VLIDORT_GissCoxMunk(SCAT%Surface,U10m(j),V10m(j),mr(i),solar_zenith (j),&
                                    sensor_zenith(j),relat_azymuth(j),scalar,rc)
          if ( rc /= 0 ) return

          SCAT%wavelength = channels(i)
          SCAT%rot => ROT(:,j,i)
          SCAT%depol_ratio => depol(i)
          SCAT%tau => tau(:,i,j)
          SCAT%ssa => ssa(:,i,j)
          SCAT%g => g(:,i,j)
          SCAT%pmom => pmom(:,i,j,:,:)

          if ( verbose > 0 ) then
            print *,'brdf output1', SCAT%Surface%Base%VIO%VBRDF_Sup_Out%BS_DBOUNCE_BRDFUNC(1,1,1,1) 
          end if

          BRDF(1,j,i) = SCAT%Surface%Base%VIO%VBRDF_Sup_Out%BS_DBOUNCE_BRDFUNC(1,1,1,1) 

          if ( verbose > 0 ) then
            print *, 'VLIDORT inputs'
            print*,channels(i),tau(1,i,j),ssa(1,i,j),g(1,i,j),pmom(1,i,j,1,1)
          end if

          call VLIDORT_Run_Scalar (SCAT, output, ier)

          radiance_VL_SURF(j,i)    = output%radiance
          reflectance_VL_SURF(j,i) = output%reflectance

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

end subroutine VLIDORT_Scalar_GissCX

!..........................................................................

subroutine VLIDORT_Vector_GissCX (km, nch, nobs,channels, nstreams, plane_parallel, nMom, &
                   nPol, ROT, depol, tau, ssa, pmom, pe, he, te, U10m, V10m, &
                   mr, solar_zenith, relat_azymuth, sensor_zenith, &
                   MISSING,verbose,radiance_VL_SURF,reflectance_VL_SURF, &
                   Q, U, BRDF, BRDF_Q, BRDF_U, rc)
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

  logical,          intent(in)  :: plane_parallel

  integer, target,  intent(in)  :: nMom  ! number of moments 
  integer, target,  intent(in)  :: nPol  ! number of components                               
  integer,          intent(in)  :: nstreams  ! number of half space streams
                  
  real*8, target,   intent(in)  :: channels(nch)    ! wavelengths [nm]

!                                                   ! --- Rayleigh Parameters ---
  real*8, target,   intent(in)  :: ROT(km,nobs,nch) ! rayleigh optical thickness
  real*8, target,   intent(in)  :: depol(nch)       ! rayleigh depolarization ratio used in phase matrix

!                                                   ! --- Mie Parameters ---
  real*8, target,   intent(in)  :: tau(km,nch,nobs) ! aerosol optical depth
  real*8, target,   intent(in)  :: ssa(km,nch,nobs) ! single scattering albedo
  real*8, target,   intent(in)  :: pmom(km,nch,nobs,nMom,nPol) !components of the scat phase matrix

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
  real*8,           intent(out) :: BRDF(nobs, nch)
  real*8,           intent(out) :: BRDF_Q(nobs, nch)
  real*8,           intent(out) :: BRDF_U(nobs, nch)  
  real*8,           intent(out) :: Q(nobs, nch)   ! Stokes parameter Q
  real*8,           intent(out) :: U(nobs, nch)   ! Stokes parameter U
!                               ---
  
  integer             :: i,j,n,p,ier 
  
  type(VLIDORT_scat) :: SCAT
  type(VLIDORT_output_vector)  :: output  

  rc = 0
  ier = 0
 
  ! set streams here
  SCAT%Surface%BASE%NSTREAMS = nstreams
  SCAT%Surface%BASE%DO_PLANE_PARALLEL = plane_parallel
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
            print*,U10m(j),V10m(j),mr(i),solar_zenith (j),&
                                    sensor_zenith(j),relat_azymuth(j)
          end if
          scalar = .false.
          call VLIDORT_GissCoxMunk(SCAT%Surface,U10m(j),V10m(j),mr(i),solar_zenith (j),&
                                    sensor_zenith(j),relat_azymuth(j),scalar,rc)

          if ( verbose > 0 ) then
            print*, 'FINISHED COX MUNK'
          end if


          if ( rc /= 0 ) return

          SCAT%wavelength = channels(i)
          SCAT%rot => ROT(:,j,i)
          SCAT%depol_ratio => depol(i)
          SCAT%tau => tau(:,i,j)
          SCAT%ssa => ssa(:,i,j)
          SCAT%pmom => pmom(:,i,j,:,:)

          if ( verbose > 0 ) then
            print *,'brdf output1', SCAT%Surface%Base%VIO%VBRDF_Sup_Out%BS_DBOUNCE_BRDFUNC(1,1,1,1) 
            print *,'brdf output2', SCAT%Surface%Base%VIO%VBRDF_Sup_Out%BS_DBOUNCE_BRDFUNC(2,1,1,1) 
            print *,'brdf output3', SCAT%Surface%Base%VIO%VBRDF_Sup_Out%BS_DBOUNCE_BRDFUNC(3,1,1,1) 
          end if

          BRDF(j,i) = SCAT%Surface%Base%VIO%VBRDF_Sup_Out%BS_DBOUNCE_BRDFUNC(1,1,1,1) 
          BRDF_Q(j,i) = SCAT%Surface%Base%VIO%VBRDF_Sup_Out%BS_DBOUNCE_BRDFUNC(2,1,1,1) 
          BRDF_U(j,i) = SCAT%Surface%Base%VIO%VBRDF_Sup_Out%BS_DBOUNCE_BRDFUNC(3,1,1,1)     

         
          call VLIDORT_Run_Vector (SCAT, output, ier)


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

end subroutine VLIDORT_Vector_GissCx

!..........................................................................

subroutine VLIDORT_Scalar_GissCX_Cloud (km, nch, nobs,channels, nstreams, plane_parallel, nMom, &
                   nPol, ROT, depol, tau, ssa, g, pmom, tauI, ssaI, gI, pmomI, tauL, ssaL, gL, pmomL, &
                   pe, he, te, U10m, V10m, &
                   mr, solar_zenith, relat_azymuth, sensor_zenith, &
                   MISSING,verbose,radiance_VL_SURF,reflectance_Vl_SURF, BRDF,rc)
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

  logical,          intent(in)  :: plane_parallel ! do plane parallel flag

  integer, target,  intent(in)  :: nMom             ! number of phase function moments    
  integer, target,  intent(in)  :: nPol  ! number of scattering matrix components                               
  integer,          intent(in)  :: nstreams  ! number of half space streams                                      
                  
  real*8, target,   intent(in)  :: channels(nch)    ! wavelengths [nm]

!                                                   ! --- Rayleigh Parameters ---
  real*8, target,   intent(in)  :: ROT(km,nobs,nch) ! rayleigh optical thickness
  real*8, target,   intent(in)  :: depol(nch)       ! rayleigh depolarization ratio used in phase matrix

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
  real*8,           intent(out) :: BRDF(nobs, nch)  
!                         ---  
  integer             :: i,j,n,p,ier

  type(VLIDORT_scat) :: SCAT
  type(VLIDORT_output_scalar)  :: output  

  rc = 0
  ier = 0

  ! set streams here
  SCAT%Surface%Base%NSTREAMS = nstreams
  SCAT%Surface%Base%DO_PLANE_PARALLEL = plane_parallel
  call VLIDORT_Init( SCAT%Surface%Base, km, rc)
  if ( rc /= 0 ) return

  SCAT%nMom    = nMom
  SCAT%NSTOKES = 1    


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
            print*,U10m(j),V10m(j),mr(i),solar_zenith (j),&
                                    sensor_zenith(j),relat_azymuth(j),scalar

          end if
          call VLIDORT_GissCoxMunk(SCAT%Surface,U10m(j),V10m(j),mr(i),solar_zenith (j),&
                                    sensor_zenith(j),relat_azymuth(j),scalar,rc)
          if ( rc /= 0 ) return

          SCAT%wavelength = channels(i)
          SCAT%rot => ROT(:,j,i)
          SCAT%depol_ratio => depol(i)
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

          if ( verbose > 0 ) then
            print *,'brdf output1', SCAT%Surface%Base%VIO%VBRDF_Sup_Out%BS_DBOUNCE_BRDFUNC(1,1,1,1) 
          end if

          BRDF(j,i) = SCAT%Surface%Base%VIO%VBRDF_Sup_Out%BS_DBOUNCE_BRDFUNC(1,1,1,1) 

          if ( verbose > 0 ) then
            print *, 'VLIDORT inputs'
            print*,channels(i),tau(1,i,j),ssa(1,i,j),g(1,i,j),pmom(1,i,j,1,1)
          end if

          call VLIDORT_Run_Scalar_Cloud (SCAT, output, ier)

          radiance_VL_SURF(j,i)    = output%radiance
          reflectance_VL_SURF(j,i) = output%reflectance

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

end subroutine VLIDORT_Scalar_GissCX_Cloud

!..........................................................................

subroutine VLIDORT_Scalar_OCIGissCX_Cloud (km, nch, nobs,channels, nstreams, plane_parallel, nMom, &
                   nPol, ROT, depol, alpha, tau, ssa, g, pmom, tauI, ssaI, gI, pmomI, tauL, ssaL, gL, pmomL, &
                   pe, he, te, U10m, V10m, &
                   mr, solar_zenith, relat_azymuth, sensor_zenith, flux_factor, &
                   MISSING,verbose,radiance_VL_SURF,reflectance_Vl_SURF, BRDF,rc)
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

  logical,          intent(in)  :: plane_parallel ! do plane parallel flag

  integer, target,  intent(in)  :: nMom             ! number of phase function moments    
  integer, target,  intent(in)  :: nPol  ! number of scattering matrix components                               
  integer,          intent(in)  :: nstreams  ! number of half space streams                                      
                  
  real*8, target,   intent(in)  :: channels(nch)    ! wavelengths [nm]

!                                                   ! --- Rayleigh Parameters ---
  real*8, target,   intent(in)  :: ROT(km,nobs,nch) ! rayleigh optical thickness
  real*8, target,   intent(in)  :: depol(nch)       ! rayleigh depolarization ratio used in phase matrix

!                                                   ! --- Trace Gas Absorption ---
  real*8, target,   intent(in)  :: alpha(km,nobs,nch) ! trace gas absoprtion optical thickness

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

  real*8, target,   intent(in)  :: U10m(nobs)   ! Wind speed components [m/s]
  real*8, target,   intent(in)  :: V10m(nobs)    
  real*8, target,   intent(in)  :: mr(nch)       ! refractive index
  
  real*8, target,   intent(in)  :: MISSING          ! MISSING VALUE                                      
  real*8, target,   intent(in)  :: solar_zenith(nobs)  
  real*8, target,   intent(in)  :: relat_azymuth(nobs) 
  real*8, target,   intent(in)  :: sensor_zenith(nobs) 

  real*8,           intent(in)  :: flux_factor(nch,nobs) ! solar flux (F0)  
  integer,          intent(in)  :: verbose

! !OUTPUT PARAMETERS:

  real*8,           intent(out) :: radiance_VL_SURF(nobs,nch)       ! TOA normalized radiance from VLIDORT
  integer,          intent(out) :: rc                          ! return code
  real*8,           intent(out) :: reflectance_VL_SURF(nobs, nch)   ! TOA reflectance from VLIDORT
  real*8,           intent(out) :: BRDF(nobs, nch)  
!                         ---  
  integer             :: i,j,n,p,ier

  type(VLIDORT_scat) :: SCAT
  type(VLIDORT_output_scalar)  :: output  

  rc = 0
  ier = 0

  ! set streams here
  SCAT%Surface%Base%NSTREAMS = nstreams
  SCAT%Surface%Base%DO_PLANE_PARALLEL = plane_parallel
  call VLIDORT_Init( SCAT%Surface%Base, km, rc)
  if ( rc /= 0 ) return

  SCAT%nMom    = nMom
  SCAT%NSTOKES = 1    


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
          ! set solar flux
          SCAT%Surface%Base%VIO%VLIDORT_FixIn%SunRays%TS_FLUX_FACTOR = flux_factor(i,j)
       
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
            print*,U10m(j),V10m(j),mr(i),solar_zenith (j),&
                                    sensor_zenith(j),relat_azymuth(j),scalar

          end if
          call VLIDORT_OCIGissCoxMunk(SCAT%Surface,U10m(j),V10m(j),mr(i),solar_zenith (j),&
                                    sensor_zenith(j),relat_azymuth(j),scalar,rc)
          if ( rc /= 0 ) return

          SCAT%wavelength = channels(i)
          SCAT%rot => ROT(:,j,i)
          SCAT%depol_ratio => depol(i)
          SCAT%alpha => alpha(:,j,i)
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

          if ( verbose > 0 ) then
            print *,'brdf output1', SCAT%Surface%Base%VIO%VBRDF_Sup_Out%BS_DBOUNCE_BRDFUNC(1,1,1,1) 
          end if

          BRDF(j,i) = SCAT%Surface%Base%VIO%VBRDF_Sup_Out%BS_DBOUNCE_BRDFUNC(1,1,1,1) 

          if ( verbose > 0 ) then
            print *, 'VLIDORT inputs'
            print*,channels(i),tau(1,i,j),ssa(1,i,j),g(1,i,j),pmom(1,i,j,1,1)
          end if

          call VLIDORT_Run_Scalar_Cloud (SCAT, output, ier)

          radiance_VL_SURF(j,i)    = output%radiance
          reflectance_VL_SURF(j,i) = output%reflectance

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

end subroutine VLIDORT_Scalar_OCIGissCX_Cloud


subroutine VLIDORT_Vector_GissCX_Cloud (km, nch, nobs,channels, nstreams, plane_parallel, nMom, &
                   nPol, ROT, depol, tau, ssa, pmom, tauI, ssaI, pmomI, tauL, ssaL, pmomL, &
                   pe, he, te, U10m, V10m, &
                   mr, solar_zenith, relat_azymuth, sensor_zenith, &
                   MISSING,verbose,radiance_VL_SURF,reflectance_VL_SURF, Q, U, BRDF, BRDF_Q, BRDF_U,rc)
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

  logical,          intent(in)  :: plane_parallel ! do plane parallel flag

  integer, target,  intent(in)  :: nMom  ! number of moments 
  integer, target,  intent(in)  :: nPol  ! number of components                               
  integer,          intent(in)  :: nstreams  ! number of half space streams
                  
  real*8, target,   intent(in)  :: channels(nch)    ! wavelengths [nm]

!                                                   ! --- Rayleigh Parameters ---
  real*8, target,   intent(in)  :: ROT(km,nobs,nch) ! rayleigh optical thickness
  real*8, target,   intent(in)  :: depol(nch)       ! rayleigh depolarization ratio used in phase matrix

!                                                   ! --- Mie Parameters ---
  real*8, target,   intent(in)  :: tau(km,nch,nobs) ! aerosol optical depth
  real*8, target,   intent(in)  :: ssa(km,nch,nobs) ! single scattering albedo
  real*8, target,   intent(in)  :: pmom(km,nch,nobs,nMom,nPol) !components of the scat phase matrix

  real*8, target,   intent(in)  :: tauI(km,nch,nobs) ! ice cloud optical depth
  real*8, target,   intent(in)  :: ssaI(km,nch,nobs) ! ice cloud single scattering albedo    
  real*8, target,   intent(in)  :: pmomI(km,nch,nobs,nMom,nPol) !components of the scat phase matrix


  real*8, target,   intent(in)  :: tauL(km,nch,nobs) ! liquid cloud optical depth
  real*8, target,   intent(in)  :: ssaL(km,nch,nobs) ! liquid cloud single scattering albedo  
  real*8, target,   intent(in)  :: pmomL(km,nch,nobs,nMom,nPol) !components of the scat phase matrix

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
  real*8,           intent(out) :: BRDF(nobs, nch)  
  real*8,           intent(out) :: BRDF_Q(nobs, nch)  
  real*8,           intent(out) :: BRDF_U(nobs, nch)      
  real*8,           intent(out) :: Q(nobs, nch)   ! Stokes parameter Q
  real*8,           intent(out) :: U(nobs, nch)   ! Stokes parameter U
!                               ---
  
  integer             :: i,j,n,p,ier 
  
  type(VLIDORT_scat) :: SCAT
  type(VLIDORT_output_vector)  :: output  

  rc = 0
  ier = 0
 
  ! set streams here
  SCAT%Surface%Base%NSTREAMS = nstreams
  SCAT%Surface%Base%DO_PLANE_PARALLEL = plane_parallel
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
            print*,U10m(j),V10m(j),mr(i),solar_zenith (j),&
                                    sensor_zenith(j),relat_azymuth(j)
          end if
          scalar = .false.
          call VLIDORT_GissCoxMunk(SCAT%Surface,U10m(j),V10m(j),mr(i),solar_zenith (j),&
                                    sensor_zenith(j),relat_azymuth(j),scalar,rc)

          if ( verbose > 0 ) then
            print*, 'FINISHED COX MUNK'
          end if


          if ( rc /= 0 ) return

          SCAT%wavelength = channels(i)
          SCAT%rot => ROT(:,j,i)
          SCAT%depol_ratio => depol(i)
          SCAT%tau => tau(:,i,j)
          SCAT%ssa => ssa(:,i,j)
          SCAT%pmom => pmom(:,i,j,:,:)
          SCAT%tauI => tauI(:,i,j)
          SCAT%ssaI => ssaI(:,i,j)
          SCAT%pmomI => pmomI(:,i,j,:,:)
          SCAT%tauL => tauL(:,i,j)
          SCAT%ssaL => ssaL(:,i,j)
          SCAT%pmomL => pmomL(:,i,j,:,:)


          if ( verbose > 0 ) then
            print *,'brdf output1', SCAT%Surface%Base%VIO%VBRDF_Sup_Out%BS_DBOUNCE_BRDFUNC(1,1,1,1) 
            print *,'brdf output2', SCAT%Surface%Base%VIO%VBRDF_Sup_Out%BS_DBOUNCE_BRDFUNC(2,1,1,1) 
            print *,'brdf output3', SCAT%Surface%Base%VIO%VBRDF_Sup_Out%BS_DBOUNCE_BRDFUNC(3,1,1,1) 
          end if

          BRDF(j,i) = SCAT%Surface%Base%VIO%VBRDF_Sup_Out%BS_DBOUNCE_BRDFUNC(1,1,1,1) 
          BRDF_Q(j,i) = SCAT%Surface%Base%VIO%VBRDF_Sup_Out%BS_DBOUNCE_BRDFUNC(2,1,1,1) 
          BRDF_U(j,i) = SCAT%Surface%Base%VIO%VBRDF_Sup_Out%BS_DBOUNCE_BRDFUNC(3,1,1,1)     

         
          call VLIDORT_Run_Vector_Cloud (SCAT, output, ier)


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

end subroutine VLIDORT_Vector_GissCx_Cloud

subroutine VLIDORT_Vector_OCIGissCX_Cloud (km, nch, nobs,channels, nstreams, plane_parallel, nMom, &
                   nPol, ROT, depol, alpha, tau, ssa, pmom, tauI, ssaI, pmomI, tauL, ssaL, pmomL, &
                   pe, he, te, U10m, V10m, &
                   mr, solar_zenith, relat_azymuth, sensor_zenith, flux_factor, &
                   MISSING,verbose,radiance_VL_SURF,reflectance_VL_SURF, Q, U, BRDF, BRDF_Q, BRDF_U,rc)
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

  logical,          intent(in)  :: plane_parallel ! do plane parallel flag

  integer, target,  intent(in)  :: nMom  ! number of moments 
  integer, target,  intent(in)  :: nPol  ! number of components                               
  integer,          intent(in)  :: nstreams  ! number of half space streams
                  
  real*8, target,   intent(in)  :: channels(nch)    ! wavelengths [nm]

!                                                   ! --- Rayleigh Parameters ---
  real*8, target,   intent(in)  :: ROT(km,nobs,nch) ! rayleigh optical thickness
  real*8, target,   intent(in)  :: depol(nch)       ! rayleigh depolarization ratio used in phase matrix

!                                                   ! --- Trace Gas Absorption ---
  real*8, target,   intent(in)  :: alpha(km,nobs,nch) ! trace gas absoprtion optical thickness


!                                                   ! --- Mie Parameters ---
  real*8, target,   intent(in)  :: tau(km,nch,nobs) ! aerosol optical depth
  real*8, target,   intent(in)  :: ssa(km,nch,nobs) ! single scattering albedo
  real*8, target,   intent(in)  :: pmom(km,nch,nobs,nMom,nPol) !components of the scat phase matrix

  real*8, target,   intent(in)  :: tauI(km,nch,nobs) ! ice cloud optical depth
  real*8, target,   intent(in)  :: ssaI(km,nch,nobs) ! ice cloud single scattering albedo    
  real*8, target,   intent(in)  :: pmomI(km,nch,nobs,nMom,nPol) !components of the scat phase matrix


  real*8, target,   intent(in)  :: tauL(km,nch,nobs) ! liquid cloud optical depth
  real*8, target,   intent(in)  :: ssaL(km,nch,nobs) ! liquid cloud single scattering albedo  
  real*8, target,   intent(in)  :: pmomL(km,nch,nobs,nMom,nPol) !components of the scat phase matrix

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

  real*8,           intent(in)  :: flux_factor(nch,nobs) ! solar flux (F0)
  integer,          intent(in)  :: verbose

! !OUTPUT PARAMETERS:

  real*8,           intent(out) :: radiance_VL_SURF(nobs,nch)       ! TOA normalized radiance from VLIDORT
  integer,          intent(out) :: rc                          ! return code
  real*8,           intent(out) :: reflectance_VL_SURF(nobs, nch)   ! TOA reflectance from VLIDORT
  real*8,           intent(out) :: BRDF(nobs, nch)  
  real*8,           intent(out) :: BRDF_Q(nobs, nch)  
  real*8,           intent(out) :: BRDF_U(nobs, nch)      
  real*8,           intent(out) :: Q(nobs, nch)   ! Stokes parameter Q
  real*8,           intent(out) :: U(nobs, nch)   ! Stokes parameter U
!                               ---
  
  integer             :: i,j,n,p,ier 
  
  type(VLIDORT_scat) :: SCAT
  type(VLIDORT_output_vector)  :: output  

  rc = 0
  ier = 0
 
  ! set streams here
  SCAT%Surface%Base%NSTREAMS = nstreams
  SCAT%Surface%Base%DO_PLANE_PARALLEL = plane_parallel
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
          ! set solar flux
          SCAT%Surface%Base%VIO%VLIDORT_FixIn%SunRays%TS_FLUX_FACTOR = flux_factor(i,j)

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
            print*,U10m(j),V10m(j),mr(i),solar_zenith (j),&
                                    sensor_zenith(j),relat_azymuth(j)
          end if
          scalar = .false.
          call VLIDORT_OCIGissCoxMunk(SCAT%Surface,U10m(j),V10m(j),mr(i),solar_zenith (j),&
                                    sensor_zenith(j),relat_azymuth(j),scalar,rc)

          if ( verbose > 0 ) then
            print*, 'FINISHED COX MUNK'
          end if


          if ( rc /= 0 ) return

          SCAT%wavelength = channels(i)
          SCAT%rot => ROT(:,j,i)
          SCAT%depol_ratio => depol(i)
          SCAT%alpha => alpha(:,j,i)
          SCAT%tau => tau(:,i,j)
          SCAT%ssa => ssa(:,i,j)
          SCAT%pmom => pmom(:,i,j,:,:)
          SCAT%tauI => tauI(:,i,j)
          SCAT%ssaI => ssaI(:,i,j)
          SCAT%pmomI => pmomI(:,i,j,:,:)
          SCAT%tauL => tauL(:,i,j)
          SCAT%ssaL => ssaL(:,i,j)
          SCAT%pmomL => pmomL(:,i,j,:,:)


          if ( verbose > 0 ) then
            print *,'brdf output1', SCAT%Surface%Base%VIO%VBRDF_Sup_Out%BS_DBOUNCE_BRDFUNC(1,1,1,1) 
            print *,'brdf output2', SCAT%Surface%Base%VIO%VBRDF_Sup_Out%BS_DBOUNCE_BRDFUNC(2,1,1,1) 
            print *,'brdf output3', SCAT%Surface%Base%VIO%VBRDF_Sup_Out%BS_DBOUNCE_BRDFUNC(3,1,1,1) 
          end if

          BRDF(j,i) = SCAT%Surface%Base%VIO%VBRDF_Sup_Out%BS_DBOUNCE_BRDFUNC(1,1,1,1) 
          BRDF_Q(j,i) = SCAT%Surface%Base%VIO%VBRDF_Sup_Out%BS_DBOUNCE_BRDFUNC(2,1,1,1) 
          BRDF_U(j,i) = SCAT%Surface%Base%VIO%VBRDF_Sup_Out%BS_DBOUNCE_BRDFUNC(3,1,1,1)     

         
          call VLIDORT_Run_Vector_Cloud (SCAT, output, ier)


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

end subroutine VLIDORT_Vector_OCIGissCx_Cloud

!..........................................................................

subroutine VLIDORT_Scalar_GissCX_NOBM_Cloud (km, nch, nobs,channels, nstreams, plane_parallel, nMom, &
                   nPol, ROT, depol, tau, ssa, g, pmom, tauI, ssaI, gI, pmomI, tauL, ssaL, gL, pmomL, &
                   pe, he, te, U10m, V10m, &
                   mr, sleave, sleave_iso, sleave_adjust, &
                   solar_zenith, relat_azymuth, sensor_zenith, &
                   MISSING,verbose,radiance_VL_SURF,reflectance_Vl_SURF, &
                   BRDF,rc, ADJUSTED_SLEAVE)
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

  logical,          intent(in)  :: plane_parallel ! do plane parallel flag

  integer, target,  intent(in)  :: nMom             ! number of phase function moments    
  integer, target,  intent(in)  :: nPol  ! number of scattering matrix components                               
  integer,          intent(in)  :: nstreams  ! number of half space streams
                                                        
  real*8, target,   intent(in)  :: channels(nch)    ! wavelengths [nm]

!                                                   ! --- Rayleigh Parameters ---
  real*8, target,   intent(in)  :: ROT(km,nobs,nch) ! rayleigh optical thickness
  real*8, target,   intent(in)  :: depol(nch)       ! rayleigh depolarization ratio used in phase matrix


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

  real*8, target,   intent(in)  :: U10m(nobs)   ! Wind speed components [m/s]
  real*8, target,   intent(in)  :: V10m(nobs)    
  real*8, target,   intent(in)  :: mr(nch)       ! refractive index
  real*8, target,   intent(in)  :: sleave(nch,nobs)       ! sun normalized water leaving radiance from NOBM
                                                          ! described and used in Gregg & Rousseaux 2017
  logical,          intent(in)  :: sleave_iso  ! do isotropic water leaving radiance 
  logical, target,  intent(in)  :: sleave_adjust  ! do water leaving adjustment for transmittance.
  
  real*8, target,   intent(in)  :: MISSING          ! MISSING VALUE                                      
  real*8, target,   intent(in)  :: solar_zenith(nobs)  
  real*8, target,   intent(in)  :: relat_azymuth(nobs) 
  real*8, target,   intent(in)  :: sensor_zenith(nobs) 
  
  integer,          intent(in)  :: verbose

! !OUTPUT PARAMETERS:

  real*8,           intent(out) :: radiance_VL_SURF(nobs,nch)       ! TOA normalized radiance from VLIDORT
  integer,          intent(out) :: rc                          ! return code
  real*8,           intent(out) :: reflectance_VL_SURF(nobs, nch)   ! TOA reflectance from VLIDORT
  real*8,           intent(out) :: BRDF(nobs, nch)  
  real*8, optional, intent(out) :: ADJUSTED_SLEAVE(nobs,nch)
!                         ---  
  integer             :: i,j,n,p,ier

  type(VLIDORT_scat) :: SCAT
  type(VLIDORT_output_scalar)  :: output  

  rc = 0
  ier = 0

  ! set streams here
  SCAT%Surface%Base%NSTREAMS = nstreams
  SCAT%Surface%Base%DO_PLANE_PARALLEL = plane_parallel
  call VLIDORT_Init( SCAT%Surface%Base, km, rc)
  if ( rc /= 0 ) return

  SCAT%nMom    = nMom
  SCAT%NSTOKES = 1    


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
            print*,U10m(j),V10m(j),mr(i),solar_zenith (j),&
                                    sensor_zenith(j),relat_azymuth(j),scalar

          end if
          call VLIDORT_GissCoxMunk(SCAT%Surface,U10m(j),V10m(j),mr(i),solar_zenith (j),&
                                    sensor_zenith(j),relat_azymuth(j),scalar,rc)
          if ( rc /= 0 ) return
          ! must be called second
          if (sleave_iso) then
             call VLIDORT_NOBM_ISO(SCAT%Surface,sleave(i,j),solar_zenith(j),&
                                    sensor_zenith(j),relat_azymuth(j),&
                                    scalar,rc,sleave_adjust)
          else
             call VLIDORT_NOBM(SCAT%Surface,sleave(i,j),solar_zenith(j),&
                                    sensor_zenith(j),relat_azymuth(j),&
                                    scalar,rc,sleave_adjust)
          end if

          if ( rc /= 0 ) return

          SCAT%wavelength = channels(i)
          SCAT%rot => ROT(:,j,i)
          SCAT%depol_ratio => depol(i)
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

          if ( verbose > 0 ) then
            print *,'brdf output1', SCAT%Surface%Base%VIO%VBRDF_Sup_Out%BS_DBOUNCE_BRDFUNC(1,1,1,1) 
          end if

          BRDF(j,i) = SCAT%Surface%Base%VIO%VBRDF_Sup_Out%BS_DBOUNCE_BRDFUNC(1,1,1,1) 

          if ( verbose > 0 ) then
            print *, 'VLIDORT inputs'
            print*,channels(i),tau(1,i,j),ssa(1,i,j),g(1,i,j),pmom(1,i,j,1,1)
          end if

          call VLIDORT_Run_Scalar_Cloud (SCAT, output, ier)

          radiance_VL_SURF(j,i)    = output%radiance
          reflectance_VL_SURF(j,i) = output%reflectance
          if (sleave_adjust) then
              ADJUSTED_SLEAVE(j,i) = output%ADJUSTED_SLEAVE
          else
              if (present(ADJUSTED_SLEAVE)) then
                ADJUSTED_SLEAVE(j,i) = 0.0
              end if
          end if
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

end subroutine VLIDORT_Scalar_GissCX_NOBM_Cloud

!..........................................................................

subroutine VLIDORT_Scalar_OCIGissCX_NOBM_Cloud (km, nch, nobs,channels, nstreams, plane_parallel, nMom, &
                   nPol, ROT, depol, alpha, tau, ssa, g, pmom, tauI, ssaI, gI, pmomI, tauL, ssaL, gL, pmomL, &
                   pe, he, te, U10m, V10m, &
                   mr, sleave, sleave_iso, sleave_adjust, &
                   solar_zenith, relat_azymuth, sensor_zenith, flux_factor, &
                   MISSING,verbose,radiance_VL_SURF,reflectance_Vl_SURF, BRDF,rc,&
                   ADJUSTED_SLEAVE)
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

  logical,          intent(in)  :: plane_parallel ! do plane parallel flag

  integer, target,  intent(in)  :: nMom             ! number of phase function moments    
  integer, target,  intent(in)  :: nPol  ! number of scattering matrix components                               
  integer,          intent(in)  :: nstreams  ! number of half space streams
                                                        
  real*8, target,   intent(in)  :: channels(nch)    ! wavelengths [nm]

!                                                   ! --- Rayleigh Parameters ---
  real*8, target,   intent(in)  :: ROT(km,nobs,nch) ! rayleigh optical thickness
  real*8, target,   intent(in)  :: depol(nch)       ! rayleigh depolarization ratio used in phase matrix

!                                                   ! --- Trace Gas Absorption ---
  real*8, target,   intent(in)  :: alpha(km,nobs,nch) ! trace gas absoprtion optical thickness


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

  real*8, target,   intent(in)  :: U10m(nobs)   ! Wind speed components [m/s]
  real*8, target,   intent(in)  :: V10m(nobs)    
  real*8, target,   intent(in)  :: mr(nch)       ! refractive index
  real*8, target,   intent(in)  :: sleave(nch,nobs)       ! sun normalized water leaving radiance from NOBM
  logical,          intent(in)  :: sleave_iso  ! do isotropic water leaving radiance
  logical, target,  intent(in)  :: sleave_adjust  ! do water leaving adjustment for transmittance.
                                                          ! described and used in Gregg & Rousseaux 2017 
  
  real*8, target,   intent(in)  :: MISSING          ! MISSING VALUE                                      
  real*8, target,   intent(in)  :: solar_zenith(nobs)  
  real*8, target,   intent(in)  :: relat_azymuth(nobs) 
  real*8, target,   intent(in)  :: sensor_zenith(nobs) 

  real*8,           intent(in)  :: flux_factor(nch,nobs) ! solar flux (F0)  
  integer,          intent(in)  :: verbose

! !OUTPUT PARAMETERS:

  real*8,           intent(out) :: radiance_VL_SURF(nobs,nch)       ! TOA normalized radiance from VLIDORT
  integer,          intent(out) :: rc                          ! return code
  real*8,           intent(out) :: reflectance_VL_SURF(nobs, nch)   ! TOA reflectance from VLIDORT
  real*8,           intent(out) :: BRDF(nobs, nch)
  real*8, optional, intent(out) :: ADJUSTED_SLEAVE(nobs,nch)  
!                         ---  
  integer             :: i,j,n,p,ier

  type(VLIDORT_scat) :: SCAT
  type(VLIDORT_output_scalar)  :: output  

  rc = 0
  ier = 0

  ! set streams here
  SCAT%Surface%Base%NSTREAMS = nstreams
  SCAT%Surface%Base%DO_PLANE_PARALLEL = plane_parallel
  call VLIDORT_Init( SCAT%Surface%Base, km, rc)
  if ( rc /= 0 ) return

  SCAT%nMom    = nMom
  SCAT%NSTOKES = 1    


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
          ! set solar flux
          SCAT%Surface%Base%VIO%VLIDORT_FixIn%SunRays%TS_FLUX_FACTOR = flux_factor(i,j)
       
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
            print*,U10m(j),V10m(j),mr(i),solar_zenith (j),&
                                    sensor_zenith(j),relat_azymuth(j),scalar

          end if
          call VLIDORT_OCIGissCoxMunk(SCAT%Surface,U10m(j),V10m(j),mr(i),solar_zenith (j),&
                                    sensor_zenith(j),relat_azymuth(j),scalar,rc)
          if ( rc /= 0 ) return
          ! must be called second
          if (sleave_iso) then
             call VLIDORT_NOBM_ISO(SCAT%Surface,sleave(i,j),solar_zenith(j),&
                                    sensor_zenith(j),relat_azymuth(j),&
                                    scalar,rc,sleave_adjust)
          else
             call VLIDORT_NOBM(SCAT%Surface,sleave(i,j),solar_zenith(j),&
                                    sensor_zenith(j),relat_azymuth(j),&
                                    scalar,rc,sleave_adjust)
          end if

          if ( rc /= 0 ) return

          SCAT%wavelength = channels(i)
          SCAT%rot => ROT(:,j,i)
          SCAT%depol_ratio => depol(i)
          SCAT%alpha => alpha(:,j,i)
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

          if ( verbose > 0 ) then
            print *,'brdf output1', SCAT%Surface%Base%VIO%VBRDF_Sup_Out%BS_DBOUNCE_BRDFUNC(1,1,1,1) 
          end if

          BRDF(j,i) = SCAT%Surface%Base%VIO%VBRDF_Sup_Out%BS_DBOUNCE_BRDFUNC(1,1,1,1) 

          if ( verbose > 0 ) then
            print *, 'VLIDORT inputs'
            print*,channels(i),tau(1,i,j),ssa(1,i,j),g(1,i,j),pmom(1,i,j,1,1)
          end if

          call VLIDORT_Run_Scalar_Cloud (SCAT, output, ier)

          radiance_VL_SURF(j,i)    = output%radiance
          reflectance_VL_SURF(j,i) = output%reflectance

          if (sleave_adjust) then
              ADJUSTED_SLEAVE(j,i) = output%ADJUSTED_SLEAVE
          else
              if (present(ADJUSTED_SLEAVE)) then
                ADJUSTED_SLEAVE(j,i) = 0.0
              end if
          end if
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

end subroutine VLIDORT_Scalar_OCIGissCX_NOBM_Cloud


subroutine VLIDORT_Vector_GissCX_NOBM_Cloud (km, nch, nobs,channels, nstreams, plane_parallel, nMom, &
                   nPol, ROT, depol, tau, ssa, pmom, tauI, ssaI, pmomI, tauL, ssaL, pmomL, &
                   pe, he, te, U10m, V10m, &
                   mr, sleave, sleave_iso, sleave_adjust, &
                   solar_zenith, relat_azymuth, sensor_zenith, &
                   MISSING,verbose,radiance_VL_SURF,reflectance_VL_SURF, Q, U, BRDF, BRDF_Q, BRDF_U,rc,ADJUSTED_SLEAVE)
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

  logical,          intent(in)  :: plane_parallel ! do plane parallel flag

  integer, target,  intent(in)  :: nMom  ! number of moments 
  integer, target,  intent(in)  :: nPol  ! number of components                               
  integer,          intent(in)  :: nstreams  ! number of half space streams
                  
  real*8, target,   intent(in)  :: channels(nch)    ! wavelengths [nm]

!                                                   ! --- Rayleigh Parameters ---
  real*8, target,   intent(in)  :: ROT(km,nobs,nch) ! rayleigh optical thickness
  real*8, target,   intent(in)  :: depol(nch)       ! rayleigh depolarization ratio used in phase matrix

!                                                   ! --- Mie Parameters ---
  real*8, target,   intent(in)  :: tau(km,nch,nobs) ! aerosol optical depth
  real*8, target,   intent(in)  :: ssa(km,nch,nobs) ! single scattering albedo
  real*8, target,   intent(in)  :: pmom(km,nch,nobs,nMom,nPol) !components of the scat phase matrix

  real*8, target,   intent(in)  :: tauI(km,nch,nobs) ! ice cloud optical depth
  real*8, target,   intent(in)  :: ssaI(km,nch,nobs) ! ice cloud single scattering albedo    
  real*8, target,   intent(in)  :: pmomI(km,nch,nobs,nMom,nPol) !components of the scat phase matrix


  real*8, target,   intent(in)  :: tauL(km,nch,nobs) ! liquid cloud optical depth
  real*8, target,   intent(in)  :: ssaL(km,nch,nobs) ! liquid cloud single scattering albedo  
  real*8, target,   intent(in)  :: pmomL(km,nch,nobs,nMom,nPol) !components of the scat phase matrix

  real*8, target,   intent(in)  :: MISSING          ! MISSING VALUE
  real*8, target,   intent(in)  :: pe(km+1,nobs)    ! pressure at layer edges [Pa]
  real*8, target,   intent(in)  :: he(km+1,nobs)    ! height above sea-level  [m]
  real*8, target,   intent(in)  :: te(km+1,nobs)    ! temperature at layer edges [K]

  real*8, target,   intent(in)  :: U10m(nobs)   ! Wind speed components [m/s]
  real*8, target,   intent(in)  :: V10m(nobs)    
  real*8, target,   intent(in)  :: mr(nch)       ! refractive index
  real*8, target,   intent(in)  :: sleave(nch,nobs)       ! sun normalized water leaving radiance from NOBM
                                                          ! described and used in Gregg & Rousseaux 2017 
  logical,          intent(in)  :: sleave_iso  ! do isotropic water leaving radiance
  logical, target,  intent(in)  :: sleave_adjust  ! do water leaving adjustment for transmittance.
                       
  real*8, target,   intent(in)  :: solar_zenith(nobs)  
  real*8, target,   intent(in)  :: relat_azymuth(nobs) 
  real*8, target,   intent(in)  :: sensor_zenith(nobs) 

  integer,          intent(in)  :: verbose

! !OUTPUT PARAMETERS:

  real*8,           intent(out) :: radiance_VL_SURF(nobs,nch)       ! TOA normalized radiance from VLIDORT
  integer,          intent(out) :: rc                          ! return code
  real*8,           intent(out) :: reflectance_VL_SURF(nobs, nch)   ! TOA reflectance from VLIDORT
  real*8,           intent(out) :: BRDF(nobs, nch)  
  real*8,           intent(out) :: BRDF_Q(nobs, nch)  
  real*8,           intent(out) :: BRDF_U(nobs, nch)      
  real*8,           intent(out) :: Q(nobs, nch)   ! Stokes parameter Q
  real*8,           intent(out) :: U(nobs, nch)   ! Stokes parameter U

  real*8, optional, intent(out) :: ADJUSTED_SLEAVE(nobs,nch)  ! do water leaving adjustment for transmittance.
!                               ---
  
  integer             :: i,j,n,p,ier 
  
  type(VLIDORT_scat) :: SCAT
  type(VLIDORT_output_vector)  :: output  

  rc = 0
  ier = 0
 
  ! set streams here
  SCAT%Surface%Base%NSTREAMS = nstreams
  SCAT%Surface%Base%DO_PLANE_PARALLEL = plane_parallel
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
            print*,U10m(j),V10m(j),mr(i),solar_zenith (j),&
                                    sensor_zenith(j),relat_azymuth(j)
          end if
          scalar = .false.
          call VLIDORT_GissCoxMunk(SCAT%Surface,U10m(j),V10m(j),mr(i),solar_zenith (j),&
                                    sensor_zenith(j),relat_azymuth(j),scalar,rc)

          if ( rc /= 0 ) return
          ! must be called second
          if (sleave_iso) then
             call VLIDORT_NOBM_ISO(SCAT%Surface,sleave(i,j),solar_zenith(j),&
                                sensor_zenith(j),relat_azymuth(j),&
                                scalar,rc,sleave_adjust)
          else
             call VLIDORT_NOBM(SCAT%Surface,sleave(i,j),solar_zenith(j),&
                                sensor_zenith(j),relat_azymuth(j),&
                                scalar,rc,sleave_adjust)
          end if

          if ( verbose > 0 ) then
            print*, 'FINISHED COX MUNK'
          end if


          if ( rc /= 0 ) return

          SCAT%wavelength = channels(i)
          SCAT%rot => ROT(:,j,i)
          SCAT%depol_ratio => depol(i)
          SCAT%tau => tau(:,i,j)
          SCAT%ssa => ssa(:,i,j)
          SCAT%pmom => pmom(:,i,j,:,:)
          SCAT%tauI => tauI(:,i,j)
          SCAT%ssaI => ssaI(:,i,j)
          SCAT%pmomI => pmomI(:,i,j,:,:)
          SCAT%tauL => tauL(:,i,j)
          SCAT%ssaL => ssaL(:,i,j)
          SCAT%pmomL => pmomL(:,i,j,:,:)


          if ( verbose > 0 ) then
            print *,'brdf output1', SCAT%Surface%Base%VIO%VBRDF_Sup_Out%BS_DBOUNCE_BRDFUNC(1,1,1,1) 
            print *,'brdf output2', SCAT%Surface%Base%VIO%VBRDF_Sup_Out%BS_DBOUNCE_BRDFUNC(2,1,1,1) 
            print *,'brdf output3', SCAT%Surface%Base%VIO%VBRDF_Sup_Out%BS_DBOUNCE_BRDFUNC(3,1,1,1) 
          end if

          BRDF(j,i) = SCAT%Surface%Base%VIO%VBRDF_Sup_Out%BS_DBOUNCE_BRDFUNC(1,1,1,1) 
          BRDF_Q(j,i) = SCAT%Surface%Base%VIO%VBRDF_Sup_Out%BS_DBOUNCE_BRDFUNC(2,1,1,1) 
          BRDF_U(j,i) = SCAT%Surface%Base%VIO%VBRDF_Sup_Out%BS_DBOUNCE_BRDFUNC(3,1,1,1)     

         
          call VLIDORT_Run_Vector_Cloud (SCAT, output, ier)


          radiance_VL_SURF(j,i)    = output%radiance
          reflectance_VL_SURF(j,i) = output%reflectance
          Q(j,i)                   = output%Q
          U(j,i)                   = output%U                
          if (sleave_adjust) then
              ADJUSTED_SLEAVE(j,i) = output%ADJUSTED_SLEAVE
          else
              if (present(ADJUSTED_SLEAVE)) then
                ADJUSTED_SLEAVE(j,i) = 0.0
              end if
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

end subroutine VLIDORT_Vector_GissCx_NOBM_Cloud

subroutine VLIDORT_Vector_OCIGissCx_NOBM_Cloud (km, nch, nobs,channels, nstreams, plane_parallel, nMom, &
                   nPol, ROT, depol, alpha, tau, ssa, pmom, tauI, ssaI, pmomI, tauL, ssaL, pmomL, &
                   pe, he, te, U10m, V10m, &
                   mr, sleave, sleave_iso, sleave_adjust, &
                   solar_zenith, relat_azymuth, sensor_zenith, flux_factor, &
                   MISSING,verbose,radiance_VL_SURF,reflectance_VL_SURF, &
                   Q, U, BRDF, BRDF_Q, BRDF_U,rc,ADJUSTED_SLEAVE)
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

  logical,          intent(in)  :: plane_parallel ! do plane parallel flag

  integer, target,  intent(in)  :: nMom  ! number of moments 
  integer, target,  intent(in)  :: nPol  ! number of components                               
  integer,          intent(in)  :: nstreams  ! number of half space streams
                  
  real*8, target,   intent(in)  :: channels(nch)    ! wavelengths [nm]

!                                                   ! --- Rayleigh Parameters ---
  real*8, target,   intent(in)  :: ROT(km,nobs,nch) ! rayleigh optical thickness
  real*8, target,   intent(in)  :: depol(nch)       ! rayleigh depolarization ratio used in phase matrix

!                                                   ! --- Trace Gas Absorption ---
  real*8, target,   intent(in)  :: alpha(km,nobs,nch) ! trace gas absoprtion optical thickness

!                                                   ! --- Mie Parameters ---
  real*8, target,   intent(in)  :: tau(km,nch,nobs) ! aerosol optical depth
  real*8, target,   intent(in)  :: ssa(km,nch,nobs) ! single scattering albedo
  real*8, target,   intent(in)  :: pmom(km,nch,nobs,nMom,nPol) !components of the scat phase matrix

  real*8, target,   intent(in)  :: tauI(km,nch,nobs) ! ice cloud optical depth
  real*8, target,   intent(in)  :: ssaI(km,nch,nobs) ! ice cloud single scattering albedo    
  real*8, target,   intent(in)  :: pmomI(km,nch,nobs,nMom,nPol) !components of the scat phase matrix


  real*8, target,   intent(in)  :: tauL(km,nch,nobs) ! liquid cloud optical depth
  real*8, target,   intent(in)  :: ssaL(km,nch,nobs) ! liquid cloud single scattering albedo  
  real*8, target,   intent(in)  :: pmomL(km,nch,nobs,nMom,nPol) !components of the scat phase matrix

  real*8, target,   intent(in)  :: MISSING          ! MISSING VALUE
  real*8, target,   intent(in)  :: pe(km+1,nobs)    ! pressure at layer edges [Pa]
  real*8, target,   intent(in)  :: he(km+1,nobs)    ! height above sea-level  [m]
  real*8, target,   intent(in)  :: te(km+1,nobs)    ! temperature at layer edges [K]

  real*8, target,   intent(in)  :: U10m(nobs)   ! Wind speed components [m/s]
  real*8, target,   intent(in)  :: V10m(nobs)    
  real*8, target,   intent(in)  :: mr(nch)       ! refractive index
  real*8, target,   intent(in)  :: sleave(nch,nobs)       ! sun normalized water leaving radiance from NOBM
                                                          ! described and used in Gregg & Rousseaux 2017 
  logical,          intent(in)  :: sleave_iso  ! do isotropic water leaving radiance
  logical, target,  intent(in)  :: sleave_adjust  ! do water leaving adjustment for transmittance.
                       
  real*8, target,   intent(in)  :: solar_zenith(nobs)  
  real*8, target,   intent(in)  :: relat_azymuth(nobs) 
  real*8, target,   intent(in)  :: sensor_zenith(nobs) 

  real*8,           intent(in)  :: flux_factor(nch,nobs) ! solar flux (F0)
  integer,          intent(in)  :: verbose

! !OUTPUT PARAMETERS:

  real*8,           intent(out) :: radiance_VL_SURF(nobs,nch)       ! TOA normalized radiance from VLIDORT
  integer,          intent(out) :: rc                          ! return code
  real*8,           intent(out) :: reflectance_VL_SURF(nobs, nch)   ! TOA reflectance from VLIDORT
  real*8,           intent(out) :: BRDF(nobs, nch)  
  real*8,           intent(out) :: BRDF_Q(nobs, nch)  
  real*8,           intent(out) :: BRDF_U(nobs, nch)      
  real*8,           intent(out) :: Q(nobs, nch)   ! Stokes parameter Q
  real*8,           intent(out) :: U(nobs, nch)   ! Stokes parameter U
  real*8, optional, intent(out) :: ADJUSTED_SLEAVE(nobs,nch)
!                               ---
  
  integer             :: i,j,n,p,ier 
  
  type(VLIDORT_scat) :: SCAT
  type(VLIDORT_output_vector)  :: output  

  rc = 0
  ier = 0
 
  ! set streams here
  SCAT%Surface%Base%NSTREAMS = nstreams
  SCAT%Surface%Base%DO_PLANE_PARALLEL = plane_parallel
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
          ! set solar flux
          SCAT%Surface%Base%VIO%VLIDORT_FixIn%SunRays%TS_FLUX_FACTOR = flux_factor(i,j)
     
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
            print*,U10m(j),V10m(j),mr(i),solar_zenith (j),&
                                    sensor_zenith(j),relat_azymuth(j)
          end if
          scalar = .false.
          call VLIDORT_OCIGissCoxMunk(SCAT%Surface,U10m(j),V10m(j),mr(i),solar_zenith (j),&
                                    sensor_zenith(j),relat_azymuth(j),scalar,rc)

          if ( rc /= 0 ) return
          ! must be called second
          if (sleave_iso) then
             call VLIDORT_NOBM_ISO(SCAT%Surface,sleave(i,j),solar_zenith(j),&
                                    sensor_zenith(j),relat_azymuth(j),&
                                    scalar,rc,sleave_adjust)
          else
             call VLIDORT_NOBM(SCAT%Surface,sleave(i,j),solar_zenith(j),&
                                    sensor_zenith(j),relat_azymuth(j),&
                                    scalar,rc,sleave_adjust)
          end if
          if ( verbose > 0 ) then
            print*, 'FINISHED COX MUNK'
          end if


          if ( rc /= 0 ) return

          SCAT%wavelength = channels(i)
          SCAT%rot => ROT(:,j,i)
          SCAT%depol_ratio => depol(i)
          SCAT%alpha => alpha(:,j,i)
          SCAT%tau => tau(:,i,j)
          SCAT%ssa => ssa(:,i,j)
          SCAT%pmom => pmom(:,i,j,:,:)
          SCAT%tauI => tauI(:,i,j)
          SCAT%ssaI => ssaI(:,i,j)
          SCAT%pmomI => pmomI(:,i,j,:,:)
          SCAT%tauL => tauL(:,i,j)
          SCAT%ssaL => ssaL(:,i,j)
          SCAT%pmomL => pmomL(:,i,j,:,:)


          if ( verbose > 0 ) then
            print *,'brdf output1', SCAT%Surface%Base%VIO%VBRDF_Sup_Out%BS_DBOUNCE_BRDFUNC(1,1,1,1) 
            print *,'brdf output2', SCAT%Surface%Base%VIO%VBRDF_Sup_Out%BS_DBOUNCE_BRDFUNC(2,1,1,1) 
            print *,'brdf output3', SCAT%Surface%Base%VIO%VBRDF_Sup_Out%BS_DBOUNCE_BRDFUNC(3,1,1,1) 
          end if

          BRDF(j,i) = SCAT%Surface%Base%VIO%VBRDF_Sup_Out%BS_DBOUNCE_BRDFUNC(1,1,1,1) 
          BRDF_Q(j,i) = SCAT%Surface%Base%VIO%VBRDF_Sup_Out%BS_DBOUNCE_BRDFUNC(2,1,1,1) 
          BRDF_U(j,i) = SCAT%Surface%Base%VIO%VBRDF_Sup_Out%BS_DBOUNCE_BRDFUNC(3,1,1,1)     

         
          call VLIDORT_Run_Vector_Cloud (SCAT, output, ier)


          radiance_VL_SURF(j,i)    = output%radiance
          reflectance_VL_SURF(j,i) = output%reflectance
          Q(j,i)                   = output%Q
          U(j,i)                   = output%U                

          if ( ier /= 0 ) then
              radiance_VL_SURF(j,i) = MISSING
              reflectance_VL_SURF(j,i) = MISSING                          
              cycle
          end if
          if (sleave_adjust) then
              ADJUSTED_SLEAVE(j,i) = output%ADJUSTED_SLEAVE
          else
              if (present(ADJUSTED_SLEAVE)) then
                ADJUSTED_SLEAVE(j,i) = 0.0
              end if
          end if

        end do ! end loop over channels
     
        if ( verbose > 0 ) then
           if ( mod(j-1,1000) == 0 ) then
              print *, '<> VLIDORT Vector: ', nint(j*100./nobs), '%'
           end if
        end if

  end do ! Loop over obs

end subroutine VLIDORT_Vector_OCIGissCx_NOBM_Cloud

end module VLIDORT_BRDF_CX
