module VLIDORT_BRDF_MODIS   
!
!  Simple f77 wrapper for the Python interface to VLIDORT for OMI aerosol
!  channels.
!
!.............................................................................
  implicit NONE

  PUBLIC Scalar_LandMODIS
  PUBLIC Vector_LandMODIS

  contains

  logical function IS_MISSING(x,MISSING)
    real*8, intent(in)     :: x
    real*8, intent(in)     :: MISSING

    IS_MISSING = abs(x/MISSING-1)<0.001
    return
  end function IS_MISSING

  subroutine Scalar_LandMODIS (km, nch, nobs,channels,        &
                     tau, ssa, g, pe, he, te, kernel_wt, param, &
                     solar_zenith, relat_azymuth, sensor_zenith, &
                     MISSING,verbose,radiance_VL_SURF,reflectance_VL_SURF, ROT, WSA, BSA, rc, &
                     albedo, reflectance_VL,radiance_VL )
  !
  ! Uses VLIDORT in scalar mode to compute OMI aerosol TOA radiances.
  !
    use VLIDORT_ScatMod

    implicit NONE

    integer, parameter            :: nkernel = 3
    integer, parameter            :: nparam  = 2
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

    real*8, target,   intent(in)  :: kernel_wt(nkernel,nch,nobs)   ! kernel weights (/fiso,fgeo,fvol/)
    real*8, target,   intent(in)  :: param(nparam,nch,nobs)        ! Li-Sparse parameters 
                                                                   ! param1 = crown relative height (h/b)
                                                                   ! param2 = shape parameter (b/r)
    
    real*8, target,   intent(in)  :: MISSING          ! MISSING VALUE                                      
    real*8, target,   intent(in)  :: solar_zenith(nobs)  
    real*8, target,   intent(in)  :: relat_azymuth(nobs) 
    real*8, target,   intent(in)  :: sensor_zenith(nobs) 

    real*8, optional, target,   intent(in)  :: albedo(nobs,nch)       ! surface albedo
    
    integer,          intent(in)            :: verbose

  ! !OUTPUT PARAMETERS:

    real*8,           intent(out) :: radiance_VL_SURF(nobs,nch)       ! TOA normalized radiance from VLIDORT
    real*8,           intent(out) :: reflectance_VL_SURF(nobs, nch)   ! TOA reflectance from VLIDORT using surface module
    integer,          intent(out) :: rc                               ! return code
    real*8,           intent(out) :: ROT(km,nobs,nch)                    ! rayleigh optical thickness  
    real*8,           intent(out) :: WSA(nobs,nch)                    ! White-sky albedo 
    real*8,           intent(out) :: BSA(nobs,nch)                    ! Black-sky albedo     
    real*8, optional, intent(out) :: reflectance_VL(nobs, nch)        ! TOA reflectance from VLIDORT using albedo
    real*8, optional, intent(out) :: radiance_VL(nobs,nch)            ! TOA normalized radiance from VLIDORT

    real*8                        :: Q(nobs, nch)                     ! Stokes parameter Q
    real*8                        :: U(nobs, nch)                     ! Stokes parameter U  
  !  real*8,           intent(out) :: BRDF(nobs, nch)  
  !                         ---  
    integer             :: i,j,n,p,ier

  
    type(VLIDORT_scat) :: SCAT

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
        ! Make sure kernel weights and parameters are defined
        do n = 1, nkernel
          if (IS_MISSING(kernel_wt(n,i,j),MISSING)) then
            radiance_VL_SURF(j,i) = MISSING
            reflectance_VL_SURF(j,i) = MISSING
            cycle
          end if
        end do
        do p = 1, nparam
          if (IS_MISSING(param(p,i,j),MISSING)) then
            radiance_VL_SURF(j,i) = MISSING
            reflectance_VL_SURF(j,i) = MISSING
            cycle
          end if
        end do

        if ( verbose > 0 ) then
          print*, 'DO MODIS BRDF'
        end if
        call VLIDORT_LANDMODIS(SCAT%Surface,solar_zenith(j),&
                               sensor_zenith(j),relat_azymuth(j),&
                               kernel_wt(1,i,j),kernel_wt(2,i,j),kernel_wt(3,i,j),&
                               reshape(param(:,i,j),(/nparam/)),&
                               .true.,rc)

        if ( rc /= 0 ) return

        SCAT%wavelength = channels(i)
        SCAT%tau => tau(:,i,j)
        SCAT%ssa => ssa(:,i,j)
        SCAT%g => g(:,i,j)
         
        call VLIDORT_Run (SCAT, radiance_VL_SURF(j,i), reflectance_VL_SURF(j,i), &
                          ROT(:,j,i), Q(j,i), U(j,i), .true., .true., ier)

        WSA(j,i) = SCAT%Surface%Base%VIO%VBRDF_Sup_Out%BS_WSA_CALCULATED
        BSA(j,i) = SCAT%Surface%Base%VIO%VBRDF_Sup_Out%BS_BSA_CALCULATED

        if ( verbose > 0 ) then
          print *, 'My radiance land modis',radiance_VL_SURF(j,i), reflectance_VL_SURF(j,i) 
        end if

        if ( ier /= 0 ) then
          radiance_VL(j,i) = MISSING
          reflectance_VL(j,i) = MISSING
          cycle
        end if

        if ( verbose > 0 ) then
          print *, 'now check for albedo'       
        end if

        ! Mare sure albedo is defined
        ! ---------------------------
        if (present(albedo)) then
          if ( verbose > 0 ) then
            print *, 'albedo present'
        end if 
        if ( IS_MISSING(albedo(j,i),MISSING) ) then
          radiance_VL(j,i) = MISSING
          reflectance_VL(j,i) = MISSING
          cycle
        end if
        else
          if ( verbose > 0 ) then
            print *, 'albedo not present'
          end if
          cycle
        end if
        if ( verbose > 0 ) then
          print*, 'DO SURFACE LAMB'
        end if
        call VLIDORT_SurfaceLamb(SCAT%Surface,albedo(j,i),solar_zenith (j),sensor_zenith(j),&
                               relat_azymuth(j),.true.)
         
        call VLIDORT_Run (SCAT, radiance_VL(j,i), reflectance_VL(j,i), &
                        ROT(:,j,i), Q(j,i), U(j,i), .true., .true., ier)
        if ( verbose > 0 ) then
          print *, 'radiance albedo',albedo(j,i),radiance_VL(j,i), reflectance_VL(j,i) 
        end if
        if ( ier /= 0 ) then
          radiance_VL(j,i) = MISSING
          reflectance_VL(j,i) = MISSING
          cycle
        end if
      end do ! end loop over channels
     
      if ( verbose > 0 ) then
        if ( mod(j-1,1000) == 0 ) then
          print *, '<> VLIDORT Scalar: ', nint(j*100./nobs), '%'
        end if
      end if

    end do ! Loop over obs

  end subroutine Scalar_LandMODIS

  !..........................................................................

  subroutine Vector_LandMODIS (km, nch, nobs, channels, nMom,  &
                     nPol,tau, ssa, g, pmom, pe, he, te, kernel_wt, param, &
                     solar_zenith, relat_azymuth, sensor_zenith, &
                     MISSING,verbose, radiance_VL_SURF,reflectance_VL_SURF, ROT, WSA, BSA, Q, U, rc, &
                     albedo, reflectance_VL,radiance_VL, Q_lamb, U_lamb)
  !
  ! Place holder.
  !
     use VLIDORT_ScatMod
   
     implicit NONE

    integer, parameter            :: nkernel = 3
    integer, parameter            :: nparam  = 2

  ! !INPUT PARAMETERS:

    integer,          intent(in)            :: km    ! number of levels on file
    integer,          intent(in)            :: nch   ! number of channels
    integer,          intent(in)            :: nobs  ! number of observations

    integer, target,  intent(in)            :: nMom  ! number of moments 
    integer, target,  intent(in)            :: nPol  ! number of components                               
                    
    real*8, target,   intent(in)            :: channels(nch)    ! wavelengths [nm]

  !                                                   ! --- Mie Parameters ---
    real*8, target,   intent(in)            :: tau(km,nch,nobs) ! aerosol optical depth
    real*8, target,   intent(in)            :: ssa(km,nch,nobs) ! single scattering albedo
    real*8, target,   intent(in)            :: g(km,nch,nobs)   ! asymmetry factor
    
    real*8, target,   intent(in)            :: pmom(km,nMom,nPol,nch,nobs) !components of the scat phase matrix

    real*8, target,   intent(in)            :: MISSING          ! MISSING VALUE
    real*8, target,   intent(in)            :: pe(km+1,nobs)    ! pressure at layer edges [Pa]
    real*8, target,   intent(in)            :: he(km+1,nobs)    ! height above sea-level  [m]
    real*8, target,   intent(in)            :: te(km+1,nobs)    ! temperature at layer edges [K]

    real*8, target,   intent(in)            :: kernel_wt(nkernel,nch,nobs)   ! kernel weights (/fiso,fgeo,fvol/)
    real*8, target,   intent(in)            :: param(nparam,nch,nobs)        ! Li-Sparse parameters 
                                                                   ! param1 = crown relative height (h/b)
                                                                   ! param2 = shape parameter (b/r)
                         
    real*8, target,   intent(in)            :: solar_zenith(nobs)  
    real*8, target,   intent(in)            :: relat_azymuth(nobs) 
    real*8, target,   intent(in)            :: sensor_zenith(nobs) 

    real*8, optional, target,   intent(in)  :: albedo(nobs,nch)       ! surface albedo

    integer,          intent(in)            :: verbose

  ! !OUTPUT PARAMETERS:
    real*8,           intent(out)           :: radiance_VL_SURF(nobs,nch)         ! TOA normalized radiance from VLIDORT
    real*8,           intent(out)           :: reflectance_VL_SURF(nobs, nch)   ! TOA reflectance from VLIDORT using surface module
    integer,          intent(out)           :: rc                            ! return code

    real*8,           intent(out)           :: ROT(km,nobs,nch)                    ! rayleigh optical thickness
    real*8,           intent(out)           :: WSA(nobs,nch)                    ! white sky albedo 
    real*8,           intent(out)           :: BSA(nobs,nch)                    ! black sky albedo 
    real*8,           intent(out)           :: Q(nobs, nch)                     ! Stokes parameter Q
    real*8,           intent(out)           :: U(nobs, nch)                     ! Stokes parameter U   

    real*8, optional, intent(out)           :: reflectance_VL(nobs, nch)     ! TOA reflectance from VLIDORT using albedo
    real*8, optional, intent(out)           :: radiance_VL(nobs,nch)         ! TOA normalized radiance from VLIDORT
    real*8, optional, intent(out)           :: Q_lamb(nobs, nch)   ! Stokes parameter Q
    real*8, optional, intent(out)           :: U_lamb(nobs, nch)   ! Stokes parameter U    
  !                               ---
    
    integer             :: i,j,n,p,ier

    
    type(VLIDORT_scat) :: SCAT

    rc = 0
    ier = 0
   
    call VLIDORT_Init( SCAT%Surface%Base, km, rc)
    if ( rc /= 0 ) return

    SCAT%nMom = nMom
    SCAT%nPol = nPol

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
          if (IS_MISSING(param(p,i,j),MISSING)) then
            radiance_VL_SURF(j,i) = MISSING
            reflectance_VL_SURF(j,i) = MISSING
            cycle
          end if
        end do

        if ( verbose > 0 ) then
          print*, 'DO MODIS BRDF'
        end if

        call VLIDORT_LANDMODIS(SCAT%Surface,solar_zenith(j),&
                               sensor_zenith(j),relat_azymuth(j),&
                               kernel_wt(1,i,j),kernel_wt(2,i,j),kernel_wt(3,i,j),&
                               reshape(param(:,i,j),(/nparam/)),&
                               .true.,rc)

        SCAT%wavelength = channels(i)        
        SCAT%tau => tau(:,i,j)
        SCAT%ssa => ssa(:,i,j)
        SCAT%g => g(:,i,j)
        SCAT%pmom => pmom(:,:,:,i,j)

        call VLIDORT_Run (SCAT, radiance_VL_SURF(j,i),reflectance_VL_SURF(j,i),&
                          ROT(:,j,i), Q(j,i), U(j,i), .false., .true., ier)

        WSA(j,i) = SCAT%Surface%Base%VIO%VBRDF_Sup_Out%BS_WSA_CALCULATED
        BSA(j,i) = SCAT%Surface%Base%VIO%VBRDF_Sup_Out%BS_BSA_CALCULATED
        
        if ( ier /= 0 ) then
          radiance_VL_SURF(j,i) = MISSING
          reflectance_VL_SURF(j,i) = MISSING               
          cycle
        end if

        ! Mare sure albedo is defined
        ! ---------------------------
        if (present(albedo)) then
          if ( verbose > 0 ) then
            print *, 'albedo present'
          end if 
          if ( IS_MISSING(albedo(j,i),MISSING) ) then
            radiance_VL(j,i) = MISSING
            reflectance_VL(j,i) = MISSING
            cycle
          end if
        else
          if ( verbose > 0 ) then
            print *, 'albedo not present'
          end if
          cycle
        end if
        
        if ( verbose > 0 ) then
          print*, 'DO SURFACE LAMB'
        end if
        call VLIDORT_SurfaceLamb(SCAT%Surface,albedo(j,i),solar_zenith (j),sensor_zenith(j),&
                                 relat_azymuth(j),.true.)
        if ( rc /= 0 ) return
           
        call VLIDORT_Run (SCAT, radiance_VL(j,i), reflectance_VL(j,i), &
                          ROT(:,j,i), Q_lamb(j,i), U_lamb(j,i), .false., .true., ier)
 
        if ( ier /= 0 ) then
          radiance_VL(j,i) = MISSING
          reflectance_VL(j,i) = MISSING
          cycle
        end if
      end do ! end loop over channels
       
      if ( verbose > 0 ) then
        if ( mod(j-1,1000) == 0 ) then
          print *, '<> VLIDORT Vector: ', nint(j*100./nobs), '%'
        end if
      end if

    end do ! Loop over obs

  end subroutine Vector_LandMODIS

end module VLIDORT_BRDF_MODIS