module LIDORT_BRDF_MODIS   
!
!  Simple f77 wrapper for the Python interface to LIDORT for OMI aerosol
!  channels.
!
!.............................................................................
  implicit NONE

  PUBLIC LIDORT_Scalar_LandMODIS
  PUBLIC LIDORT_Vector_LandMODIS

  contains

  logical function IS_MISSING(x,MISSING)
    real*8, intent(in)     :: x
    real*8, intent(in)     :: MISSING

    IS_MISSING = abs(x/MISSING-1)<0.001
    return
  end function IS_MISSING

  subroutine LIDORT_Scalar_LandMODIS (km, nch, nobs,channels,        &
                     tau, ssa, g, pe, he, te, kernel_wt, param, &
                     solar_zenith, relat_azymuth, sensor_zenith, &
                     MISSING,verbose,radiance_VL_SURF,reflectance_VL_SURF, ROT, BR, rc, &
                     albedo, reflectance_L,radiance_L )
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
                                        
                    
    real*8, target,   intent(in)  :: channels(nch)    ! wavelengths [nm]

  !                                                   ! --- Aerosol Optical Properties ---
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

    real*8, optional, target,   intent(in)  :: albedo(nobs,nch)       ! surface albedo optional
    
    integer,          intent(in)            :: verbose

  ! !OUTPUT PARAMETERS:

    real*8,           intent(out) :: radiance_L_SURF(nobs,nch)       ! TOA normalized radiance from LIDORT using surface module
    real*8,           intent(out) :: reflectance_L_SURF(nobs, nch)   ! TOA reflectance from LIDORT using surface module
    integer,          intent(out) :: rc                               ! return code
    real*8,           intent(out) :: ROT(km,nobs,nch)                 ! rayleigh optical thickness  
    real*8,           intent(out) :: BR(nobs,nch)                     ! bi-directional reflectance 
    real*8, optional, intent(out) :: reflectance_L(nobs, nch)        ! TOA reflectance from LIDORT using albedo
    real*8, optional, intent(out) :: radiance_L(nobs,nch)            ! TOA normalized radiance from LIDORT using albedo

  !  real*8,           intent(out) :: BRDF(nobs, nch)  
  !                         ---  
    integer             :: i,j,n,p,ier

  
    type(LIDORT_scat) :: SCAT

    rc = 0
    ier = 0

    call LIDORT_Init( SCAT%Surface%Base, km, rc)
    if ( rc /= 0 ) return

    do j = 1,nobs

       ! Make sure albedo and angles are available
       ! -----------------------------------------
      if ( IS_MISSING(solar_zenith(j),MISSING)  .OR. & 
           IS_MISSING(sensor_zenith(j),MISSING) .OR. &
           IS_MISSING(relat_azymuth(j),MISSING)  )  then

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
         
        call LIDORT_Run (SCAT, radiance_L_SURF(j,i), reflectance_L_SURF(j,i), &
                          ROT(:,j,i), Q(j,i), U(j,i), aerosol, ier)

        BR(j,i) = SCAT%Surface%Base%VIO%VBRDF_Sup_Out%BS_DBOUNCE_BRDFUNC(1,1,1,1)

        if ( verbose > 0 ) then
          print *, 'My radiance land modis',radiance_L_SURF(j,i), reflectance_L_SURF(j,i) 
        end if

        if ( ier /= 0 ) then
          radiance_L(j,i) = MISSING
          reflectance_L(j,i) = MISSING
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
          radiance_L(j,i) = MISSING
          reflectance_L(j,i) = MISSING
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
        call LIDORT_SurfaceLamb(SCAT%Surface,albedo(j,i),solar_zenith (j),sensor_zenith(j),&
                               relat_azymuth(j))
         
        call LIDORT_Run (SCAT, radiance_L(j,i), reflectance_L(j,i), &
                        ROT(:,j,i), Q(j,i), U(j,i), aerosol, ier)
        if ( verbose > 0 ) then
          print *, 'radiance albedo',albedo(j,i),radiance_L(j,i), reflectance_L(j,i) 
        end if
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

  end subroutine LIDORT_Scalar_LandMODIS

  !..........................................................................
end module LIDORT_BRDF_MODIS