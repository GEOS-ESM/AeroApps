module VLIDORT_SURFACE   
!
!  F90 Driver to run VLIDROT surface supplement 
!
!.............................................................................
  implicit NONE

  PUBLIC VLIDORT_SURFACE_LandMODIS

  contains

  logical function IS_MISSING(x,MISSING)
    real*8, intent(in)     :: x
    real*8, intent(in)     :: MISSING

    IS_MISSING = abs(x/MISSING-1)<0.001
    return
  end function IS_MISSING

  subroutine VLIDORT_SURFACE_LandMODIS (km, nch, nobs, channels,                    &
                                kernel_wt, param,                           &
                                solar_zenith, relat_azymuth, sensor_zenith, &
                                MISSING, verbose, BR, WSA, BSA, rc )
  !
  ! Uses VLIDORT surface supplement to compute surface reflectance
  ! from MODIS RTLS kernel weights
  !
    use VLIDORT_ScatMod

    implicit NONE

    logical, parameter            :: scalar = .true.
    integer, parameter            :: nkernel = 3
    integer, parameter            :: nparam  = 2
  ! !INPUT PARAMETERS:
    integer,          intent(in)  :: km    ! number of vertical levels
    integer,          intent(in)  :: nch   ! number of channels
    integer,          intent(in)  :: nobs  ! number of observations
                                        
                    
    real*8, target,   intent(in)  :: channels(nch)    ! wavelengths [nm]


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

    integer,          intent(out) :: rc                               ! return code
    real*8,           intent(out) :: WSA(nobs,nch)                    ! White-sky albedo
    real*8,           intent(out) :: BSA(nobs,nch)                    ! Black-sky albedo        
    real*8,           intent(out) :: BR(nobs,nch)                     ! bi-directional reflectance 
    real*8, optional, intent(out) :: reflectance_VL(nobs, nch)        ! TOA reflectance from VLIDORT using albedo
  !                         ---  
    integer             :: i,j,n,p,ier

  
    type(VLIDORT_scat) :: SCAT

    rc = 0
    ier = 0

    call VLIDORT_Init( SCAT%Surface%Base, km, rc)
    if ( rc /= 0 ) return

    do j = 1,nobs

       ! Make sure angles are available
       ! -----------------------------------------
      if ( IS_MISSING(solar_zenith(j),MISSING)  .OR. & 
           IS_MISSING(sensor_zenith(j),MISSING) .OR. &
           IS_MISSING(relat_azymuth(j),MISSING)  )  then

        WSA(j,:) = MISSING
        BSA(j,:) = MISSING
        BR(j,:)  = MISSING
        cycle
      end if
           
      ! Loop over channels
      ! ------------------
      do i = 1, nch 
        ! Make sure kernel weights and parameters are defined
        do n = 1, nkernel
          if (IS_MISSING(kernel_wt(n,i,j),MISSING)) then
            WSA(j,:) = MISSING
            BSA(j,:) = MISSING
            BR(j,:)  = MISSING
            cycle
          end if
        end do
        do p = 1, nparam
          if (IS_MISSING(param(p,i,j),MISSING)) then
            WSA(j,:) = MISSING
            BSA(j,:) = MISSING
            BR(j,:)  = MISSING
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
                               scalar,rc)

        if ( rc /= 0 ) return       

        if ( ier /= 0 ) then
          WSA(j,:) = MISSING
          BSA(j,:) = MISSING
          BR(j,:)  = MISSING
          cycle
        else
          BR(j,i)  = SCAT%Surface%Base%VIO%VBRDF_Sup_Out%BS_DBOUNCE_BRDFUNC(1,1,1,1)
          WSA(j,i) = SCAT%Surface%Base%VIO%VBRDF_Sup_Out%BS_WSA_CALCULATED
          BSA(j,i) = SCAT%Surface%Base%VIO%VBRDF_Sup_Out%BS_BSA_CALCULATED
        end if

        if ( verbose > 0 ) then
          print *, 'My radiance land modis',BR(j,i), WSA(j,i), BSA(j,i) 
        end if

      end do ! end loop over channels
     
      if ( verbose > 0 ) then
        if ( mod(j-1,1000) == 0 ) then
          print *, '<> VLIDORT Scalar: ', nint(j*100./nobs), '%'
        end if
      end if

    end do ! Loop over obs

  end subroutine VLIDORT_SURFACE_LandMODIS

  !..........................................................................


end module VLIDORT_SURFACE