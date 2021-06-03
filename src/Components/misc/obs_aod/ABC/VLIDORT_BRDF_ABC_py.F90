!
!  Simple f90 wrapper for the Python interface to VLIDORT Surface Reflectance
!  Calculators
!
!.............................................................................

subroutine CoxMunk (km, nch, nobs,channels, U10m,V10m, &
                   mr, solar_zenith, relat_azymuth, sensor_zenith, &
                   MISSING,verbose, BRDF,rc)
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

  real*8, target,   intent(in)  :: U10m(nobs)   ! Wind speed components [m/s]
  real*8, target,   intent(in)  :: V10m(nobs)    
  real*8, target,   intent(in)  :: mr(nch)       ! refractive index
  
  real*8, target,   intent(in)  :: MISSING          ! MISSING VALUE                                      
  real*8, target,   intent(in)  :: solar_zenith(nobs)  
  real*8, target,   intent(in)  :: relat_azymuth(nobs) 
  real*8, target,   intent(in)  :: sensor_zenith(nobs) 
  
  integer,          intent(in)  :: verbose

! !OUTPUT PARAMETERS:

  integer,          intent(out) :: rc                          ! return code
  real*8,           intent(out) :: BRDF(nobs, nch)  
!                         ---  
  integer             :: i,j, ier

  type(VLIDORT_scat) :: SCAT

#define IS_MISSING(x) (abs(x/MISSING-1)<0.001)
  rc = 0
  ier = 0

  call VLIDORT_Init( SCAT%Surface%Base, km, rc)
  if ( rc /= 0 ) return

  do j = 1,nobs

     ! Make sure albedo and angles are available
     ! -----------------------------------------
     if ( IS_MISSING(solar_zenith(j))  .OR. & 
          IS_MISSING(sensor_zenith(j)) .OR. &
          IS_MISSING(relat_azymuth(j))  )  then

        BRDF(j,:) = MISSING
        cycle

      end if
      
 
      ! Loop over channels
      ! ------------------
      do i = 1, nch 


           call VLIDORT_CoxMunk(SCAT%Surface,U10m(j),V10m(j),mr(i),solar_zenith (j),&
                                    sensor_zenith(j),relat_azymuth(j),.true.,rc)

           BRDF(j,i) = SCAT%Surface%Base%VIO%VBRDF_Sup_Out%BS_DBOUNCE_BRDFUNC(1,1,1,1)                                
           if ( rc /= 0 ) then
             BRDF(j,i) = MISSING
             return
           end if

      end do ! end loop over channels

      if ( verbose > 0 ) then
        if ( mod(j-1,1000) == 0 ) then
          print *, '<> VLIDORT BRDF: ', nint(j*100./nobs), '%'
        end if
      end if
     
  end do ! Loop over obs

end subroutine CoxMunk
