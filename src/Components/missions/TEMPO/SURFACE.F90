module SURFACE   
!
!  Calculates surface reflectance for different models
!
!.............................................................................
  implicit NONE

  PUBLIC SURFACE_LandMODIS

  contains

  logical function IS_MISSING(x,MISSING)
    real*8, intent(in)     :: x
    real*8, intent(in)     :: MISSING

    IS_MISSING = abs(x/MISSING-1)<0.001
    return
  end function IS_MISSING

  subroutine SURFACE_LandMODIS (km, nch, nobs, channels,                    &
                                kernel_wt, param,                           &
                                solar_zenith, relat_azymuth, sensor_zenith, &
                                MISSING, verbose, BR, Kvol, Kgeo, rc )
  !
  ! Computes surface reflectance
  ! from MODIS RTLS kernel weights
  !


    implicit NONE

    logical, parameter            :: scalar = .true.
    integer, parameter            :: nkernel = 3
    integer, parameter            :: nparam  = 2
    real*8, parameter             :: pi = 4.*atan(1)
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
    real*8,           intent(out) :: Kvol(nobs,nch)                    ! White-sky albedo
    real*8,           intent(out) :: Kgeo(nobs,nch)                    ! Black-sky albedo        
    real*8,           intent(out) :: BR(nobs,nch)                     ! bi-directional reflectance 
  !                         ---  
    integer                       :: i,j,n,p,ier
    real*8                        :: sza, vza, raa
    real*8                        :: fiso, fgeo, fvol
    real*8                        :: h_b, b_r
    real*8                        :: delta, D
    real*8                        :: cost, t, O


  
    rc = 0
    ier = 0

    do j = 1,nobs

       ! Make sure angles are available
       ! -----------------------------------------
      if ( IS_MISSING(solar_zenith(j),MISSING)  .OR. & 
           IS_MISSING(sensor_zenith(j),MISSING) .OR. &
           IS_MISSING(relat_azymuth(j),MISSING)  )  then

        Kvolj,:) = MISSING
        Kgeo(j,:) = MISSING
        BR(j,:)  = MISSING
        cycle
      end if

      vza = sensor_zenith(j)
      sza = solar_zenith(j)
      raa = relat_azymuth(j)
           
      ! Loop over channels
      ! ------------------
      do i = 1, nch 
        ! Make sure kernel weights and parameters are defined
        do n = 1, nkernel
          if (IS_MISSING(kernel_wt(n,i,j),MISSING)) then
            Kvol(j,:) = MISSING
            Kgeo(j,:) = MISSING
            BR(j,:)   = MISSING
            cycle
          end if
        end do
        do p = 1, nparam
          if (IS_MISSING(param(p,i,j),MISSING)) then
            BR(j,:)    = MISSING
            Kvol(j,:)  = MISSING
            Kgeo(j,:)  = MISSING
            cycle
          end if
        end do

        fiso = kernel_wt(1,i,j)
        fgeo = kernel_wt(2,i,j)
        fvol = kernel_wt(3,i,j)

        h_b  = param(1,i,j)
        b_r  = param(2,i,j)

        if ( verbose > 0 ) then
          print*, 'DO MODIS BRDF'
        end if


        delta = acos(cos(sza)*cos(vza) + sin(sza)*sin(vza)*cos(raa))
        D     = sqrt(tan(sza)*tan(sza) + tan(vza)*tan(vza) - 2*tan(sza)*tan(vza)*cos(raa))
        ! ROSS-THICK (KVOL)
        !-----------------------
        Kvol(j,i) = ( ((0.5*pi - delta)*cos(delta) + sin(delta) ) / (cos(sza) + cos(vza)) ) - 0.25*pi

        ! Li-Sparse (KGEO)
        !-----------------------
        sza = atan(b_r*tan(sza))
        vza = atan(b_r*tan(vza))

        delta = acos(cos(sza)*cos(vza) + sin(sza)*sin(vza)*cos(raa))
        D     = sqrt(tan(sza)*tan(sza) + tan(vza)*tan(vza) - 2*tan(sza)*tan(vza)*cos(raa))

        cost = h_b*sqrt( D*D + (tan(sza)*tan(vza)*sin(raa))^2 ) / ( sec(sza) + sec(vza) )
        cost = min(1.0,cost)

        t    = acos(cost)

        O    = ( t - sin(t)*cost )*( sec(sza) + sec(vza) ) / pi

        Kgeo(j,i) = O - sec(sza) - sec(vza) + 0.5*sec(vza)*( 1 + cos(delta) )

        BR(j,i) = fiso + fgeo*Kgeo(j,i) + fvol*Kvol(j,i)       

        if ( verbose > 0 ) then
          print *, 'My reflectance land modis',BR(j,i), Kgeo(j,i), Kvol(j,i) 
        end if

      end do ! end loop over channels
     
      if ( verbose > 0 ) then
        if ( mod(j-1,1000) == 0 ) then
          print *, '<> VLIDORT Scalar: ', nint(j*100./nobs), '%'
        end if
      end if

    end do ! Loop over obs

  end subroutine SURFACE_LandMODIS

  !..........................................................................


end module SURFACE