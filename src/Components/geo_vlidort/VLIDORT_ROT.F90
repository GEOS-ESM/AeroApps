module VLIDORT_ROT
!
!  Subroutines to call the VLIDORT ROT calculator
!.............................................................................
  implicit NONE

  PUBLIC VLIDORT_ROT_CALC

  contains

  logical function IS_MISSING(x,MISSING)
    real*8, intent(in)     :: x
    real*8, intent(in)     :: MISSING

    IS_MISSING = abs(x/MISSING-1)<0.001
    return
  end function IS_MISSING

  subroutine VLIDORT_ROT_CALC(km,nch,nobs,channels,pe,he,te,MISSING,verbose,ROT,rc)  
  !
  ! Uses VLIDORT_Rayleigh subroutine to calc ROT
  !
    use VLIDORT_ScatMod

    implicit NONE

    !INPUT PARAMETERS:

    integer,          intent(in)  :: km    ! number of levels on file
    integer,          intent(in)  :: nch   ! number of channels
    integer,          intent(in)  :: nobs  ! number of observations
                                        
    real*8, target,   intent(in)  :: channels(nch)    ! wavelengths [nm]

    real*8, target,   intent(in)  :: pe(km+1,nobs)    ! pressure at layer edges [Pa]
    real*8, target,   intent(in)  :: he(km+1,nobs)    ! height above sea-level  [m]
    real*8, target,   intent(in)  :: te(km+1,nobs)    ! temperature at layer edges [K]
    
    real*8, target,   intent(in)  :: MISSING          ! MISSING VALUE                                      
    integer,          intent(in)  :: verbose

    !OUTPUT PARAMETERS
    integer,          intent(out) :: rc                          ! return code
    real*8,           intent(out) :: ROT(km,nobs,nch)                 ! rayleigh optical thickness 

    !

    integer                       :: i,j, ier
   
    type(VLIDORT_scat)            :: SCAT
    type(VLIDORT_output_scalar)   :: output  

    rc = 0
    ier = 0

    call VLIDORT_Init( SCAT%Surface%Base, km, rc)
    if ( rc /= 0 ) return

    do j = 1, nobs
        SCAT%pe => pe(:,j)
        SCAT%ze => he(:,j)
        SCAT%te => te(:,j) 
        ! Loop over channels
        ! ------------------
        do i = 1, nch 
          SCAT%wavelength = channels(i)
          call VLIDORT_Rayleigh(SCAT,ier)
          ROT(:,j,i) = SCAT%rot 

          if ( ier /= 0 ) then
            ROT(:,j,i) = MISSING
            cycle
          end if
        end do ! end loop over channels

        if ( verbose > 0 ) then
          if ( mod(j-1,1000) == 0 ) then
            print *, '<> VLIDORT ROT Calculator: ', nint(j*100./nobs), '%'
          end if
        end if

    end do ! Loop over obs
  end subroutine VLIDORT_ROT_CALC
end module VLIDORT_ROT