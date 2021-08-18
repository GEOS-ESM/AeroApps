module cirrus_MieMod

  use netcdf_mod
  implicit none

  private

  public cirrusMieTable
  public get_cirrusOP

! Cirrus LUT
  type cirrusMieTable
    real*8, pointer :: channels(:)
    real*8, pointer :: reff(:)
    real*8, pointer :: ssa(:,:)
    real*8, pointer :: ext(:,:)
    real*8, pointer :: pmom(:,:,:,:,:)
    integer         :: nch
    integer         :: nreff

  end type cirrusMieTable


  contains

  subroutine get_cirrusOP(filename,nobs,nch,nMom,nPol,channels,tau532,reff,tau,ssa,pmom)
    ! INPUTS
    character(len=*), intent(in)         :: filename  ! cirrus properties table file name
    integer, intent(in)                  :: nobs      ! number of pixels
    integer, intent(in)                  :: nch       ! number of channels
    integer, intent(in)                  :: nMom      ! number of phase function moments
    integer, intent(in)                  :: nPol      ! number of components of the scattering matrix
    real*8, intent(in)                   :: channels(nch)  ! channels to interpolate to
    real*8, intent(in)                   :: tau532(nobs)   ! cirrus OD @ 532 nm
    real*8, intent(in)                   :: reff(nobs)     ! effective radius

    ! OUTPUTS
    real*8, intent(out)                  :: tau(nch,nobs)
    real*8, intent(out)                  :: ssa(nch,nobs)
    real*8, intent(out)                  :: pmom(nch,nobs,nMom,nPol)

    integer                              :: tnMom, nreff, nlam
    real*8,allocatable                   :: ssaT(:,:)
    real*8,allocatable                   :: extT(:,:)
    real*8,allocatable                   :: beta(:,:)
    real*8                               :: ext532
    real*8,allocatable                   :: pmomT(:,:,:,:,:)
    real*8,allocatable                   :: reffT(:)
    real*8,allocatable                   :: channelsT(:)

    real*8,allocatable                   :: ssaCH(:)
    real*8,allocatable                   :: betaCH(:)
    real*8,allocatable                   :: pmomCH(:,:,:)

    integer                              :: r,idX,o,p,s,c

    ! get dimensions of Table
    call readDim('nlam',filename,nlam)
    call readDim('nreff',filename,nreff)
    call readDim('nmommax',filename,tnMom)

    ! allocate arrays
    allocate(channelsT(nlam))
    allocate(reffT(nreff))
    allocate(ssaT(nlam,nreff))
    allocate(extT(nlam,nreff))
    allocate(beta(nlam,nreff))
    allocate(pmomT(nPol,tnMom,nreff,1,nlam))

    allocate(ssaCH(nlam))
    allocate(betaCH(nlam))
    allocate(pmomCH(nPol,tnMom,nlam))
    ssaCH = 0
    betaCH = 0
    pmomCH = 0

    call readvar1D("wavelen",filename,channelsT)
    call readvar1D("reff",filename,reffT)
    call readvar2D("ssa",filename,ssaT)
    call readvar2D("ext",filename,extT)
    call readvar5D("pmom",filename,pmomT)

    ! convert microns to nm
    channelsT = 1d3*channelsT

    ! calculate beta = ext/ext532
    do r = 1, nreff
        ext532 = nn_interp(channelsT,extT(:,r),532d0)
        beta(:,r) = extT(:,r)/ext532
    end do

    ! interpolate to effective raidus and wavelength
    tau = 0.0
    ssa = 0.0
    pmom = 0.0
    do o=1,nobs
        do idX = 1, nlam
            ! first interpolate along effective radius
            ssaCH(idX) = nn_interp(reffT,ssaT(idX,:),reff(o))
            betaCH(idX) = nn_interp(reffT,beta(idX,:),reff(o))
            do p = 1, nPol
                do s = 1, tnMom
                    pmomCH(p,s,idX) = nn_interp(reffT,pmomT(p,s,:,1,idX),reff(o))
                end do
            end do
        end do

        ! interpolate along channel
        do c = 1, nch
            ssa(c,o) = nn_interp(channelsT,ssaCH,channels(c))
            tau(c,o) = tau532(o)*nn_interp(channelsT,betaCH,channels(c))
            do p = 1, nPol
                do s = 1, tnMom
                    pmom(c,o,s,p) = nn_interp(channelsT,pmomCH(p,s,:),channels(c))
                end do
            end do
        end do
    end do
  end subroutine get_cirrusOP

!;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
! NAME
!     nn_interp
! PURPOSE
!     linear interpolation between nearest neighbors
! INPUT
!     x     : x-values
!     y     : y-values
!     xint  : x-value to interpolate to
! OUTPUT
!     nn_interp  : interpolated y-value
!  HISTORY
!     10 June 2015 P. Castellanos
!;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

function nn_interp(x,y,xint)
    real*8,intent(in),dimension(:)    :: x, y
    real*8,intent(in)                 :: xint
    real*8                            :: nn_interp
    integer                           :: below, above
    real*8                            :: top, bottom

    if (xint .GE. maxval(x)) then
        below = maxloc(x, dim = 1)
        nn_interp = y(below)
    else if (xint .LT. minval(x)) then
        above = minloc(x, dim = 1)
        nn_interp = y(above)
    else
        above = minloc(abs(xint - x), dim = 1, mask = (xint - x) .LT. 0)
        below = minloc(abs(xint - x), dim = 1, mask = (xint - x) .GE. 0)

        top = y(above) - y(below)
        bottom = x(above) - x(below)
        nn_interp = y(below) + (xint-x(below)) * top / bottom
    end if

end function nn_interp

end module cirrus_MieMod






