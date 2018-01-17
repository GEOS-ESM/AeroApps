module cloud_MieMod

use netcdf_Mod
implicit none

private

public :: getCOPvector

contains

!;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
! NAME
!    getCOPvector
! PURPOSE
!     reads the cloud optical properties file
! INPUT
!     dimname   : string of dimension name
!     filename  : file to be read
!     dimvar    : the variable to be read to
! OUTPUT
!     pmom for cloud particles
!  HISTORY
!
!;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
subroutine getCOPvector(filename, km, nobs, nch, nMom, nPol, idX, re, ssa, g, pmom, beta, truncation_factor)

    ! INPUT PARAMETERS:
    character(len=*), intent(in)          :: filename  ! cloud optics table file name
    integer, intent(in)                   :: km        ! number of vertical levels
    integer, intent(in)                   :: nobs      ! number of pixels
    integer, intent(in)                   :: nch       ! number of channels
    integer, intent(in)                   :: nMom      ! number of phase function moments
    integer, intent(in)                   :: nPol      ! number of components of the scattering matrix
    integer, intent(in)                   :: idX       ! channel index
    real, intent(in)                      :: re(km)    ! effective raidus [microns]

    ! OUTPUT PARAMETERS:
    real, intent(out)                     :: pmom(km,nch,nobs,nMom,nPol)
    real, intent(out)                     :: g(km,nch,nobs)
    real, intent(out)                     :: ssa(km,nch,nobs)
    real, intent(out)                     :: beta(km,nch,nobs)
    real, intent(out)                     :: truncation_factor(km,nch,nobs)

    integer                               :: maxStreams, maxRe, maxChannels
    real, allocatable                     :: pmomTable(:,:,:)
    real*8, allocatable                   :: gTable(:,:)
    real*8, allocatable                   :: ssaTable(:,:)
    real*8, allocatable                   :: radius(:)
    real*8, allocatable                   :: betaTable(:,:)
    real*8, allocatable                   :: tfTable(:,:)
    integer                               :: k, c, o, s 

    call readDim("streams", filename, maxStreams)
    call readDim("radius", filename, maxRe)
    call readDim("channels", filename, maxChannels)

    allocate( pmomTable(maxChannels, maxRe, maxStreams) )
    allocate( radius(maxRe) )
    allocate( gTable(maxChannels, maxRe) )
    allocate( ssaTable(maxChannels, maxRe) )
    allocate( betaTable(maxChannels, maxRe) )
    allocate( tfTable(maxChannels, maxRe) )

    call readvar3D("Legendre_Coefficients", filename, pmomTable)
    call readvar1D("Effective_Radius", filename, radius)
    call readvar2D("Asymmetry_Factor", filename, gTable)
    call readvar2D("Single_Scatter_Albedo", filename, ssaTable)
    call readvar2D("Extinction_Efficiency", filename, betaTable)
    call readvar2D("Truncation_Factor", filename, tfTable)

    ! convert legendre coefficients to be consistent with ours
    do s = 1, maxStreams
        pmomTable(:,:,s) = pmomTable(:,:,s)*(2*(s-1) + 1)
    end do


    ! Interpolate to effective radius
    pmom = 0
    do k = 1, km
        do c = 1, nch
            do o = 1, nobs
                g(k,c,o) = nn_interp(sngl(radius),sngl(gTable(idX,:)),re(k))
                ssa(k,c,o) = nn_interp(sngl(radius),sngl(ssaTable(idX,:)),re(k))
                beta(k,c,o) = nn_interp(sngl(radius),sngl(betaTable(idX,:)),re(k))/nn_interp(sngl(radius),sngl(betaTable(1,:)),re(k))
                truncation_factor(k,c,o) = nn_interp(sngl(radius),sngl(tfTable(idX,:)),re(k))
                do s = 1, maxStreams
                    pmom(k,c,o,s,1) = nn_interp(sngl(radius),pmomTable(idX,:,s),re(k))
                end do
            end do
        end do
    end do

end subroutine getCOPvector

!;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
! NAME
!     nn_interp
! PURPOSE
!     nearest neighbor interpolation
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
    real,intent(in),dimension(:)    :: x, y
    real,intent(in)                 :: xint
    real                            :: nn_interp
    integer                         :: below, above
    real                            :: top, bottom

    above = minloc(abs(xint - x), dim = 1, mask = (xint - x) .LT. 0)
    below = minloc(abs(xint - x), dim = 1, mask = (xint - x) .GE. 0)

    top = y(above) - y(below)
    bottom = x(above) - x(below)
    nn_interp = y(below) + (xint-x(below)) * top / bottom

end function nn_interp   

end module cloud_MieMod