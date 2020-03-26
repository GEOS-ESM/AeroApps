!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !MODULE: m_llgrad - compute gradient of a lat.-long. gridded field
!
! !DESCRIPTION:
!
! !INTERFACE:

    module m_llgrad
      implicit none
      private	! except

      public :: llgrad

      interface llgrad; module procedure grad_; end interface

! !REVISION HISTORY:
! 	24Apr00	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- Documented for further development.
!		- Modified the algorithm for polar vectors.
! 	24Mar98	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname='m_llgrad'

contains
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: grad_ - compute gridded gradient
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine grad_(h,hm,hl)
      implicit none
      real,dimension(:,:),intent(in)  :: h
      real,dimension(:,:),intent(out) :: hm,hl	! the gradient

! !REVISION HISTORY:
! 	25Mar98	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::grad_'
  real :: deg
  real :: dlat,dlon
  real :: clat,delm,dell
  real :: wx,wy
  integer :: i,j
  integer :: im,jm

	! Preconditions: [shape(hm), shape(hl)] >= shape(h), etc.

  im=size(h,1)
  jm=size(h,2)

  deg=4.*atan(1.)/180.

  dlat=180./(jm-1)	! this is not a periodic grid
  dlon=360./im		! this is a periodic grid

	! Finite differencing everywhere except at the South and the
	! North Poles.

  do j=2,jm-1
    clat=cos(((j-1)*dlat-90.)*deg)
    delm=2.*dlat*deg
    dell=2.*dlon*deg*clat

    do i=1,im
      hm(i,j)=(h(i,j+1)-h(i,j-1))/delm
    end do

    hl(  1,j)=(h(2  ,j)-h(im  ,j))/dell
    do i=2,im-1
      hl(i,j)=(h(i+1,j)-h(i-1 ,j))/dell
    end do
    hl( im,j)=(h(1  ,j)-h(im-1,j))/dell
  end do

	! At a pole,
	!
	!   1) compute the averaged "angular-vector" at the poles,
	!     using vectors at the pole's "ring" grid points;
	!   2) compute the vectors at the polar grid points using the
	!     averaged "angular-vector".

	! At the South Pole

  call meanvect_( hl(1:im,2),hm(1:im,2),-90.+dlat,dlon,deg,wx,wy)
  call polarvect_(hl(1:im,1),hm(1:im,1),-90.     ,dlon,deg,wx,wy)

	! At the North Pole

  call meanvect_( hl(1:im,jm-1),hm(1:im,jm-1),+90.-dlat,dlon,deg,wx,wy)
  call polarvect_(hl(1:im,jm  ),hm(1:im,jm  ),+90.     ,dlon,deg,wx,wy)

end subroutine grad_
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: meanvect_ - compute mean-vector arround a pole
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine meanvect_(vl,vm,rlat,dlon,deg,wx,wy)

      implicit none
      real,dimension(:),intent(in) :: vl,vm
      real,intent(in)  :: rlat,dlon,deg
      real,intent(out) :: wx,wy	! ,wz is not needed for a polar-mean

! !REVISION HISTORY:
! 	24Apr00	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::meanvect_'

  real :: sinlat
  real :: rlon
  real :: coslon,sinlon
  real :: emx,emy
  real :: elx,ely
  integer :: i,im

	! Precondition: size(vl) == size(vm)

  im=size(vl)

  sinlat=sin(rlat*deg)
  wx=0.
  wy=0.
  do i=1,im
    rlon=(i-1)*dlon-180.
    coslon=cos(rlon*deg)
    sinlon=sin(rlon*deg)

    emx=-coslon*sinlat
    emy=-sinlon*sinlat
    elx=-sinlon
    ely= coslon

    wx=wx + vl(i)*emx - vm(i)*elx
    wy=wy + vl(i)*emy - vm(i)*ely

  end do
  wx=wx/im
  wy=wy/im

end subroutine meanvect_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: polarvect_ - compute a vector at a pole
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine polarvect_(vl,vm,rlat,dlon,deg,wx,wy)

      implicit none
      real,dimension(:),intent(out) :: vl,vm
      real,intent(in) :: rlat,dlon,deg
      real,intent(in) :: wx,wy	! ,wz is not needed for a polar-mean

! !REVISION HISTORY:
! 	24Apr00	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::polarvect_'

  real :: sinlat
  real :: rlon
  real :: coslon,sinlon
  real :: emx,emy
  real :: elx,ely
  integer :: i,im

  im=size(vl)

	! Note that the z component of w, wz, is not included
	! because emz = cos(lat) == 0 at the poles.

  sinlat=sign(1.,rlat)	! sin(lat) is either 1. (at the North pole) or
			! -1. (at the South pole)
  do i=1,im
    rlon=(i-1)*dlon-180.
    coslon=cos(rlon*deg)
    sinlon=sin(rlon*deg)

    emx=-coslon*sinlat
    emy=-sinlon*sinlat
    elx=-sinlon
    ely= coslon

    vl(i) =   wx*emx + wy*emy
    vm(i) = -(wx*elx + wy*ely)

  end do
end subroutine polarvect_
end module m_llgrad
