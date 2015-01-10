!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !MODULE: m_sigmaPsi -
!
! !DESCRIPTION:
!
! !INTERFACE:

    module m_sigmaPsi
      implicit none
      private	! except

      public :: sigmaPsi_matvec
      public :: sigmaPsi_adjvec

      interface sigmaPsi_matvec; module procedure	&
	matvec_; end interface
      interface sigmaPsi_adjvec; module procedure	&
	adjvec_; end interface

! !REVISION HISTORY:
! 	25May00	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname='m_sigmaPsi'

contains

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
!
! !IROUTINE: matvec_ -  
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine matvec_(nav,krNav,ktNav,	&
	nsize,ktab,jtab,nvecs,Cvec,Rvec		)

      use m_Navigator,only : Navigator
      use m_Navigator,only : lsize
      use m_Navigator,only : get
      use m_die,  only : die,perr
      use FEsigW_imat,only : FEsigS_imat

      implicit none
      type(Navigator),intent(in) :: nav
      integer,dimension(:),intent(in) :: krNav
      integer,dimension(:),intent(in) :: ktNav

      integer,intent(in) :: nsize
      integer,dimension(:),intent(in) :: ktab
      integer,dimension(:),intent(in) :: jtab

      integer,intent(in) :: nvecs
      real   ,dimension(:,:),intent(in ) :: Cvec
      real   ,dimension(:,:),intent(out) :: Rvec

! !REVISION HISTORY:
! 	25May00	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code to replace mv_diag()
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::matvec_'
  include "ktmax.h"

  real,allocatable,dimension(:) :: sigPsi

  integer :: kt,i,ivec
  integer :: ier
  integer :: lc,le
  integer :: inav

  logical :: invalid	! a condition flag

!-----------------------------------------------------------------------
  
	! Dimension consistency checking

  invalid=.false.
  if(nsize>size(Cvec,1)) then
    call perr(myname_,'nsize',nsize,'size(Cvec,1)',size(Cvec,1))
    invalid=.true.
  endif

  if(nsize>size(Rvec,1)) then
    call perr(myname_,'nsize',nsize,'size(Rvec,1)',size(Rvec,1))
    invalid=.true.
  endif

  if(nvecs>size(Cvec,2)) then
    call perr(myname_,'nvecs',nvecs,'size(Cvec,2)',size(Cvec,2))
    invalid=.true.
  endif

  if(nvecs>size(Rvec,2)) then
    call perr(myname_,'nvecs',nvecs,'size(Rvec,2)',size(Rvec,2))
    invalid=.true.
  endif

  if(invalid) call die(myname_)

!-----------------------------------------------------------------------
	! loop over all Rvec partitions

	allocate(sigPsi(nsize),stat=ier)
		if(ier/=0) call die(myname_,'allocate()',ier)

  do inav=1,lsize(nav)
    call get(nav,inav,lc=lc,le=le)
    kt=ktNav(inav)

    select case(kt)
    case (ktUU,ktVV,ktUs,ktVs)

      sigPsi(lc:le)=(/ (FEsigS_imat(ktab(i),jtab(i)), i=lc,le) /)

      do ivec=1,nvecs
	Rvec(lc:le,ivec)=sigPsi(lc:le)*Cvec(lc:le,ivec)
      end do

    case default
      do ivec=1,nvecs
	Rvec(lc:le,ivec)=0.
      end do
    end select
  end do

	deallocate(sigPsi,stat=ier)
		if(ier/=0) call die(myname_,'allocate()',ier)

end subroutine matvec_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
!
! !IROUTINE: adjvec_ -
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine adjvec_(nav,krNav,ktNav,	&
	nsize,ktab,jtab,nvecs,Rvec,Cvec		)

      use m_Navigator,only : Navigator
      use m_Navigator,only : lsize
      use m_Navigator,only : get
      use m_die,  only : die,perr
      use FEsigW_imat,only : FEsigS_imat

      implicit none
      type(Navigator),intent(in) :: nav
      integer,dimension(:),intent(in) :: krNav
      integer,dimension(:),intent(in) :: ktNav

      integer,intent(in) :: nsize
      integer,dimension(:),intent(in) :: ktab
      integer,dimension(:),intent(in) :: jtab

      integer,intent(in) :: nvecs
      real   ,dimension(:,:),intent(in ) :: Rvec
      real   ,dimension(:,:),intent(out) :: Cvec

! !REVISION HISTORY:
! 	25May00	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::adjvec_'
  include "ktmax.h"

  real,allocatable,dimension(:) :: sigPsi

  integer :: kt,i,ivec
  integer :: ier
  integer :: lc,le
  integer :: inav

  logical :: invalid	! a condition flag

!-----------------------------------------------------------------------
  
	! Dimension consistency checking

  invalid=.false.
  if(nsize>size(Cvec,1)) then
    call perr(myname_,'nsize',nsize,'size(Cvec,1)',size(Cvec,1))
    invalid=.true.
  endif

  if(nsize>size(Rvec,1)) then
    call perr(myname_,'nsize',nsize,'size(Rvec,1)',size(Rvec,1))
    invalid=.true.
  endif

  if(nvecs>size(Cvec,2)) then
    call perr(myname_,'nvecs',nvecs,'size(Cvec,2)',size(Cvec,2))
    invalid=.true.
  endif

  if(nvecs>size(Rvec,2)) then
    call perr(myname_,'nvecs',nvecs,'size(Rvec,2)',size(Rvec,2))
    invalid=.true.
  endif

  if(invalid) call die(myname_)

!-----------------------------------------------------------------------
	! loop over all Rvec partitions

	allocate(sigPsi(nsize),stat=ier)
		if(ier/=0) call die(myname_,'allocate()',ier)

  do inav=1,lsize(nav)
    call get(nav,inav,lc=lc,le=le)
    kt=ktNav(inav)

    select case(kt)
    case (ktUU,ktVV,ktUs,ktVs)

      sigPsi(lc:le)=(/ (FEsigS_imat(ktab(i),jtab(i)), i=lc,le) /)

      do ivec=1,nvecs
	Cvec(lc:le,ivec)=sigPsi(lc:le)*Rvec(lc:le,ivec)
      end do

    case default
      do ivec=1,nvecs
	Cvec(lc:le,ivec)=0.
      end do
    end select
  end do

	deallocate(sigPsi,stat=ier)
		if(ier/=0) call die(myname_,'allocate()',ier)

end subroutine adjvec_

end module m_sigmaPsi
