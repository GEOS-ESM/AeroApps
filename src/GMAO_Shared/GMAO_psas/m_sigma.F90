!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !MODULE: m_sigma -
!
! !DESCRIPTION:
!
! !INTERFACE:

    module m_sigma
      implicit none
      private	! except

      public :: sigma_matvecpy
      public :: sigma_adjvec

      interface sigma_matvecpy; module procedure	&
	matvecpy_; end interface
      interface sigma_adjvec; module procedure	&
	adjvec_; end interface

! !REVISION HISTORY:
! 	25May00	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname='m_sigma'

contains

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
!
! !IROUTINE: matvecpy_ -  
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine matvecpy_(nav,krNav,ktNav,	&
	nsize,sigma,nvecs,Cvec,Rvec		)

      use m_Navigator,only : Navigator
      use m_Navigator,only : lsize
      use m_Navigator,only : get
      use m_die,  only : die,perr

      implicit none
      type(Navigator),intent(in) :: nav
      integer,dimension(:),intent(in) :: krNav
      integer,dimension(:),intent(in) :: ktNav

      integer,intent(in) :: nsize
      real   ,dimension(:),intent(in) :: sigma

      integer,intent(in) :: nvecs
      real   ,dimension(:,:),intent(in   ) :: Cvec
      real   ,dimension(:,:),intent(inout) :: Rvec

! !REVISION HISTORY:
! 	25May00	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code to replace mv_diag()
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::matvecpy_'

  integer :: i,ivec
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

  do inav=1,lsize(nav)
    call get(nav,inav,lc=lc,le=le)
    do ivec=1,nvecs
      Rvec(lc:le,ivec)=sigma(lc:le)*Cvec(lc:le,ivec)+Rvec(lc:le,ivec)
    end do
  end do

end subroutine matvecpy_
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
	nsize,sigma,nvecs,Cvec,Rvec		)

      use m_Navigator,only : Navigator
      use m_Navigator,only : lsize
      use m_Navigator,only : get
      use m_die,  only : die,perr

      implicit none
      type(Navigator),intent(in) :: nav
      integer,dimension(:),intent(in) :: krNav
      integer,dimension(:),intent(in) :: ktNav

      integer,intent(in) :: nsize
      real   ,dimension(:),intent(in) :: sigma

      integer,intent(in) :: nvecs
      real   ,dimension(:,:),intent(in ) :: Cvec
      real   ,dimension(:,:),intent(out) :: Rvec

! !REVISION HISTORY:
! 	25May00	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code to replace mv_diag()
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::adjvec_'

  integer :: i,ivec
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

  do inav=1,lsize(nav)
    call get(nav,inav,lc=lc,le=le)
    do ivec=1,nvecs
      Rvec(lc:le,ivec)=sigma(lc:le)*Cvec(lc:le,ivec)
    end do
  end do

end subroutine adjvec_

end module m_sigma
