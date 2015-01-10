!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !MODULE: m_ObsErrCovMatx -
!
! !DESCRIPTION:
!
! !INTERFACE:

    module m_ObsErrCovMatx
      use m_SparseComm, only : SparseComm
      use m_AttrVect,only : AttrVect
      implicit none
      private	! except

      public :: ObsErrCovMatx		! The class data structure
      public :: ObsErrCovMatx_Cxpy
      public :: ObsErrCovMatx_Cx
      public :: ObsErrCovMatx_init
      public :: ObsErrCovMatx_clean
      public :: ObsErrCovMatx_showCosts
      public :: kind_OECov
      public :: kind_covOu
      public :: kind_covOc

    type ObsErrCovMatx
      integer :: kind_cov		! a place holder for now
      logical :: symmetric
      real :: myCostOc
      real :: myCostOu
      type(AttrVect) :: mNavOc
      type(AttrVect) :: mNavOu

      type(SparseComm) :: spGathOc
      type(SparseComm) :: spScatOc
      type(SparseComm) :: spGathOu
      type(SparseComm) :: spScatOu

    end type ObsErrCovMatx

      interface ObsErrCovMatx_Cxpy; module procedure	&
	symCxpy_; end interface
      interface ObsErrCovMatx_Cx; module procedure	&
	symCx_; end interface
      interface ObsErrCovMatx_init; module procedure	&
	symSched_; end interface
      interface ObsErrCovMatx_clean; module procedure	&
	clean_; end interface
      interface ObsErrCovMatx_showCosts; module procedure	&
	showCosts_; end interface

! !REVISION HISTORY:
! 	25May00	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname='m_ObsErrCovMatx'
  include "kind_covs.h"
  integer,parameter :: kind_covOu=kind_covU
  integer,parameter :: kind_covOc=kind_covC
  integer,parameter :: kind_OECov=ior(kind_covOu,kind_covOc)

#include "assert.H"
contains

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: symCx_ - operator y=R(x), R is Obs.Err.Cov.
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine symCx_(oecov,comm, nav,krNav,ktNav,	&
	sigOu,sigOc,corObs, x, y)

      use m_Navigator,only : Navigator
      use m_CorAttrX ,only : CorAttrX
      implicit none

      type(ObsErrCovMatx),intent(inout) :: oecov
      integer,intent(in) :: comm	! message-passing communicator

      type(Navigator),intent(in) :: nav
      integer,dimension(:),intent(in) :: krNav
      integer,dimension(:),intent(in) :: ktNav

      real,dimension(:),intent(in) :: sigOu
      real,dimension(:),intent(in) :: sigOc
      type(CorAttrX)   ,intent(in) :: corObs	! now for bothr Cc + Cu

      real   ,dimension(:,:),intent(in ) :: x
      real   ,dimension(:,:),intent(out) :: y

! !REVISION HISTORY:
! 	25May00	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::symCx_'
  integer :: inav
  integer :: lc,le

	! Initialize all vector segments of y to zero

  y(:,:)=0.

	! y = Cx + y

  call symCxpy_(oecov,comm, nav,krNav,ktNav,	&
	sigOu,sigOc,corObs, x, y)

end subroutine symCx_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: symCxpy_ - operator y=R(x)+y, R is Obs.Err.Cov.
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine symCxpy_(oecov,comm, nav,krNav,ktNav,	&
	sigOu,sigOc,corObs, x, y)

      use m_sigma,only : sigma_matvecpy
      use m_sigma,only : sigma_adjvec
      use m_CorAttrX,only : CorAttrX
      use m_CorAttrX,only : ptr_Navigator
      use m_CorAttrX,only : ptr_krNav
      use m_CorAttrX,only : ptr_ktNav
      use m_CorAttrX,only : ptr_kx
      use m_CorAttrX,only : ptr_ks
      use m_CorAttrX,only : ptr_qr
      use m_CorAttrX,only : ptr_kl
      use m_CorAttrX,only : ptr_collVec

      use m_CorMatxO,only : symCx
      use m_CorMatxU,only : symCx

      use m_Collector,only : Collector
      use m_Collector,only : localSize

      use m_die  ,only : die,assert_
      use m_mall ,only : mall_ison,mall_mci,mall_mco
      use m_Navigator,only : Navigator
      use m_Navigator,only : lsize

      implicit none

      type(ObsErrCovMatx),intent(inout) :: oecov
      integer,intent(in) :: comm

      type(Navigator),intent(in) :: nav
      integer,dimension(:),intent(in) :: krNav
      integer,dimension(:),intent(in) :: ktNav

      real,dimension(:),intent(in) :: sigOu
      real,dimension(:),intent(in) :: sigOc
      type(CorAttrX)   ,intent(in) :: corObs	! now for bothr Cc + Cu

      real   ,dimension(:,:),intent(in   ) :: x
      real   ,dimension(:,:),intent(inout) :: y

! !REVISION HISTORY:
! 	25May00	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::symCxpy_'
  integer :: ier
  real,allocatable,dimension(:,:) :: tx,ty
  integer :: nvecs,nvi,nvseg
  type(Collector),pointer :: collvec

	collvec => ptr_collvec(corObs)

  nvi=localSize(collvec)
	ASSERT( nvi<=size(sigOu) )
	ASSERT( nvi<=size(sigOc) )
	ASSERT( nvi<=size(x,1) )
	ASSERT( nvi==size(y,1) )

	nullify(collvec)

  nvecs=size(x,2)
	ASSERT( nvecs==size(y,2) )

  nvseg=lsize(nav)
	ASSERT( nvseg<=size(krNav) )
	ASSERT( nvseg<=size(ktNav) )

	allocate(tx(nvi,nvecs),ty(nvi,nvecs),stat=ier)
		if(ier/=0) call die(myname_,'allocate()',ier)
		if(mall_ison()) then
		  call mall_mci(tx,myname)
		  call mall_mci(ty,myname)
		endif

  if(iand(oecov%kind_cov,kind_covOu)/=0) then

	! tx = sigOu'(x)

    call sigma_adjvec(nav,krNav,ktNav, nvi,sigOu,nvecs,x,tx)

	! ty = corOu(tx)

    call symCx(oecov%mNavOu,oecov%spGathOu,oecov%spScatOu, &
	 kind_covOu, comm,	&
	ptr_Navigator(corObs,whole=.true.),	&
	    ptr_krNav(corObs,whole=.true.),	&
	    ptr_ktNav(corObs,whole=.true.),	&
	       ptr_kx(corObs,whole=.true.),	&
	       ptr_ks(corObs,whole=.true.),	&
	       ptr_kl(corObs,whole=.true.),	&
	ptr_collvec(corObs),tx,ty, stat=ier)

	if(ier/=0) call die(myname_,'symCx(corOu)',ier)

	! y = y+sigOu(ty)

    call sigma_matvecpy(nav,krNav,ktNav, nvi,sigOu,nvecs,ty,y)

  endif

  if(iand(oecov%kind_cov,kind_covOc)/=0) then

	! tx = sigOc'(x)

    call sigma_adjvec(nav,krNav,ktNav, nvi,sigOc,nvecs,x,tx)

	! ty = corOc(tx)

    call symCx(oecov%mNavOc,oecov%spGathOc,oecov%spScatOc, &
	 kind_covOc, comm,	&
	ptr_Navigator(corObs,whole=.true.),	&
	    ptr_krNav(corObs,whole=.true.),	&
	    ptr_ktNav(corObs,whole=.true.),	&
	       ptr_kx(corObs,whole=.true.),	&
	       ptr_qr(corObs,whole=.true.),	&
	       ptr_kl(corObs,whole=.true.),	&
	ptr_collvec(corObs),tx,ty, stat=ier)

	if(ier/=0) call die(myname_,'symCx(corOc)',ier)

	! y = y+sigOc(ty)

    call sigma_matvecpy(nav,krNav,ktNav, nvi,sigOc,nvecs,ty,y)

  endif

		if(mall_ison()) then
		  call mall_mco(tx,myname)
		  call mall_mco(ty,myname)
		endif
	deallocate(tx,ty,stat=ier)
		if(ier/=0) call die(myname_,'deallocate()',ier)

end subroutine symCxpy_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: symSched_ - initialize an object
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine symSched_(oecov,corObs,comm,kind_mat,kind_cov)
      use m_symMatx ,only : symMatx_schedule
      use m_CorAttrX,only : CorAttrX
      use m_CorAttrX,only : ptr_Navigator
      use m_CorAttrX,only : ptr_krNav
      use m_CorAttrX,only : ptr_ktNav
      use m_CorAttrX,only : ptr_collNav
      use m_Collector,only : Collector
      use m_sparse, only : sparse

      implicit none

      type(ObsErrCovMatx),intent(out) :: oecov
      type(CorAttrX),intent(in) :: corObs
      integer,intent(in) :: comm
      integer,intent(in) :: kind_mat
      integer,optional,intent(in) :: kind_cov

! !REVISION HISTORY:
!	06Sep00	- Jing Guo
!		. modified from init_() to symSched_()
! 	25May00	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::symSched_'


  Type(Collector), Pointer :: collnav

  oecov%kind_cov = kind_OECov
  if(present(kind_cov)) oecov%kind_cov=iand(kind_cov,kind_OECov)

  collnav => ptr_collNav(corObs)

  if(iand(oecov%kind_cov,kind_covOu) /=0) then

    call symMatx_schedule(oecov%mNavOu,		  &
	 oecov%spGathOu, oecov%spScatOu, collnav, &
	ptr_Navigator(corObs,whole=.true.),	  &
	    ptr_krNav(corObs,whole=.true.),	  &
	    ptr_ktNav(corObs,whole=.true.),	  &
	comm,kind_covOu,kind_mat,sparse,oecov%myCostOu)
  endif

  if(iand(oecov%kind_cov,kind_covOc) /=0) then

    call symMatx_schedule(oecov%mNavOc,		  &
	 oecov%spGathOc, oecov%spScatOc, collnav, &
	ptr_Navigator(corObs,whole=.true.),	  &
	    ptr_krNav(corObs,whole=.true.),	  &
	    ptr_ktNav(corObs,whole=.true.),	  &
	comm,kind_covOc,kind_mat,sparse,oecov%myCostOc)
  endif

  Nullify(collnav)

  oecov%symmetric=.true.

end subroutine symSched_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: clean_ - clean an object
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine clean_(oecov)
      use m_AttrVect,only : clean
      implicit none
      type(ObsErrCovMatx),intent(inout) :: oecov

! !REVISION HISTORY:
! 	25May00	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::clean_'

  if(iand(oecov%kind_cov,kind_covOu)/=0) call clean(oecov%mNavOu)
  if(iand(oecov%kind_cov,kind_covOc)/=0) call clean(oecov%mNavOc)

  oecov%myCostOu=0.
  oecov%myCostOc=0.
  oecov%kind_cov=0
  oecov%symmetric=.false.

end subroutine clean_
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: showCosts_ - show the estimated costs of matrix operators
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine showCosts_(oecov,lu,root,comm,listall)
      use m_showDistrib,only : showDistrib
      implicit none
      type(ObsErrCovMatx),intent(in) :: oecov	! An oprator
      integer         ,intent(in) :: lu		! significant on root
      integer         ,intent(in) :: root	! root PE ID
      integer         ,intent(in) :: comm	! communicator
      logical,optional,intent(in) :: listall


! !REVISION HISTORY:
! 	05Mar01	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::showCosts_'

  call showDistrib(lu,myname//'::Ru',oecov%myCostOu,root,comm,	&
	listall=listall)
  call showDistrib(lu,myname//'::Rc',oecov%myCostOc,root,comm,	&
	listall=listall)

end subroutine showCosts_
end module m_ObsErrCovMatx
