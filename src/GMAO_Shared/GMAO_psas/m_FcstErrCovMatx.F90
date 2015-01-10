!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !MODULE: m_FcstErrCovMatx -
!
! !DESCRIPTION:
!
!     Note that bue to the symmetry of a covariance matrix operator,
!   the adjoint of the operator is abtained simply by switching the
!   row arguments and the column arguments of this subroutine.
!
! !INTERFACE:

    module m_FcstErrCovMatx
      use m_AttrVect,only : AttrVect
      use m_SparseComm, only: SparseComm
      implicit none
      private	! except

      public :: FcstErrCovMatx		! The class data structure
      public :: FcstErrCovMatx_Cxpy
      public :: FcstErrCovMatx_Cx
      public :: FcstErrCovMatx_init
      public :: FcstErrCovMatx_clean
      public :: FcstErrCovMatx_showCosts
      public :: kind_FECov
      public :: kind_covPhi
      public :: kind_covPsi
      public :: kind_covChi

    type FcstErrCovMatx
      integer :: kind_cov		! a place holder for now
      logical :: symmetric
      real :: myCostPhi
      real :: myCostPsi
      real :: myCostChi
      type(AttrVect) :: mNavPhi
      type(AttrVect) :: mNavPsi
      type(AttrVect) :: mNavChi
      type(SparseComm) :: spGathPhi
      type(SparseComm) :: spScatPhi
      type(SparseComm) :: spGathPsi
      type(SparseComm) :: spScatPsi
      type(SparseComm) :: spGathChi
      type(SparseComm) :: spScatChi
    end type FcstErrCovMatx

      interface FcstErrCovMatx_Cxpy; module procedure	&
	symCxpy_,	&
	recCxpy_
      end interface

      interface FcstErrCovMatx_Cx; module procedure	&
	symCx_,		&
	recCx_
      end interface

      interface FcstErrCovMatx_init; module procedure	&
	symSched_,	&
	recSched_
      end interface

      interface FcstErrCovMatx_clean; module procedure	&
	clean_; end interface
      interface FcstErrCovMatx_showCosts; module procedure	&
	showCosts_; end interface

! !REVISION HISTORY:
! 	25May00	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname='m_FcstErrCovMatx'
  include "kind_covs.h"
  integer,parameter :: kind_covPhi=kind_covF
  integer,parameter :: kind_covPsi=kind_covS
  integer,parameter :: kind_covChi=kind_covV
  integer,parameter :: kind_FECov=ior(kind_covPhi,	&
				      ior(kind_covPsi,kind_covChi))
#include "assert.H"
contains

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: symCx_ - operator y=P(x), P is Fcst.Err.Cov.
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine symCx_(fecov,comm, nav,krNav,ktNav,	&
	n_x,ktab,jtab,sigPhi,corPhi,corPsi,		&
	x, y						)

      use m_Navigator,only : Navigator
      use m_Navigator,only : lsize
      use m_Navigator,only : get
      use m_sigmaPhi ,only : sigmaPhi
      use m_CorAttrF ,only : CorAttrF
      use m_die,only : die
      implicit none

      type(FcstErrCovMatx),intent(inout) :: fecov
      integer,intent(in) :: comm	! message-passing communicator

      type(Navigator),intent(in) :: nav
      integer,dimension(:),intent(in) :: krNav
      integer,dimension(:),intent(in) :: ktNav

      integer,intent(in) :: n_x
      integer,dimension(:),intent(in) :: ktab
      integer,dimension(:),intent(in) :: jtab
      type(sigmaPhi),intent(in) :: sigPhi
      type(CorAttrF),intent(in) :: corPhi
      type(CorAttrF),intent(in) :: corPsi

      real   ,dimension(:,:),intent(in ) :: x
      real   ,dimension(:,:),intent(out) :: y

! !REVISION HISTORY:
! 	25May00	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::symCx_'

  if(.not.fecov%symmetric) call die(myname_,'fecov is not symmetric')

	! Initialize all vector segments of y to zero

  y(:,:)=0.

	! y = Cx

  call symCxpy_(fecov,comm,		&
	nav,krNav,ktNav,		&
	n_x,ktab,jtab,sigPhi,		&
	corPhi,corPsi,			&
	x, y)

end subroutine symCx_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: symCxpy_ - operator y=P(x)+y, P is Fcst.Err.Cov.
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine symCxpy_(fecov,comm,				&
	nav,krNav,ktNav, n_x, ktab, jtab, sigPhi,corPhi,corPsi,	&
	x, y)

      use m_sigmaPhi,only : sigmaPhi
      use m_sigmaPhi,only : col_size
      use m_sigmaPhi,only : sigmaPhi_matvec
      use m_sigmaPhi,only : sigmaPhi_adjvec
      use m_sigmaPsi,only : sigmaPsi_matvec
      use m_sigmaPsi,only : sigmaPsi_adjvec
      use m_sigmaChi,only : sigmaChi_matvec
      use m_sigmaChi,only : sigmaChi_adjvec
      use m_gammaPhi,only : gammaPhi_matvecpy
      use m_gammaPhi,only : gammaPhi_adjvec
      use m_gammaPsi,only : gammaPsi_matvecpy
      use m_gammaPsi,only : gammaPsi_adjvec
      use m_gammaChi,only : gammaChi_matvecpy
      use m_gammaChi,only : gammaChi_adjvec

      use m_CorAttrF,only : CorAttrF
      use m_CorAttrF,only : ptr_Navigator
      use m_CorAttrF,only : ptr_krNav
      use m_CorAttrF,only : ptr_ktNav
      use m_CorAttrF,only : ptr_qr
      use m_CorAttrF,only : ptr_qd
      use m_CorAttrF,only : ptr_kl
      use m_CorAttrF,only : ptr_collVec

      use m_CorMatxF,only : symCx
      use m_die  ,only : die,assert_
      use m_mall ,only : mall_ison,mall_mci,mall_mco
      use m_Navigator,only : Navigator

      use m_zeit, only : zeit_ci, zeit_co
      implicit none

      type(FcstErrCovMatx),intent(inout) :: fecov
      integer,intent(in) :: comm	! message-passing communicator

      type(Navigator),intent(in) :: nav
      integer,dimension(:),intent(in) :: krNav
      integer,dimension(:),intent(in) :: ktNav

      integer,intent(in) :: n_x
      integer,dimension(:),intent(in) :: ktab
      integer,dimension(:),intent(in) :: jtab
      type(sigmaPhi),intent(in) :: sigPhi
      type(CorAttrF),intent(in) :: corPhi
      type(CorAttrF),intent(in) :: corPsi

      real   ,dimension(:,:),intent(in   ) :: x
      real   ,dimension(:,:),intent(inout) :: y

! !REVISION HISTORY:
! 	25May00	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::symCxpy_'
  integer :: ier
  real,allocatable,dimension(:,:) :: tx,ty
  real,allocatable,dimension(:,:) :: qx,qy
      Interface 
         logical function sparse(kind_mat,kind_cov, kr_i,kt_i, kr_j,kt_j)
           integer, intent(in)	:: kind_mat	! which matrix
           integer, intent(in)	:: kind_cov	! which covariance
           integer, intent(in)	:: kr_i,kt_i	! row block indices
           integer, intent(in)	:: kr_j,kt_j	! column block indices
         End function sparse
      End Interface
  real,pointer,dimension(:,:) :: qr,qd
  integer :: ncor
  integer :: nvecs

  if(.not.fecov%symmetric) call die(myname_,'fecov is not symmetric')

  nvecs=size(x,2)
	ASSERT( nvecs==size(y,2) )

  ncor = col_size(sigPhi)

	allocate( tx(n_x ,nvecs),ty(n_x ,nvecs),	&
		  qx(ncor,nvecs),qy(ncor,nvecs),stat=ier)
		if(ier/=0) call die(myname_,'allocate()',ier)

		if(mall_ison()) then
		  call mall_mci(tx,myname)
		  call mall_mci(ty,myname)
		  call mall_mci(qx,myname)
		  call mall_mci(qy,myname)
		endif

  if(iand(fecov%kind_cov,kind_covPhi)/=0) then

	! tx = gamPhi'(x)

    call gammaPhi_adjvec(nav,krNav,ktNav, n_x,ktab,jtab,nvecs,x,tx)

	! qx = sigPhi'(tx)

    call sigmaPhi_adjvec(sigPhi, nvecs,tx,qx)

	! qy = corPhi(qx)

	Call zeit_ci('fcst Phi')
    call symCx(fecov%mNavPhi,fecov%spGathPhi,fecov%spScatPhi, &
	 kind_covPhi,comm,	&
	ptr_Navigator(corPhi,whole=.true.),	&
	    ptr_krNav(corPhi,whole=.true.),	&
	    ptr_ktNav(corPhi,whole=.true.),	&
	       ptr_qr(corPhi,whole=.true.),	&
	       ptr_qd(corPhi,whole=.true.),	&
	       ptr_kl(corPhi,whole=.true.),	&
	ptr_collvec(corPhi),			&
	  qx(1:ncor,1:nvecs),			&
	  qy(1:ncor,1:nvecs),	stat=ier	)
	Call zeit_co('fcst Phi')

	if(ier/=0) call die(myname_,'symCx(corPhi)',ier)

	! ty = sigPhi(qy)

    call sigmaPhi_matvec(sigPhi, nvecs,qy,ty)

	! y = y + gamPhi(ty)

    call gammaPhi_matvecpy(nav,krNav,ktNav, n_x,ktab,jtab,nvecs,ty,y)

  endif

  if(iand(fecov%kind_cov,kind_covPsi)/=0) then

	! ty = gamPsi'(x)

    call gammaPsi_adjvec(nav,krNav,ktNav, n_x,ktab,jtab,nvecs,x,ty)

	! tx = sigPsi'(ty)

    call sigmaPsi_adjvec(nav,krNav,ktNav, n_x,ktab,jtab,nvecs,ty,tx)

	! ty = corPsi(tx)
	Call zeit_ci('fcst Psi')

    call symCx(fecov%mNavPsi,fecov%spGathPsi,fecov%spScatPsi, &
	 kind_covPsi,comm,	&
	ptr_Navigator(corPsi,whole=.true.),	&
	    ptr_krNav(corPsi,whole=.true.),	&
	    ptr_ktNav(corPsi,whole=.true.),	&
	       ptr_qr(corPsi,whole=.true.),	&
	       ptr_qd(corPsi,whole=.true.),	&
	       ptr_kl(corPsi,whole=.true.),	&
	ptr_collvec(corPsi),			&
	  tx(1:n_x,1:nvecs),			&
	  ty(1:n_x,1:nvecs),	stat=ier	)
	Call zeit_co('fcst Psi')

	if(ier/=0) call die(myname_,'sym_Cxpy(corPsi)',ier)

	! tx = sigPsi(ty)

    call sigmaPsi_matvec(nav,krNav,ktNav, n_x,ktab,jtab,nvecs,ty,tx)

	! y = y + gamPsi(tx)

    call gammaPsi_matvecpy(nav,krNav,ktNav, n_x,ktab,jtab,nvecs,tx,y)

  endif

  if(iand(fecov%kind_cov,kind_covChi)/=0) then

	! ty = gamChi'(x)

    call gammaChi_adjvec(nav,krNav,ktNav, n_x,ktab,jtab,nvecs,x,ty)

	! tx = sigChi'(ty)

    call sigmaChi_adjvec(nav,krNav,ktNav, n_x,ktab,jtab,nvecs,ty,tx)

	! ty = corChi(tx)

	Call zeit_ci('fcst Chi')

    call symCx(fecov%mNavChi,fecov%spGathChi,fecov%spScatChi, &
	 kind_covChi,comm,	&
	ptr_Navigator(corPsi,whole=.true.),	&
	    ptr_krNav(corPsi,whole=.true.),	&
	    ptr_ktNav(corPsi,whole=.true.),	&
	       ptr_qr(corPsi,whole=.true.),	&
	       ptr_qd(corPsi,whole=.true.),	&
	       ptr_kl(corPsi,whole=.true.),	&
	ptr_collvec(corPsi),			&
	  tx(1:n_x,1:nvecs),			&
	  ty(1:n_x,1:nvecs),	stat=ier	)
	Call zeit_co('fcst Chi')

	if(ier/=0) call die(myname_,'sym_Cxpy(corChi)',ier)

	! tx = sigChi(ty)

    call sigmaChi_matvec(nav,krNav,ktNav, n_x,ktab,jtab,nvecs,ty,tx)

	! y = y + gamChi(tx)

    call gammaChi_matvecpy(nav,krNav,ktNav, n_x,ktab,jtab,nvecs,tx,y)

  endif

		if(mall_ison()) then
		  call mall_mco(tx,myname)
		  call mall_mco(ty,myname)
		  call mall_mco(qx,myname)
		  call mall_mco(qy,myname)
		endif

	deallocate(tx,ty,qx,qy,stat=ier)
		if(ier/=0) call die(myname_,'deallocate()',ier)

end subroutine symCxpy_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: recCx_ - operator y=P(x), P is Fcst.Err.Cov.
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine recCx_(fecov,comm,				&
	navi,krNavi,ktNavi, n_i, ktabi,jtabi,			&
		sigPhii,corPhii,corPsii,			&
	navj,krNavj,ktNavj, n_j, ktabj,jtabj,			&
		sigPhij,corPhij,corPsij,			&
	x, y)

      use m_Navigator,only : Navigator
      use m_Navigator,only : lsize
      use m_Navigator,only : get

      use m_sigmaPhi ,only : sigmaPhi
      use m_CorAttrF ,only : CorAttrF

      use m_die,only : die

      implicit none

      type(FcstErrCovMatx),intent(inout) :: fecov
      integer,intent(in) :: comm	! message passing communicator

		! Attributes of vector y

      type(Navigator),intent(in) :: navi
      integer,dimension(:),intent(in) :: krNavi
      integer,dimension(:),intent(in) :: ktNavi

      integer,intent(in) :: n_i
      integer,dimension(:),intent(in) :: ktabi
      integer,dimension(:),intent(in) :: jtabi
      type(sigmaPhi),intent(in) :: sigPhii
      type(CorAttrF),intent(in) :: corPhii
      type(CorAttrF),intent(in) :: corPsii

		! Attributes of vector x

      type(Navigator),intent(in) :: navj
      integer,dimension(:),intent(in) :: krNavj
      integer,dimension(:),intent(in) :: ktNavj

      integer,intent(in) :: n_j
      integer,dimension(:),intent(in) :: ktabj
      integer,dimension(:),intent(in) :: jtabj
      type(sigmaPhi),intent(in) :: sigPhij
      type(CorAttrF),intent(in) :: corPhij
      type(CorAttrF),intent(in) :: corPsij

      real   ,dimension(:,:),intent(in ) :: x	! (n_j,nvecs)
      real   ,dimension(:,:),intent(out) :: y	! (n_i,nvecs)

! !REVISION HISTORY:
! 	25May00	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::recCx_'

  if(fecov%symmetric) call die(myname_,'fecov is symmetric')

	! Initialize vector y to zero

  y(:,:)=0.

	! y = Cx + y

  call recCxpy_(fecov,comm,			&
	navi,krNavi,ktNavi, n_i,ktabi,jtabi,	&
		sigPhii,corPhii,corPsii,	&
	navj,krNavj,ktNavj, n_j,ktabj,jtabj,	&
		sigPhij,corPhij,corPsij,	&
	x, y)

end subroutine recCx_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: recCxpy_ - operator y=P(x)+y, P is Fcst.Err.Cov.
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine recCxpy_(fecov,comm,		&
	navi,krNavi,ktNavi, n_i, ktabi, jtabi,	&
		sigPhii,corPhii,corPsii,	&
	navj,krNavj,ktNavj, n_j, ktabj, jtabj,	&
		sigPhij,corPhij,corPsij,	&
	x, y)

      use m_sigmaPhi,only : sigmaPhi
      use m_sigmaPhi,only : col_size
      use m_sigmaPhi,only : sigmaPhi_matvec
      use m_sigmaPhi,only : sigmaPhi_adjvec
      use m_sigmaPsi,only : sigmaPsi_matvec
      use m_sigmaPsi,only : sigmaPsi_adjvec
      use m_sigmaChi,only : sigmaChi_matvec
      use m_sigmaChi,only : sigmaChi_adjvec
      use m_gammaPhi,only : gammaPhi_matvecpy
      use m_gammaPhi,only : gammaPhi_adjvec
      use m_gammaPsi,only : gammaPsi_matvecpy
      use m_gammaPsi,only : gammaPsi_adjvec
      use m_gammaChi,only : gammaChi_matvecpy
      use m_gammaChi,only : gammaChi_adjvec

      use m_CorAttrF,only : CorAttrF
      use m_CorAttrF,only : ptr_Navigator
      use m_CorAttrF,only : ptr_krNav
      use m_CorAttrF,only : ptr_ktNav
      use m_CorAttrF,only : ptr_qr
      use m_CorAttrF,only : ptr_qd
      use m_CorAttrF,only : ptr_kl
      use m_CorAttrF,only : ptr_collVec

      use m_CorMatxF,only : recCx
      use m_die  ,only : die,assert_
      use m_mall ,only : mall_ison,mall_mci,mall_mco
      use m_Navigator,only : Navigator
      !use m_parDOT,only : parDOT

      implicit none

      type(FcstErrCovMatx),intent(inout) :: fecov
      integer,intent(in) :: comm	! message-passing communicator

		! attributes of vector y

      type(Navigator),intent(in) :: navi
      integer,dimension(:),intent(in) :: krNavi
      integer,dimension(:),intent(in) :: ktNavi

      integer,intent(in) :: n_i
      integer,dimension(:),intent(in) :: ktabi
      integer,dimension(:),intent(in) :: jtabi
      type(sigmaPhi),intent(in) :: sigPhii
      type(CorAttrF),intent(in) :: corPhii
      type(CorAttrF),intent(in) :: corPsii

		! attributes of vector x

      type(Navigator),intent(in) :: navj
      integer,dimension(:),intent(in) :: krNavj
      integer,dimension(:),intent(in) :: ktNavj

      integer,intent(in) :: n_j
      integer,dimension(:),intent(in) :: ktabj
      integer,dimension(:),intent(in) :: jtabj
      type(sigmaPhi),intent(in) :: sigPhij
      type(CorAttrF),intent(in) :: corPhij
      type(CorAttrF),intent(in) :: corPsij

      real   ,dimension(:,:),intent(in   ) :: x	! (n_j,nvecs)
      real   ,dimension(:,:),intent(inout) :: y	! (n_i,nvecs)

! !REVISION HISTORY:
! 	25May00	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::recCxpy_'
  integer :: ier
  real,allocatable,dimension(:,:) :: tx,ty
  real,allocatable,dimension(:,:) :: px,py
  real,allocatable,dimension(:,:) :: qx,qy
      Interface 
         logical function sparse(kind_mat,kind_cov, kr_i,kt_i, kr_j,kt_j)
           integer, intent(in)	:: kind_mat	! which matrix
           integer, intent(in)	:: kind_cov	! which covariance
           integer, intent(in)	:: kr_i,kt_i	! row block indices
           integer, intent(in)	:: kr_j,kt_j	! column block indices
         End function sparse
      End Interface
  real,pointer,dimension(:,:) :: qri,qdi
  real,pointer,dimension(:,:) :: qrj,qdj
  integer :: ncori
  integer :: ncorj
  integer :: nvecs

  if(fecov%symmetric) call die(myname_,'fecov is symmetric')

  nvecs=size(x,2)
	ASSERT( nvecs == size(y,2) )

  ncori=col_size(sigPhii)
  ncorj=col_size(sigPhij)

	allocate( tx(n_j  ,nvecs),ty(n_i  ,nvecs),	&
		  px(n_j  ,nvecs),py(n_i  ,nvecs),	&
		  qx(ncorj,nvecs),qy(ncori,nvecs),	stat=ier)
		if(ier/=0) call die(myname_,'allocate()',ier)

		if(mall_ison()) then
		  call mall_mci(tx,myname)
		  call mall_mci(ty,myname)
		  call mall_mci(px,myname)
		  call mall_mci(py,myname)
		  call mall_mci(qx,myname)
		  call mall_mci(qy,myname)
		endif

  if(iand(fecov%kind_cov,kind_covPhi)/=0) then
    !print*,myname_,'::dot( x1)=',parDOT( x(1:n_j,1), x(1:n_j,1),comm)

	! tx = gamPhi'(x)

    call gammaPhi_adjvec(navj,krNavj,ktNavj,	&
	n_j,ktabj,jtabj,nvecs,x,tx)
    !print*,myname_,'::dot(tx1)=',parDOT(tx(1:n_j,1),tx(1:n_j,1),comm)

	! qx = sigPhi'(tx)

    call sigmaPhi_adjvec(sigPhij, nvecs,tx,qx)
    !print*,myname_,'::dot(qx1)=',parDOT(qx(1:ncorj,1),qx(1:ncorj,1),comm)

	! qy = corPhi(qx)

    call recCx(fecov%mNavPhi,fecov%spGathPhi,fecov%spScatPhi, &
	kind_covPhi,comm,			&
	ptr_Navigator(corPhii,whole=.true.),	&
	    ptr_krNav(corPhii,whole=.true.),	&
	    ptr_ktNav(corPhii,whole=.true.),	&
	       ptr_qr(corPhii,whole=.true.),	&
	       ptr_qd(corPhii,whole=.true.),	&
	       ptr_kl(corPhii,whole=.true.),	&
	ptr_Navigator(corPhij,whole=.true.),	&
	    ptr_krNav(corPhij,whole=.true.),	&
	    ptr_ktNav(corPhij,whole=.true.),	&
	       ptr_qr(corPhij,whole=.true.),	&
	       ptr_qd(corPhij,whole=.true.),	&
	       ptr_kl(corPhij,whole=.true.),	&
	ptr_collvec(corPhij),			&
	  qx(1:ncorj,1:nvecs),			&
	ptr_collvec(corPhii),			&
	  qy(1:ncori,1:nvecs),	stat=ier	)

	if(ier/=0) call die(myname_,'rec_Cxpy(corPhi)',ier)
    !print*,myname_,'::dot(qy1)=',parDOT(qy(1:ncori,1),qy(1:ncori,1),comm)

	! ty = sigPhi(qy)

    call sigmaPhi_matvec(sigPhii, nvecs,qy,ty)
    !print*,myname_,'::dot(ty1)=',parDOT(ty(1:n_i,1),ty(1:n_i,1),comm)

	! y = y + gamPhi(ty)

    call gammaPhi_matvecpy(navi,krNavi,ktNavi,	&
	n_i,ktabi,jtabi,nvecs,ty,y)
    !print*,myname_,'::dot( y1)=',parDOT( y(1:n_i,1), y(1:n_i,1),comm)

  endif

  if(iand(fecov%kind_cov,kind_covPsi)/=0) then
    !print*,myname_,'::dot( x2)=',parDOT( x(1:n_j,1), x(1:n_j,1),comm)

	! px = gamPsi'(x)

    call gammaPsi_adjvec(navj,krNavj,ktNavj,	&
	n_j,ktabj,jtabj,nvecs,x,px)
    !print*,myname_,'::dot(px2)=',parDOT(px(1:n_j,1),px(1:n_j,1),comm)

	! tx = sigPsi'(px)

    call sigmaPsi_adjvec(navj,krNavj,ktNavj,	&
	n_j,ktabj,jtabj,nvecs,px,tx)
    !print*,myname_,'::dot(tx2)=',parDOT(tx(1:n_j,1),tx(1:n_j,1),comm)

	! ty = corPsi(tx)

    call recCx(fecov%mNavPsi,fecov%spGathPsi,fecov%spScatPsi,&
	kind_covPsi,comm,			&
	ptr_Navigator(corPsii,whole=.true.),	&
	    ptr_krNav(corPsii,whole=.true.),	&
	    ptr_ktNav(corPsii,whole=.true.),	&
	       ptr_qr(corPsii,whole=.true.),	&
	       ptr_qd(corPsii,whole=.true.),	&
	       ptr_kl(corPsii,whole=.true.),	&
	ptr_Navigator(corPsij,whole=.true.),	&
	    ptr_krNav(corPsij,whole=.true.),	&
	    ptr_ktNav(corPsij,whole=.true.),	&
	       ptr_qr(corPsij,whole=.true.),	&
	       ptr_qd(corPsij,whole=.true.),	&
	       ptr_kl(corPsij,whole=.true.),	&
	ptr_collvec(corPsij),			&
	  tx(1:n_j,1:nvecs),			&
	ptr_collvec(corPsii),			&
	  ty(1:n_i,1:nvecs),	stat=ier	)

	if(ier/=0) call die(myname_,'rec_Cxpy(corPsi)',ier)
    !print*,myname_,'::dot(ty2)=',parDOT(ty(1:n_i,1),ty(1:n_i,1),comm)

	! py = sigPsi(ty)

    call sigmaPsi_matvec(navi,krNavi,ktNavi,	&
	n_i,ktabi,jtabi,nvecs,ty,py)
    !print*,myname_,'::dot(py2)=',parDOT(py(1:n_i,1),py(1:n_i,1),comm)

	! y = y + gamPsi(py)

    call gammaPsi_matvecpy(navi,krNavi,ktNavi,	&
	n_i,ktabi,jtabi,nvecs,py,y)
    !print*,myname_,'::dot( y2)=',parDOT( y(1:n_i,1), y(1:n_i,1),comm)

  endif

  if(iand(fecov%kind_cov,kind_covChi)/=0) then
    !print*,myname_,'::dot( x3)=',parDOT( x(1:n_j,1), x(1:n_j,1),comm)

	! px = gamChi'(x)

    call gammaChi_adjvec(navj,krNavj,ktNavj,	&
	n_j,ktabj,jtabj,nvecs,x,px)
    !print*,myname_,'::dot(px3)=',parDOT(px(1:n_j,1),px(1:n_j,1),comm)

	! tx = sigChi'(px)

    call sigmaChi_adjvec(navj,krNavj,ktNavj,	&
	n_j,ktabj,jtabj,nvecs,px,tx)
    !print*,myname_,'::dot(tx3)=',parDOT(tx(1:n_j,1),tx(1:n_j,1),comm)

	! ty = corChi(tx)

    call recCx(fecov%mNavChi,fecov%spGathChi,fecov%spScatChi,&
	kind_covChi,comm,			&
	ptr_Navigator(corPsii,whole=.true.),	&
	    ptr_krNav(corPsii,whole=.true.),	&
	    ptr_ktNav(corPsii,whole=.true.),	&
	       ptr_qr(corPsii,whole=.true.),	&
	       ptr_qd(corPsii,whole=.true.),	&
	       ptr_kl(corPsii,whole=.true.),	&
	ptr_Navigator(corPsij,whole=.true.),	&
	    ptr_krNav(corPsij,whole=.true.),	&
	    ptr_ktNav(corPsij,whole=.true.),	&
	       ptr_qr(corPsij,whole=.true.),	&
	       ptr_qd(corPsij,whole=.true.),	&
	       ptr_kl(corPsij,whole=.true.),	&
	ptr_collvec(corPsij),			&
	  tx(1:n_j,1:nvecs),			&
	ptr_collvec(corPsii),			&
	  ty(1:n_i,1:nvecs),	stat=ier	)

	if(ier/=0) call die(myname_,'rec_Cxpy(corChi)',ier)
    !print*,myname_,'::dot(ty3)=',parDOT(ty(1:n_i,1),ty(1:n_i,1),comm)

	! py = sigChi(ty)

    call sigmaChi_matvec(navi,krNavi,ktNavi,	&
	n_i,ktabi,jtabi,nvecs,ty,py)
    !print*,myname_,'::dot(py3)=',parDOT(py(1:n_i,1),py(1:n_i,1),comm)

	! y = y + gamChi(py)

    call gammaChi_matvecpy(navi,krNavi,ktNavi,	&
	n_i,ktabi,jtabi,nvecs,py,y)
    !print*,myname_,'::dot( y3)=',parDOT( y(1:n_i,1), y(1:n_i,1),comm)

  endif

		if(mall_ison()) then
		  call mall_mco(tx,myname)
		  call mall_mco(px,myname)
		  call mall_mco(qx,myname)
		  call mall_mco(ty,myname)
		  call mall_mco(py,myname)
		  call mall_mco(qy,myname)
		endif
	deallocate(tx,px,qx,ty,py,qy,stat=ier)
		if(ier/=0) call die(myname_,'deallocate()',ier)

end subroutine recCxpy_


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: symSched_ - initialize an object
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine symSched_(fecov,corPhi,corPsi,comm,kind_mat,kind_cov)
      use m_symMatx ,only : symMatx_schedule
      use m_CorAttrF,only : CorAttrF
      use m_CorAttrF,only : ptr_Navigator
      use m_CorAttrF,only : ptr_krNav
      use m_CorAttrF,only : ptr_ktNav
      use m_CorAttrF,only : ptr_collNav
      use m_Collector,only : Collector
      use m_sparse, only : sparse
      implicit none

      type(FcstErrCovMatx),intent(out) :: fecov
      type(CorAttrF),intent(in) :: corPhi
      type(CorAttrF),intent(in) :: corPsi
      integer,intent(in) :: comm
      integer,intent(in) :: kind_mat
      integer,optional,intent(in) :: kind_cov

! !REVISION HISTORY:
! 	25May00	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::symSched_'

  Type(Collector), Pointer :: collnav

  fecov%kind_cov = kind_FECov
  if(present(kind_cov)) fecov%kind_cov=iand(kind_cov,kind_FECov)

 
  if(iand(fecov%kind_cov,kind_covPhi) /=0) then
     collnav => ptr_collNav(corPhi)

    call symMatx_schedule(fecov%mNavPhi,            &
	 fecov%spGathPhi, fecov%spScatPhi, collnav, &
	ptr_Navigator(corPhi,whole=.true.),	    &
	    ptr_krNav(corPhi,whole=.true.),	    &
	    ptr_ktNav(corPhi,whole=.true.),	    &
	comm,kind_covPhi,kind_mat,sparse,myCost=fecov%myCostPhi)
  endif

  if(iand(fecov%kind_cov,kind_covPsi) /=0) then
     collnav => ptr_collNav(corPsi)

    call symMatx_schedule(fecov%mNavPsi,	    &
	 fecov%spGathPsi, fecov%spScatPsi, collnav, &
	ptr_Navigator(corPsi,whole=.true.),	    &
	    ptr_krNav(corPsi,whole=.true.),	    &
	    ptr_ktNav(corPsi,whole=.true.),	    &
	comm,kind_covPsi,kind_mat,sparse,myCost=fecov%myCostPsi)
  endif

  if(iand(fecov%kind_cov,kind_covChi) /=0) then
     collnav => ptr_collNav(corPsi)

    call symMatx_schedule(fecov%mNavChi,	    &
	 fecov%spGathChi, fecov%spScatChi, collnav, &
	ptr_Navigator(corPsi,whole=.true.),	    &
	    ptr_krNav(corPsi,whole=.true.),	    &
	    ptr_ktNav(corPsi,whole=.true.),	    &
	comm,kind_covChi,kind_mat,sparse,myCost=fecov%myCostChi)
  endif

  Nullify(collnav)

  fecov%symmetric=.true.

end subroutine symSched_
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: recSched_ - initialize an object
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine recSched_(fecov,row_corPhi,row_corPsi,	&
	col_corPhi,col_corPsi,comm,kind_mat,kind_cov)

      use m_recMatx ,only : recMatx_schedule

      use m_CorAttrF,only : CorAttrF
      use m_CorAttrF,only : ptr_Navigator
      use m_CorAttrF,only : ptr_krNav
      use m_CorAttrF,only : ptr_ktNav
      use m_CorAttrF,only : ptr_collNav
      use m_sparse,  only : sparse

      implicit none

      type(FcstErrCovMatx),intent(out) :: fecov
      type(CorAttrF),intent(in) :: row_corPhi
      type(CorAttrF),intent(in) :: row_corPsi
      type(CorAttrF),intent(in) :: col_corPhi
      type(CorAttrF),intent(in) :: col_corPsi
      integer,intent(in) :: comm
      integer,intent(in) :: kind_mat
      integer,optional,intent(in) :: kind_cov

! !REVISION HISTORY:
! 	25May00	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::recSched_'

  fecov%kind_cov = kind_FECov
  if(present(kind_cov)) fecov%kind_cov=iand(kind_cov,kind_FECov)

  if(iand(fecov%kind_cov,kind_covPhi) /=0) then

    call recMatx_schedule(fecov%mNavPhi,fecov%spGathPhi,fecov%spScatPhi,	&
	ptr_collNav(row_corPhi), ptr_collNav(col_corPhi),  	&
	ptr_Navigator(row_corPhi,whole=.true.),	&
	    ptr_krNav(row_corPhi,whole=.true.),	&
	    ptr_ktNav(row_corPhi,whole=.true.),	&
	ptr_Navigator(col_corPhi,whole=.true.),	&
	    ptr_krNav(col_corPhi,whole=.true.),	&
	    ptr_ktNav(col_corPhi,whole=.true.),	&
	comm,kind_covPhi,kind_mat,sparse,myCost=fecov%myCostPhi)
  endif

  if(iand(fecov%kind_cov,kind_covPsi) /=0) then

    call recMatx_schedule(fecov%mNavPsi,fecov%spGathPsi,fecov%spScatPsi,	&
	ptr_collNav(row_corPsi), ptr_collNav(col_corPsi), 	&
	ptr_Navigator(row_corPsi,whole=.true.),	&
	    ptr_krNav(row_corPsi,whole=.true.),	&
	    ptr_ktNav(row_corPsi,whole=.true.),	&
	ptr_Navigator(col_corPsi,whole=.true.),	&
	    ptr_krNav(col_corPsi,whole=.true.),	&
	    ptr_ktNav(col_corPsi,whole=.true.),	&
	comm,kind_covPsi,kind_mat,sparse,myCost=fecov%myCostPsi)
  endif

  if(iand(fecov%kind_cov,kind_covChi) /=0) then

    call recMatx_schedule(fecov%mNavChi,fecov%spGathChi,fecov%spScatChi,	&
	ptr_collNav(row_corPsi), ptr_collNav(col_corPsi), &
	ptr_Navigator(row_corPsi,whole=.true.),	&
	    ptr_krNav(row_corPsi,whole=.true.),	&
	    ptr_ktNav(row_corPsi,whole=.true.),	&
	ptr_Navigator(col_corPsi,whole=.true.),	&
	    ptr_krNav(col_corPsi,whole=.true.),	&
	    ptr_ktNav(col_corPsi,whole=.true.),	&
	comm,kind_covChi,kind_mat,sparse,myCost=fecov%myCostChi)
  endif

  fecov%symmetric=.false.

end subroutine recSched_
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: clean_ - clean an object
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine clean_(fecov)
      use m_AttrVect,only : clean
      use m_SparseComm,only : clean_spComm => clean
      implicit none
      type(FcstErrCovMatx),intent(inout) :: fecov

! !REVISION HISTORY:
! 	25May00	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::clean_'

  if(iand(fecov%kind_cov,kind_covPhi)/=0) Then
	call clean(fecov%mNavPhi)
	call clean_spComm(fecov%spGathPhi)
	call clean_spComm(fecov%spScatPhi)
  end if
  if(iand(fecov%kind_cov,kind_covPsi)/=0) then
	call clean(fecov%mNavPsi)
	call clean_spComm(fecov%spGathPsi)
	call clean_spComm(fecov%spScatPsi)
  end if
  if(iand(fecov%kind_cov,kind_covChi)/=0) then
	call clean(fecov%mNavChi)
	call clean_spComm(fecov%spGathChi)
	call clean_spComm(fecov%spScatChi)
  end if

  fecov%myCostPhi=0.
  fecov%myCostPsi=0.
  fecov%myCostChi=0.

  fecov%kind_cov=0
  fecov%symmetric=.false.

end subroutine clean_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: showCosts_ - show myCost values of all matrix operators
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine showCosts_(fecov,lu,root,comm,listall)
      use m_showDistrib,only : showDistrib
      implicit none
      type(FcstErrCovMatx),intent(in) :: fecov
      integer,intent(in) :: lu		! only significant on root
      integer,intent(in) :: root	! root PE ID
      integer,intent(in) :: comm	! communicator
      logical,optional,intent(in) :: listall

! !REVISION HISTORY:
! 	05Mar01	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::showCosts_'

  call showDistrib(lu,myname//'::Phi',fecov%myCostPhi,	&
	root,comm,listall=listall)
  call showDistrib(lu,myname//'::Psi',fecov%myCostPsi,	&
	root,comm,listall=listall)
  call showDistrib(lu,myname//'::Chi',fecov%myCostChi,	&
	root,comm,listall=listall)

end subroutine showCosts_
end module m_FcstErrCovMatx
