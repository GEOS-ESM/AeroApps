!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !MODULE: m_InnovCovMatx - an initialized innov. covariance operator
!
! !DESCRIPTION:
!
!     Note that bue to the symmetry of a covariance matrix operator,
!   the adjoint of the operator is abtained simply by switching the
!   row arguments and the column arguments of this subroutine.
!
! !INTERFACE:

    module m_InnovCovMatx
      use m_FcstErrCovMatx,only : FcstErrCovMatx
      use m_FcstErrCovMatx,only : kind_FECov
      use m_ObsErrCovMatx ,only :  ObsErrCovMatx
      use m_ObsErrCovMatx ,only : kind_OECov
      implicit none
      private	! except

      public :: InnovCovMatx		! The class data structure
      public :: InnovCovMatx_Cxpy
      public :: InnovCovMatx_Cx
      public :: InnovCovMatx_init
      public :: InnovCovMatx_clean
      public :: InnovCovMatx_showCosts
      public :: kind_INCov

    type InnovCovMatx
      private
      type(FcstErrCovMatx) :: fecov
      type( ObsErrCovMatx) :: oecov
    end type InnovCovMatx

      interface InnovCovMatx_Cxpy; module procedure	&
	symCxpy_
      end interface

      interface InnovCovMatx_Cx; module procedure	&
	symCx_
      end interface

      interface InnovCovMatx_init; module procedure	&
	symSched_; end interface
      interface InnovCovMatx_clean; module procedure	&
	clean_; end interface
      interface InnovCovMatx_showCosts; module procedure	&
	showCosts_; end interface

! !REVISION HISTORY:
! 	25May00	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname='m_InnovCovMatx'

  integer,parameter :: kind_INCov = ior(kind_FECov,kind_OECov)

contains

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: symCx_ - operator y=P(x), P is Innov.Err.Cov.
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine symCx_(incov,comm, nav,krNav,ktNav,		&
	n_x,sigOu,sigOc,corObs,ktab,jtab,sigPhi,corPhi,corPsi,	&
	x, y)

      use m_Navigator,only : Navigator
      use m_sigmaPhi ,only : sigmaPhi
      use m_CorAttrF ,only : CorAttrF
      use m_CorAttrX ,only : CorAttrX
      implicit none

      type(InnovCovMatx),intent(inout) :: incov
      integer,intent(in) :: comm	! message-passing communicator

      type(Navigator),intent(in) :: nav
      integer,dimension(:),intent(in) :: krNav
      integer,dimension(:),intent(in) :: ktNav

      integer,intent(in) :: n_x
      real   ,dimension(:),intent(in) :: sigOu
      real   ,dimension(:),intent(in) :: sigOc
      type(CorAttrX) ,intent(in) :: corObs

      integer,dimension(:),intent(in) :: ktab
      integer,dimension(:),intent(in) :: jtab

      type(sigmaPhi),intent(in) :: sigPhi
      type(CorAttrF) ,intent(in) :: corPhi
      type(CorAttrF) ,intent(in) :: corPsi

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

  call symCxpy_(incov,comm, nav,krNav,ktNav,	&
	n_x,sigOu,sigOc,corObs,			&
	ktab,jtab,sigPhi,corPhi,corPsi,		&
	x, y)

end subroutine symCx_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: symCxpy_ - operator y=P(x)+y, P is Innov.Err.Cov.
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine symCxpy_(incov,comm,		&
	nav,krNav,ktNav,			&
	n_x,sigOu,sigOc,corObs,			&
	ktab,jtab, sigPhi,corPhi,corPsi,	&
	x, y)

      use m_ObsErrCovMatx ,only :  ObsErrCovMatx_Cxpy
      use m_FcstErrCovMatx,only : FcstErrCovMatx_Cxpy

      use m_die  ,only : die
      use m_Navigator,only : Navigator
      use m_sigmaPhi ,only : sigmaPhi
      use m_CorAttrF ,only : CorAttrF
      use m_CorAttrX ,only : CorAttrX
      implicit none

      type(InnovCovMatx),intent(inout) :: incov
      integer,intent(in) :: comm	! message-passing communicator

      type(Navigator),intent(in) :: nav
      integer,dimension(:),intent(in) :: krNav
      integer,dimension(:),intent(in) :: ktNav

      integer,intent(in) :: n_x
      real   ,dimension(:),intent(in) :: sigOu
      real   ,dimension(:),intent(in) :: sigOc
      type(CorAttrX),intent(in) :: corObs

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

  call ObsErrCovMatx_Cxpy(incov%oecov,comm,	&
	nav,krNav,ktNav, sigOu,sigOc,corObs,	&
	x, y)

  call FcstErrCovMatx_Cxpy(incov%fecov,comm,	&
	nav,krNav,ktNav,			&
	n_x,ktab,jtab,sigPhi,corPhi,corPsi,	&
	x, y)

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

    subroutine symSched_(incov,corObs,corPhi,corPsi,	&
	comm,kind_mat,kind_cov)
      use m_CorAttrF ,only : CorAttrF
      use m_CorAttrX ,only : CorAttrX
      use m_FcstErrCovMatx,only : FcstErrCovMatx_init
      use m_ObsErrCovMatx ,only :  ObsErrCovMatx_init
      implicit none
      type(InnovCovMatx),intent(out) :: incov
      type(CorAttrX),intent(in) :: corObs
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

  call  ObsErrCovMatx_init(incov%oecov,corObs,		&
	comm,kind_mat,kind_cov=kind_cov)

  call FcstErrCovMatx_init(incov%fecov,corPhi,corPsi,	&
	comm,kind_mat,kind_cov=kind_cov)

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

    subroutine clean_(incov)
      use m_FcstErrCovMatx,only : FcstErrCovMatx_clean
      use m_ObsErrCovMatx ,only : ObsErrCovMatx_clean
      implicit none
      type(InnovCovMatx),intent(inout) :: incov

! !REVISION HISTORY:
! 	25May00	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::clean_'

  call  ObsErrCovMatx_clean(incov%oecov)
  call FcstErrCovMatx_clean(incov%fecov)

end subroutine clean_
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: showCosts_ - show costs of all operators
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine showCosts_(incov,lu,root,comm,listall)
      use m_FcstErrCovMatx,only : FcstErrCovMatx_showCosts
      use m_ObsErrCovMatx ,only : ObsErrCovMatx_showCosts
      implicit none
      type(InnovCovMatx),intent(in) :: incov	! a matrix operator
      integer,intent(in) :: lu		! unit significant on root
      integer,intent(in) :: root	! root PE ID
      integer,intent(in) :: comm	! communicator
      logical,optional,intent(in) :: listall	! if list all values

! !REVISION HISTORY:
! 	05Mar01	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::showCosts_'

  call FcstErrCovMatx_showCosts(incov%fecov,lu,	&
	root,comm,listall=listall)
  call ObsErrCovMatx_showCosts(incov%oecov,lu,	&
	root,comm,listall=listall)

end subroutine showCosts_
end module m_InnovCovMatx
