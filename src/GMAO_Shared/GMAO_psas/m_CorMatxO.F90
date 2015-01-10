!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !MODULE: m_CorMatxO - blocked matrix attributes
!
! !DESCRIPTION:
!
! !INTERFACE:

    module m_CorMatxO
      implicit none
      private	! except

      public :: symCxpy		! y=(symMatx)x+y
      public :: symCx		! y=(symMatx)x
      interface symCxpy; module procedure symCxpy_; end interface
      interface symCx  ; module procedure symCx_  ; end interface

      public :: recCxpy		! y=(recMatx)x+y
      public :: recCx		! y=(recMatx)x
      interface recCxpy; module procedure recCxpy_; end interface
      interface recCx  ; module procedure recCx_  ; end interface

! !REVISION HISTORY:
!	23Aug00	- Jing Guo
!		. Revised from m_symMatx to m_CorMatxF
!        6Oct99 - J.W. Larson <jlarson@dao.gsfc.nasa.gov> - replaced 
!                 occurances of MP_REAL with MP_type(message) for 
!                 better support of 32 and 64-bit platforms.  This 
!                 change affects only the routine Cxpy_().
! 	15Dec98 - Jing Guo <guo@thunder> - initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname='m_CorMatxO'

#include "assert.H"
contains
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: symCx_ - an interface of message-passing symmetric Cx
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine symCx_(mNav,spGath,spScat,kind_C,comm,			&
	nav,krNav,ktNav, kx,qr,kl,collVec,xi,yi, stat	)

      use m_AttrVect ,only : AttrVect
      use m_SparseComm,only: SparseComm
      use m_Collector,only : Collector
      use m_Collector,only : localSize
      use m_Collector,only : globalSize

      use m_Navigator,only : Navigator
      use m_Navigator,only : lsize
      use m_Checks,only : checknorm

      use m_die,only : perr,die,assert_
      use m_mall,only : mall_ison,mall_mci,mall_mco

      implicit none
				! About the matrix
      type(AttrVect),intent(inout) :: mNav
      integer,       intent(in)    :: kind_C
      integer,       intent(in)    :: comm

				! Attributes: parameters for the matrix
      type(Navigator),intent(in) :: nav	
      type(SparseComm),intent(in) :: spGath,spScat
      integer,dimension(:),intent(in) :: krNav		! (nvseg)
      integer,dimension(:),intent(in) :: ktNav		! (nvseg)
      integer,dimension(  :),intent(in) :: kx		! (  nvp)
      real   ,dimension(:,:),intent(in) :: qr		! (3,nvp)
      integer,dimension(  :),intent(in) :: kl		! (  nvp)

				! About the local vectors
      type(Collector),intent(in) :: collVec
      real,dimension(:,:),intent(in ) :: xi		! (nvi,nvecs)
      real,dimension(:,:),intent(out) :: yi		! (nvi,nvecs)

      integer,optional,intent(out) :: stat

! !REVISION HISTORY:
! 	23Aug00	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::symCx_'
  integer :: nvi
  integer :: nvp
  integer :: nvseg
  integer :: nvecs
  integer :: ier
  real,allocatable,dimension(:,:) :: ty

	! Argument consistency verifications

  nvi= localSize(collVec)
	ASSERT( nvi<=size(xi,1) )
	ASSERT( nvi==size(yi,1) )	! yi is an output array

  nvp=globalSize(collVec)
	ASSERT( nvp<=size(kx,1) )
	ASSERT( nvp<=size(qr,2) )
	ASSERT( nvp<=size(kl,1) )

  nvseg=lsize(nav)
	ASSERT( nvseg<=size(krNav) )
	ASSERT( nvseg<=size(ktNav) )

  nvecs=size(xi,2)
	ASSERT( nvecs==size(yi,2) )	! yi is an output array

  if(present(stat)) stat=0

		! define the workspace

	allocate(ty(nvecs,nvi),stat=ier)
		if(ier/=0) then
		  call perr(myname_,'allocate()',ier)
		  if(.not.present(stat)) call die(myname_)
		  stat=ier
		  return
		endif
		if(mall_ison()) call mall_mci(ty,myname)

#ifdef VERIFY_
  call checknorm(myname_//'(xi)',xi(1:nvi,1:nvecs))
  call checknorm(myname_//'(xi)',xi(1:nvi,1:nvecs),comm)
#endif

  call symCx_tr_(mNav,spGath,spScat,kind_C,comm,	&
	nav,krNav,ktNav,kx,qr,kl,collVec,transpose(xi),ty,stat=ier)

	if(ier/=0) then
	  call perr(myname_,'symCx_tr_()',ier)
	  if(.not.present(stat)) call die(myname_)
	  stat=ier
	  return
	endif

  yi(1:nvi,1:nvecs) = transpose(ty)

#ifdef VERIFY_
  call checknorm(myname_//'(yi)',yi(1:nvi,1:nvecs))
  call checknorm(myname_//'(yi)',yi(1:nvi,1:nvecs),comm)
#endif
		if(mall_ison()) call mall_mco(ty,myname)

	deallocate(ty,stat=ier)
		if(ier/=0) then
		  call perr(myname_,'deallocate()',ier)
		  if(.not.present(stat)) call die(myname_)
		  stat=ier
		  return
		endif

end subroutine symCx_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: symCxpy_ - an interface of message-passing symmetric Cx+y
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine symCxpy_(mNav,spGath,spScat,kind_C,comm,		&
	nav,krNav,ktNav, kx,qr,kl,collVec,xi,yi, stat	)

      use m_AttrVect ,only : AttrVect
      use m_SparseComm, only : SparseComm

      use m_Collector,only : Collector
      use m_Collector,only : localSize
      use m_Collector,only : globalSize

      use m_Navigator,only : Navigator
      use m_Navigator,only : lsize
      use m_Checks,only : checknorm

      use m_die,only : perr,die,assert_
      use m_mall,only : mall_ison,mall_mci,mall_mco

      implicit none
				! About the matrix
      type(AttrVect),intent(inout) :: mNav
      integer,       intent(in)    :: kind_C
      integer,       intent(in)    :: comm

				! Attributes: parameters for the matrix
      type(Navigator),intent(in) :: nav	
      type(SparseComm),intent(in)  :: spGath, spScat
      integer,dimension(:),intent(in) :: krNav		! (nvseg)
      integer,dimension(:),intent(in) :: ktNav		! (nvseg)
      integer,dimension(  :),intent(in) :: kx		! (  nvp)
      real   ,dimension(:,:),intent(in) :: qr		! (3,nvp)
      integer,dimension(  :),intent(in) :: kl		! (  nvp)

				! About the local vectors
      type(Collector),intent(in) :: collVec
      real,dimension(:,:),intent(in   ) :: xi		! (nvi,nvecs)
      real,dimension(:,:),intent(inout) :: yi		! (nvi,nvecs)

      integer,optional,intent(out) :: stat

! !REVISION HISTORY:
! 	23Aug00	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::symCxpy_'
  integer :: nvi
  integer :: nvp
  integer :: nvseg
  integer :: nvecs
  integer :: ier
  real,allocatable,dimension(:,:) :: ty

	! Argument consistency verifications

  nvi= localSize(collVec)
	ASSERT( nvi<=size(xi,1) )
	ASSERT( nvi==size(yi,1) )	! yi is an output array

  nvp=globalSize(collVec)
	ASSERT( nvp<=size(kx,1) )
	ASSERT( nvp<=size(qr,2) )
	ASSERT( nvp<=size(kl,1) )

  nvseg=lsize(nav)
	ASSERT( nvseg<=size(krNav) )
	ASSERT( nvseg<=size(ktNav) )

  nvecs=size(xi,2)
	ASSERT( nvecs==size(yi,2) )	! yi is an output array

  if(present(stat)) stat=0

		! define the workspace

	allocate(ty(nvecs,nvi),stat=ier)
		if(ier/=0) then
		  call perr(myname_,'allocate()',ier)
		  if(.not.present(stat)) call die(myname_)
		  stat=ier
		  return
		endif
		if(mall_ison()) call mall_mci(ty,myname)

#ifdef VERIFY_
  call checknorm(myname_//'(xi)',xi(1:nvi,1:nvecs))
  call checknorm(myname_//'(xi)',xi(1:nvi,1:nvecs),comm)
  call checknorm(myname_//'(yi)',yi(1:nvi,1:nvecs))
  call checknorm(myname_//'(yi)',yi(1:nvi,1:nvecs),comm)
#endif

  call symCx_tr_(mNav,spGath,spScat,kind_C,comm,	&
	nav,krNav,ktNav,kx,qr,kl,collVec,transpose(xi),ty,stat=ier)

	if(ier/=0) then
	  call perr(myname_,'symCx_tr_()',ier)
	  if(.not.present(stat)) call die(myname_)
	  stat=ier
	  return
	endif

  yi(1:nvi,1:nvecs) = transpose(ty) + yi(1:nvi,1:nvecs)

#ifdef VERIFY_
  call checknorm(myname_//'(yi)',yi(1:nvi,1:nvecs))
  call checknorm(myname_//'(yi)',yi(1:nvi,1:nvecs),comm)
#endif
		if(mall_ison()) call mall_mco(ty,myname)

	deallocate(ty,stat=ier)
		if(ier/=0) then
		  call perr(myname_,'deallocate()',ier)
		  if(.not.present(stat)) call die(myname_)
		  stat=ier
		  return
		endif

end subroutine symCxpy_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !ROUTINE: symCx_tr_ - a message-passing symmetric Cx kernal
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine symCx_tr_(mNav,spGath,spScat,kind_C,comm,	&
	nav,krNav,ktNav, kx,qr,kl,		&
	collVec,xi,yi	,stat			)

      use m_AttrVect ,only : AttrVect
      use m_SparseComm, only: SparseComm
      use m_SparseComm, only: PartialGather
      use m_SparseComm, only: PartialReduceScatter

      use m_Collector,only : Collector
      use m_Collector,only : localSize
      use m_Collector,only : globalSize
      use m_Collector,only : ptr_counts
      use m_CollectorComm,only : allgatherv

      use m_Navigator,only : Navigator
      use m_Navigator,only : lsize
      use m_Checks,only : checknorm

      use m_CorOxpy,   only : symCxpy

      use m_die   ,only : perr,die
      use m_mall  ,only : mall_ison,mall_mci,mall_mco
      use m_mpif90,only : MP_perr,MP_type,MP_SUM
      use m_zeit  ,only : zeit_ci,zeit_co

      use m_mpout, only : mpout_log, mpout
      implicit none
				! About the matrix
      type(AttrVect),intent(inout) :: mNav
      type(SparseComm),intent(in)  :: spGath, spScat
      integer,       intent(in)    :: kind_C
      integer,       intent(in)    :: comm

				! Attributes: parameters for the matrix
      type(Navigator),intent(in) :: nav
      integer,dimension(:),intent(in) :: krNav		! (nvseg)
      integer,dimension(:),intent(in) :: ktNav		! (nvseg)
      integer,dimension(  :),intent(in) :: kx		! (  nvp)
      real   ,dimension(:,:),intent(in) :: qr		! (3,nvp)
      integer,dimension(  :),intent(in) :: kl		! (  nvp)

				! About the local transposed vectors
      type(Collector),intent(in) :: collVec
      real,dimension(:,:),intent(in ) :: xi	! (nvecs,nvi)
      real,dimension(:,:),intent(out) :: yi	! (nvecs,nvi)

      integer,optional,intent(out) :: stat

! !REVISION HISTORY:
! 	12Mar99 - Jing Guo <guo@thunder> - initial prototype/prolog/code
!	09Apr99 - J. Larson
!		- Fixed a bug in arguments call to MPI_reduce_scatter()
!	22Apr99 - Jing Guo
!		- Reordered the arguments to make them be consistent
!		  with recMatx_cxpy().
!       06Oct99 - J.W. Larson <jlarson@dao.gsfc.nasa.gov> - replaced 
!                 occurances of MP_REAL with MP_type(message) for 
!                 better support of 32 and 64-bit platforms.
!	29Oct99	- Jing Guo
!		. Added protection on out-of-range lbSeg (see lcPE).
!	10Nov99	- J.W.  Larson
!                 Added SHMEM_reduce_scatter call.
!EOP ___________________________________________________________________


  character(len=*),parameter :: myname_=myname//'::symCx_tr_'

  integer :: ier
  integer :: nvi,nvp,nvseg,nvecs
  integer :: lc,le

  real,allocatable,dimension(:,:) :: xp,yp

	! Assume arguments already verified

  nvi= localSize(collVec)
  nvp=globalSize(collVec)
  nvseg=lsize(nav)
  nvecs=size(xi,1)

  if(present(stat)) stat=0

	! Allocate the y' and x'.

  allocate(xp(nvecs,nvp),yp(nvecs,nvp),stat=ier)
	if(ier/=0) then
	  call perr(myname_,'allocate()',ier)
	  if(.not.present(stat)) call die(myname_)
	  stat=ier
	  return
	endif

	if(mall_ison()) then
	  call mall_mci(xp,myname)
	  call mall_mci(yp,myname)
	endif

#ifdef VERIFY_
  call checknorm(myname_//'(xi)',transpose(xi))
  call checknorm(myname_//'(xi)',transpose(xi),comm)
#endif
!________________________________________

	! At this point, a "local" vector x' and a "partial-product"
	! vector y' have been created, with their shared aNav.

	! The performance of this procedure relies on the assumption
	! that a log(nPEs) communication pattern is used.  If for a
	! given platform, the system does not support a log(nPEs)
	! communication routine, one may consider to implement an own
	! version.  This applies to both MPI_allGather() and
	! MPI_reduce_scatter()

	call zeit_ci("sMatxO_agat")

	Call PartialGather(xi, xp, spGath)

	call zeit_co("sMatxO_agat")

#ifdef VERIFY_
  call checknorm(myname_//'(xp)',transpose(xp))
  call checknorm(myname_//'(xp)',transpose(xp),comm)
#endif
!________________________________________

	! Initialize the full vector to zeroes.

  yp(:,:)=0.

	! A all-local matrix-vector multiplication

	call zeit_ci("sMatxOxpy")

  call symCxpy(mNav,kind_C,			&
	nav,nvseg,krNav,ktNav,nvp,kx,qr,kl,	&
	nvecs,xp,yp, stat=ier)
        
	if(ier/=0) then
	  call perr(myname_,'symCxpy()',ier)
	  if(.not.present(stat)) call die(myname_)
	  stat=ier
	  return
	endif

	call zeit_co("sMatxOxpy")

#ifdef VERIFY_
  call checknorm(myname_//'(yp)',transpose(yp))
  call checknorm(myname_//'(yp)',transpose(yp),comm)
#endif
!________________________________________

	call zeit_ci("sMatxO_rscat")

		! Should expand m_CollectorComm to allow this class
		! of operations (reduce_scatter()).  Later ..

	call PartialReduceScatter(yp, yi, spScat)

	call zeit_co("sMatxO_rscat")

#ifdef VERIFY_
  call checknorm(myname_//'(yi)',transpose(yi))
  call checknorm(myname_//'(yi)',transpose(yi),comm)
#endif
!________________________________________

	! At this point, one can remove the data objects created at the
	! beginning of the precedure.

	if(mall_ison()) then
	  call mall_mco(xp,myname)
	  call mall_mco(yp,myname)
	endif

  deallocate(xp,yp,stat=ier)
	if(ier/=0) then
	  call perr(myname_,'deallocate()',ier)
	  if(.not.present(stat)) call die(myname_)
	  stat=ier
	  return
	endif

end subroutine symCx_tr_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: recCx_ - an interface of message-passing rectangular Cx
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine recCx_(mNav,kind_C,comm,				&
	row_nav,row_krNav,row_ktNav,row_kx,row_qr,row_kl,	&
	col_nav,col_krNav,col_ktNav,col_kx,col_qr,col_kl,	&
	x_collVec,xi, y_collVec,yi, stat			)

      use m_Collector,only : Collector
      use m_Collector,only : localSize
      use m_Collector,only : globalSize

      use m_Navigator,only : Navigator
      use m_Navigator,only : lsize
      use m_Checks,only : checknorm

      use m_AttrVect ,only : AttrVect
      use m_die,only : perr,die,assert_
      use m_mall,only : mall_ison,mall_mci,mall_mco

      implicit none

				! About the matrix
      type(AttrVect),intent(inout) :: mNav
      integer,       intent(in)    :: kind_C
      integer,       intent(in)    :: comm

				! Attributes: parameters for the matrix
      type(Navigator),intent(in) :: row_nav
      integer,dimension(:),intent(in) :: row_krNav	! (nySeg)
      integer,dimension(:),intent(in) :: row_ktNav	! (nySeg)
      integer,dimension(  :),intent(in) :: row_kx	! (  nyp)
      real   ,dimension(:,:),intent(in) :: row_qr	! (3,nyp)
      integer,dimension(  :),intent(in) :: row_kl	! (  nyp)

      type(Navigator),intent(in) :: col_nav
      integer,dimension(:),intent(in) :: col_krNav	! (nxSeg)
      integer,dimension(:),intent(in) :: col_ktNav	! (nxSeg)
      integer,dimension(  :),intent(in) :: col_kx	! (  nxp)
      real   ,dimension(:,:),intent(in) :: col_qr	! (3,nxp)
      integer,dimension(  :),intent(in) :: col_kl	! (  nxp)

				! About the local vectors
      type(Collector),intent(in) :: x_collVec
      real,dimension(:,:),intent(in ) :: xi		! (nxi,nvecs)
      type(Collector),intent(in) :: y_collVec
      real,dimension(:,:),intent(out) :: yi		! (nyi,nvecs)

      integer,optional,intent(out) :: stat

! !REVISION HISTORY:
! 	23Aug00	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::recCx_'
  integer :: nxi,nyi
  integer :: nxp,nyp
  integer :: nxSeg,nySeg
  integer :: nvecs
  integer :: ier
  real,allocatable,dimension(:,:) :: ty

	! Argument consistency verifications

  nxi= localSize(x_collVec)
	ASSERT( nxi<=size(xi,1) )

  nyi= localSize(y_collVec)
	ASSERT( nyi==size(yi,1) )	! yi is an output array

  nxp=globalSize(x_collVec)
	ASSERT( nxp<=size(col_kx,1) )
	ASSERT( nxp<=size(col_qr,2) )
	ASSERT( nxp<=size(col_kl,1) )

  nyp=globalSize(y_collVec)
	ASSERT( nyp<=size(row_kx,1) )
	ASSERT( nyp<=size(row_qr,2) )
	ASSERT( nyp<=size(row_kl,1) )

  nxSeg=lsize(col_nav)
	ASSERT( nxSeg<=size(col_krNav) )
	ASSERT( nxSeg<=size(col_ktNav) )

  nySeg=lsize(row_nav)
	ASSERT( nxSeg<=size(row_krNav) )
	ASSERT( nxSeg<=size(row_ktNav) )

  nvecs=size(xi,2)
	ASSERT( nvecs==size(yi,2) )	! yi is an output array

  if(present(stat)) stat=0

	allocate(ty(nvecs,nyi),stat=ier)
		if(ier/=0) then
		  call perr(myname_,'allocate()',ier)
		  if(.not.present(stat)) call die(myname_)
		  stat=ier
		  return
		endif

		if(mall_ison()) call mall_mci(ty,myname)

#ifdef VERIFY_
  call checknorm(myname_//'(xi)',xi(1:nxi,1:nvecs))
  call checknorm(myname_//'(xi)',xi(1:nxi,1:nvecs),comm)
#endif

  call recCx_tr_(mNav,kind_C,comm,				&
	row_nav,row_krNav,row_ktNav,row_kx,row_qr,row_kl,	&
	col_nav,col_krNav,col_ktNav,col_kx,col_qr,col_kl,	&
	x_collVec,transpose(xi), y_collVec,ty, stat=ier		)

	if(ier/=0) then
	  call perr(myname_,'recCx_tr_()',ier)
	  if(.not.present(stat)) call die(myname_)
	  stat=ier
	  return
	endif

  yi(1:nyi,1:nvecs) = transpose(ty(:,:))

#ifdef VERIFY_
  call checknorm(myname_//'(yi)',yi(1:nyi,1:nvecs))
  call checknorm(myname_//'(yi)',yi(1:nyi,1:nvecs),comm)
#endif
		if(mall_ison()) call mall_mco(ty,myname)

	deallocate(ty,stat=ier)
		if(ier/=0) then
		  call perr(myname_,'deallocate()',ier)
		  if(.not.present(stat)) call die(myname_)
		  stat=ier
		  return
		endif

end subroutine recCx_
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: recCxpy_ - an interface of message-passing rectangular Cx+y
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine recCxpy_(mNav,kind_C,comm,			&
	row_nav,row_krNav,row_ktNav,row_kx,row_qr,row_kl,	&
	col_nav,col_krNav,col_ktNav,col_kx,col_qr,col_kl,	&
	x_collVec,xi, y_collVec,yi, stat			)

      use m_Collector,only : Collector
      use m_Collector,only : localSize
      use m_Collector,only : globalSize

      use m_Navigator,only : Navigator
      use m_Navigator,only : lsize
      use m_Checks,only : checknorm

      use m_AttrVect ,only : AttrVect
      use m_die,only : perr,die,assert_
      use m_mall,only : mall_ison,mall_mci,mall_mco

      implicit none

				! About the matrix
      type(AttrVect),intent(inout) :: mNav
      integer,       intent(in)    :: kind_C
      integer,       intent(in)    :: comm

				! Attributes: parameters for the matrix
      type(Navigator),intent(in) :: row_nav
      integer,dimension(:),intent(in) :: row_krNav	! (nySeg)
      integer,dimension(:),intent(in) :: row_ktNav	! (nySeg)
      integer,dimension(  :),intent(in) :: row_kx	! (  nyp)
      real   ,dimension(:,:),intent(in) :: row_qr	! (3,nyp)
      integer,dimension(  :),intent(in) :: row_kl	! (  nyp)

      type(Navigator),intent(in) :: col_nav
      integer,dimension(:),intent(in) :: col_krNav	! (nxSeg)
      integer,dimension(:),intent(in) :: col_ktNav	! (nxSeg)
      integer,dimension(  :),intent(in) :: col_kx	! (  nxp)
      real   ,dimension(:,:),intent(in) :: col_qr	! (3,nxp)
      integer,dimension(  :),intent(in) :: col_kl	! (  nxp)

				! About the local vectors
      type(Collector),intent(in) :: x_collVec
      real,dimension(:,:),intent(in ) :: xi		! (nxi,nvecs)
      type(Collector),intent(in) :: y_collVec
      real,dimension(:,:),intent(inout) :: yi		! (nyi,nvecs)

      integer,optional,intent(out) :: stat

! !REVISION HISTORY:
! 	23Aug00	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::recCxpy_'
  integer :: nxi,nyi
  integer :: nxp,nyp
  integer :: nxSeg,nySeg
  integer :: nvecs
  integer :: ier
  real,allocatable,dimension(:,:) :: ty

	! Argument consistency verifications

  nxi= localSize(x_collVec)
	ASSERT( nxi<=size(xi,1) )

  nyi= localSize(y_collVec)
	ASSERT( nyi==size(yi,1) )	! yi is an output array

  nxp=globalSize(x_collVec)
	ASSERT( nxp<=size(col_kx,1) )
	ASSERT( nxp<=size(col_qr,2) )
	ASSERT( nxp<=size(col_kl,1) )

  nyp=globalSize(y_collVec)
	ASSERT( nyp<=size(row_kx,1) )
	ASSERT( nyp<=size(row_qr,2) )
	ASSERT( nyp<=size(row_kl,1) )

  nxSeg=lsize(col_nav)
	ASSERT( nxSeg<=size(col_krNav) )
	ASSERT( nxSeg<=size(col_ktNav) )

  nySeg=lsize(row_nav)
	ASSERT( nxSeg<=size(row_krNav) )
	ASSERT( nxSeg<=size(row_ktNav) )

  nvecs=size(xi,2)
	ASSERT( nvecs==size(yi,2) )	! yi is an output array

  if(present(stat)) stat=0

	allocate(ty(nvecs,nyi),stat=ier)
		if(ier/=0) then
		  call perr(myname_,'allocate()',ier)
		  if(.not.present(stat)) call die(myname_)
		  stat=ier
		  return
		endif

		if(mall_ison()) call mall_mci(ty,myname)
#ifdef VERIFY_
  call checknorm(myname_//'(xi)',xi(1:nxi,1:nvecs))
  call checknorm(myname_//'(xi)',xi(1:nxi,1:nvecs),comm)
  call checknorm(myname_//'(yi)',yi(1:nyi,1:nvecs))
  call checknorm(myname_//'(yi)',yi(1:nyi,1:nvecs),comm)
#endif

  call recCx_tr_(mNav,kind_C,comm,				&
	row_nav,row_krNav,row_ktNav,row_kx,row_qr,row_kl,	&
	col_nav,col_krNav,col_ktNav,col_kx,col_qr,col_kl,	&
	x_collVec,transpose(xi), y_collVec,ty, stat=ier		)

	if(ier/=0) then
	  call perr(myname_,'recCx_tr_()',ier)
	  if(.not.present(stat)) call die(myname_)
	  stat=ier
	  return
	endif

  yi(1:nyi,1:nvecs) = transpose(ty) + yi(1:nyi,1:nvecs)

#ifdef VERIFY_
  call checknorm(myname_//'(yi)',yi(1:nyi,1:nvecs))
  call checknorm(myname_//'(yi)',yi(1:nyi,1:nvecs),comm)
#endif
		if(mall_ison()) call mall_mco(ty,myname)

	deallocate(ty,stat=ier)
		if(ier/=0) then
		  call perr(myname_,'deallocate()',ier)
		  if(.not.present(stat)) call die(myname_)
		  stat=ier
		  return
		endif

end subroutine recCxpy_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !ROUTINE: recCx_tr_ - a message-passing rectangular Cx kernal
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine recCx_tr_(mNav,kind_C,comm,			&
	row_nav,row_krNav,row_ktNav,row_kx,row_qr,row_kl,	&
	col_nav,col_krNav,col_ktNav,col_kx,col_qr,col_kl,	&
	x_collVec,xi, y_collVec,yi, stat			)

      use m_Collector, only : Collector
      use m_Collector, only : globalSize
      use m_Collector, only : localSize
      use m_Collector, only : ptr_counts
      use m_CollectorComm,only : allgatherv

      use m_Navigator, only : Navigator
      use m_Navigator, only : lsize
      use m_Checks,only : checknorm

      use m_AttrVect,  only : AttrVect

      use m_CorOxpy,   only : recCxpy

      use m_die   ,only : perr,die,assert_
      use m_mall  ,only : mall_ison,mall_mci,mall_mco
      use m_mpif90,only : MP_perr,MP_type,MP_SUM
      use m_zeit  ,only : zeit_ci,zeit_co

      implicit none
				! About the matrix
      type(AttrVect),intent(inout) :: mNav
      integer,       intent(in)    :: kind_C
      integer,       intent(in)    :: comm

				! Attributes: parameters for the matrix
      type(Navigator),intent(in) :: row_nav
      integer,dimension(:),intent(in) :: row_krNav	! (nySeg)
      integer,dimension(:),intent(in) :: row_ktNav	! (nySeg)
      integer,dimension(  :),intent(in) :: row_kx	! (  nyp)
      real   ,dimension(:,:),intent(in) :: row_qr	! (3,nyp)
      integer,dimension(  :),intent(in) :: row_kl	! (  nyp)

      type(Navigator),intent(in) :: col_nav
      integer,dimension(:),intent(in) :: col_krNav	! (nxSeg)
      integer,dimension(:),intent(in) :: col_ktNav	! (nxSeg)
      integer,dimension(  :),intent(in) :: col_kx	! (  nxp)
      real   ,dimension(:,:),intent(in) :: col_qr	! (3,nxp)
      integer,dimension(  :),intent(in) :: col_kl	! (  nxp)

				! About the local vectors
      type(Collector),intent(in) :: x_collVec
      real,dimension(:,:),intent(in   ) :: xi		! (nvecs,nxi)
      type(Collector),intent(in) :: y_collVec
      real,dimension(:,:),intent(out) :: yi		! (nvecs,nyi)

      integer,optional,intent(out) :: stat

! !REVISION HISTORY:
! 	12Mar99 - Jing Guo <guo@thunder> - initial prototype/prolog/code
!	09Apr99 - J. Larson
!		- Fixed a bug in arguments call to MPI_reduce_scatter()
!       06Oct99 - J.W. Larson <jlarson@dao.gsfc.nasa.gov> - replaced 
!                 occurances of MP_REAL with MP_type(message) for 
!                 better support of 32 and 64-bit platforms.
!	29Oct99	- Jing Guo
!		. Added protections on out-of-range lbSeg (see lcPE).
!EOP ___________________________________________________________________


  character(len=*),parameter :: myname_=myname//'::recCx_tr_'

  integer :: ier
  integer :: nxi,nxp,nxSeg,nvecs
  integer :: nyi,nyp,nySeg
  integer :: ln,lc,le

  real,dimension(:,:),allocatable :: xp,yp
!________________________________________

	! Assume all arguments have been verified

  nxi= localSize(x_collVec)
  nyi= localSize(y_collVec)

  nxp=globalSize(x_collVec)
  nyp=globalSize(y_collVec)

  nxSeg=lsize(col_nav)
  nySeg=lsize(row_nav)

  nvecs=size(xi,1)

	! Allocate the y' and x'.

  allocate(xp(nvecs,nxp),yp(nvecs,nyp),stat=ier)
	if(ier/=0) then
	  call perr(myname_,'allocate()',ier)
	  if(.not.present(stat)) call die(myname_)
	  stat=ier
	  return
	endif

	if(mall_ison()) then
	   call mall_mci(xp,myname)
	   call mall_mci(yp,myname)
	endif

#ifdef VERIFY_
  call checknorm(myname_//'(xi)',transpose(xi))
  call checknorm(myname_//'(xi)',transpose(xi),comm)
#endif
!________________________________________

	! At this point, a "local" vector x' and a "partial-product"
	! vector y' have been created, with their shared aNav.

	! The performance of this procedure relies on the assumption
	! that a log(nPE) communication pattern is used.  If for a
	! given platform, the system does not support a log(nPE)
	! communication routine, one may consider to implement an own
	! version.  This applies to both MPI_allGather() and
	! MPI_reduce_scatter()

	call zeit_ci("rMatxO_agat")

  call allgatherv(xi,xp,x_collvec,comm,stat=ier)
	if(ier/=0) then
	  call perr(myname_,'allgatherv()',ier)
	  if(.not.present(stat)) call die(myname_)
	  stat=ier
	  return
	endif

	call zeit_co("rMatxO_agat")

#ifdef VERIFY_
  call checknorm(myname_//'(xp)',transpose(xp))
  call checknorm(myname_//'(xp)',transpose(xp),comm)
#endif
!________________________________________

	! Initialize the full vector to all zeroes.

  yp(:,:)=0.

	! The local matrix-vector multiplication

	call zeit_ci("rMatxOxpy")

  call recCxpy(mNav,kind_C,			&
	row_nav,nySeg,row_krNav,row_ktNav,	&
		nyp,row_kx,row_qr,row_kl,	&
	col_nav,nxSeg,col_krNav,col_ktNav,	&
		nxp,col_kx,col_qr,col_kl,	&
	nvecs,xp,yp, stat=ier)
        
	if(ier/=0) then
	  call perr(myname_,'recCxpy()',ier)
	  if(.not.present(stat)) call die(myname_)
	  stat=ier
	  return
	endif

	call zeit_co("rMatxOxpy")

#ifdef VERIFY_
  call checknorm(myname_//'(yp)',transpose(yp))
  call checknorm(myname_//'(yp)',transpose(yp),comm)
#endif
!________________________________________

	call zeit_ci("rMatxO_rscat")


  call MPI_reduce_scatter(yp,yi,nvecs*ptr_counts(y_collVec),	&
	MP_type(yi(1,1)),MP_SUM,comm,ier)

	if(ier/=0) then
	  call MP_perr(myname_,'MPI_reduce_scatter()',ier)
	  if(.not.present(stat)) call die(myname_)
	  stat=ier
	  return
	endif

	call zeit_co("rMatxO_rscat")

#ifdef VERIFY_
  call checknorm(myname_//'(yi)',transpose(yi))
  call checknorm(myname_//'(yi)',transpose(yi),comm)
#endif
!________________________________________

	! At this point, one can remove the data objects created at the
	! beginning of the precedure.

	if(mall_ison()) then
	  call mall_mco(xp,myname)
	  call mall_mco(yp,myname)
	endif

  deallocate(xp,yp,stat=ier)
	if(ier/=0) then
	  call perr(myname_,'deallocate()',ier)
	  if(.not.present(stat)) call die(myname_)
	  stat=ier
	  return
	endif

end subroutine recCx_tr_

end module m_CorMatxO
!.
